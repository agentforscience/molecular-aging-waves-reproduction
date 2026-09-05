"""Vectorised, R-compatible LOESS for the iPOP aging pipeline.

Why this module exists
----------------------
`code/loess_py/rloess.py` is a faithful but pure-Python reimplementation of R's
`stats::loess(family="gaussian", degree=2, surface="direct")` plus the
leave-one-out span search used by jaspershen-lab/ipop_aging. It fits every
variable independently, at ~1 min per 300 variables. The 8,556-gene iPOP
transcriptome takes ~25 min, which makes any permutation study (hundreds of
re-smoothings) impossible.

The key observation that makes this tractable
---------------------------------------------
A local polynomial fit is *linear in y*. The prediction at a query point x0 is

    yhat(x0) = e_1^T (X^T W X)^{-1} X^T W y  ==  l(x0)^T y

where X is the local Vandermonde design centred at x0 and W the tricube weight
matrix. Neither X nor W depends on y. So for a fixed design (the vector of
subject ages) and a fixed span, the entire smoother is a single matrix

    L  of shape (n_grid, n_samples),        Yhat = Y @ L.T

which we build once and reuse for every one of the p variables. The same trick
applies to the leave-one-out CV used for span selection: the LOO prediction at
sample k is also a fixed linear functional of y, giving a matrix C of shape
(n_loo, n_samples) with a structural zero in column k.

A second observation makes permutation testing nearly free
----------------------------------------------------------
Permuting the age labels among subjects leaves the *set* of design points
unchanged. In sorted-age coordinates the operator matrices L and C are therefore
identical across permutations; only the assignment of y-values to sorted
positions changes. A permutation is thus a column shuffle of Y followed by the
same matmul -- no re-derivation of any smoother. This is what turns the ~40-hour
permutation study estimated in STATE.md into minutes.

Numerical agreement with the reference implementation is asserted in
`validate_fastloess.py` (exact to floating point against `code/loess_py`, and
compared against the authors' own published R output).
"""
from __future__ import annotations

import numpy as np

DEFAULT_SPANS = (0.3, 0.4, 0.5, 0.6)
DEFAULT_GRID = np.arange(26.0, 75.5, 0.5)  # Shen et al.: half-year grid, 26-75 => 99 points


def _tricube(u: np.ndarray) -> np.ndarray:
    """Tricube kernel (1 - |u|^3)^3, clipped at 0. Matches R's loess weighting."""
    w = np.clip(1.0 - np.abs(u) ** 3, 0.0, None)
    return w ** 3


def _local_weight_row(x: np.ndarray, x0: float, span: float, degree: int) -> np.ndarray:
    """Row l(x0) of the smoother operator: yhat(x0) = l(x0) . y.

    Reproduces `code/loess_py/rloess.py:loess_fit` for a single query point:
    take the q = floor(n*span) nearest design points, set the bandwidth to the
    largest of those distances (inflated by `span` when span > 1, as R does),
    weight by tricube, and read off the intercept of a weighted degree-`degree`
    polynomial fit centred at x0.

    Returns a length-n row of weights, or a row of NaN if the neighbourhood is
    degenerate (fewer than degree+1 points with positive weight).
    """
    n = x.size
    q = int(np.floor(n * span)) if span <= 1 else n
    q = max(q, degree + 1)

    d = np.abs(x - x0)
    idx = np.argpartition(d, q - 1)[:q] if q < n else np.arange(n)
    dq = d[idx]
    h = dq.max()
    if span > 1:
        h *= span

    w = _tricube(dq / h) if h > 0 else np.ones_like(dq)
    keep = w > 0
    row = np.zeros(n)
    if keep.sum() < degree + 1:
        row[:] = np.nan
        return row

    xi, wi = x[idx][keep], w[keep]
    X = np.vander(xi - x0, degree + 1, increasing=True)
    sw = np.sqrt(wi)
    Xw = X * sw[:, None]
    # Least-squares solution operator: beta = pinv(Xw) @ (y*sw); we need beta[0].
    # Row 0 of pinv(Xw), rescaled by sw, is the linear functional we want.
    pinv0 = np.linalg.pinv(Xw)[0]
    row[idx[keep]] = pinv0 * sw
    return row


def loess_operator(x: np.ndarray, xout: np.ndarray, span: float, degree: int = 2) -> np.ndarray:
    """Smoother matrix L with L @ y == LOESS fit of y on x evaluated at xout.

    Parameters
    ----------
    x : (n,) design points (subject ages). Need not be sorted.
    xout : (m,) query points (the prediction grid).
    span, degree : LOESS parameters, as in R.

    Returns
    -------
    (m, n) ndarray.
    """
    x = np.asarray(x, float)
    xout = np.atleast_1d(np.asarray(xout, float))
    return np.vstack([_local_weight_row(x, x0, span, degree) for x0 in xout])


def loo_operator(x: np.ndarray, span: float, degree: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Leave-one-out prediction operator for span selection.

    Mirrors `optimize_loess_span` in ipop_aging / artifactual-waves-of-aging,
    which leaves out interior points only (R's `2:(nrow-1)` on sorted data).

    Returns
    -------
    C : (n-2, n) ndarray with C[k, held_out[k]] == 0 by construction.
    held_out : (n-2,) indices (into sorted-x order) of the held-out samples.
    """
    x = np.asarray(x, float)
    n = x.size
    held_out = np.arange(1, n - 1)
    C = np.zeros((held_out.size, n))
    for k, idx in enumerate(held_out):
        mask = np.ones(n, bool)
        mask[idx] = False
        sub_row = _local_weight_row(x[mask], x[idx], span, degree)
        C[k, mask] = sub_row
    return C, held_out


class LoessEngine:
    """Precomputed LOESS machinery for one fixed design (one age vector).

    Build once, then smooth arbitrarily many variable matrices -- including
    permuted ones -- with matrix multiplications only.

    Notes
    -----
    All operators are expressed in *sorted-age* coordinates. `Y` passed to the
    methods below is in the caller's original sample order; the engine reorders
    internally via `self.order`.
    """

    def __init__(self, ages, grid=None, spans=DEFAULT_SPANS, degree: int = 2):
        ages = np.asarray(ages, float)
        self.order = np.argsort(ages, kind="stable")
        self.x = ages[self.order]
        self.n = self.x.size
        self.grid = np.asarray(DEFAULT_GRID if grid is None else grid, float)
        self.spans = tuple(spans)
        self.degree = degree

        # Grid smoother and LOO operator, one per candidate span.
        self.L = {s: loess_operator(self.x, self.grid, s, degree) for s in self.spans}
        self.C, self.held = {}, None
        for s in self.spans:
            C, held = loo_operator(self.x, s, degree)
            self.C[s] = C
            self.held = held

    # ---------------------------------------------------------------- helpers
    def _sorted(self, Y: np.ndarray) -> np.ndarray:
        """Reorder a (p, n) variable-by-sample matrix into sorted-age order."""
        Y = np.asarray(Y, float)
        if Y.ndim == 1:
            Y = Y[None, :]
        return Y[:, self.order]

    def select_spans(self, Y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Per-variable LOO-RMSE span selection, vectorised over variables.

        Returns
        -------
        best : (p,) chosen span per variable
        rmse : (p, n_spans) LOO RMSE per variable per candidate span
        """
        Ys = self._sorted(Y)
        truth = Ys[:, self.held]
        rmse = np.empty((Ys.shape[0], len(self.spans)))
        for j, s in enumerate(self.spans):
            pred = Ys @ self.C[s].T                      # (p, n_loo)
            resid = truth - pred
            with np.errstate(invalid="ignore"):
                rmse[:, j] = np.sqrt(np.nanmean(resid ** 2, axis=1))
        # np.argmin ties break to the first (smallest) span, matching the
        # reference loop's strict `<` update rule.
        best_idx = np.nanargmin(np.where(np.isfinite(rmse), rmse, np.inf), axis=1)
        best = np.asarray(self.spans, float)[best_idx]
        return best, rmse

    def smooth(self, Y: np.ndarray, spans=None) -> np.ndarray:
        """LOESS-smooth every row of Y onto `self.grid`.

        Parameters
        ----------
        Y : (p, n) variables x samples, in the caller's original sample order.
        spans : optional (p,) per-variable spans. If None, they are selected by
            leave-one-out CV -- i.e. the published procedure.

        Returns
        -------
        (p, len(grid)) ndarray of smoothed values.
        """
        Ys = self._sorted(Y)
        if spans is None:
            spans, _ = self.select_spans(Y)
        spans = np.asarray(spans, float)
        out = np.empty((Ys.shape[0], self.grid.size))
        for s in self.spans:
            sel = spans == s
            if sel.any():
                out[sel] = Ys[sel] @ self.L[s].T
        return out

    def smooth_permuted(self, Y: np.ndarray, rng: np.random.Generator,
                        reselect_spans: bool = True) -> np.ndarray:
        """Shuffle age labels *before* smoothing, then smooth (Carbonneau Alg. 1).

        Because the design points are unchanged by a relabelling, this is a
        column permutation of Y in sorted-age space followed by the usual
        matmul. `reselect_spans=True` re-runs CV on the permuted data, which is
        what re-running the full published pipeline on shuffled ages would do.
        """
        Ys = self._sorted(Y)
        perm = rng.permutation(self.n)
        Yp = Ys[:, perm]
        if reselect_spans:
            truth = Yp[:, self.held]
            rmse = np.empty((Yp.shape[0], len(self.spans)))
            for j, s in enumerate(self.spans):
                resid = truth - Yp @ self.C[s].T
                with np.errstate(invalid="ignore"):
                    rmse[:, j] = np.sqrt(np.nanmean(resid ** 2, axis=1))
            best_idx = np.nanargmin(np.where(np.isfinite(rmse), rmse, np.inf), axis=1)
            spans = np.asarray(self.spans, float)[best_idx]
        else:
            spans = np.full(Yp.shape[0], self.spans[0], float)
        out = np.empty((Yp.shape[0], self.grid.size))
        for s in self.spans:
            sel = spans == s
            if sel.any():
                out[sel] = Yp[sel] @ self.L[s].T
        return out
