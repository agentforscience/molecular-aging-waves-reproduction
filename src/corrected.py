"""Corrected inference for age-localised molecular change.

Carbonneau et al. recommend but do not supply a replacement for LOESS+DE-SWAN.
This module implements one, built from well-understood estimators applied
correctly rather than from a new method that would itself need validating.

The design is shared by every molecule (all molecules are measured on the same
subjects, hence the same age vector), so every model here is fitted for all
molecules simultaneously by a single least-squares solve against a precomputed
pseudo-inverse. That also makes the permutation calibration cheap: permuting age
labels is, in sorted-age coordinates, a column permutation of the data matrix,
so the design and all its pseudo-inverses are reused unchanged.

Three questions, deliberately kept separate because the literature conflates
them:

  Q1  Is there ANY age association?         -> linear F-test
  Q2  Is the age association NONLINEAR?      -> natural cubic spline vs linear,
                                                nested F-test
  Q3  Is there a DISCRETE TRANSITION, and at what age? -> segmented (broken-stick)
                                                regression, sup-F over candidate
                                                breakpoints, permutation-calibrated
                                                because the breakpoint is a
                                                nuisance parameter unidentified
                                                under the null (Davies 1987)

All p-values are calibrated against a permutation null built by shuffling ages
*before* any modelling, and controlled with Benjamini-Hochberg (independence /
PRDS) and Benjamini-Yekutieli (arbitrary dependence, appropriate here because
molecules are correlated).
"""
from __future__ import annotations

import numpy as np


# ------------------------------------------------------------------ utilities
def _rss_all(Y: np.ndarray, X: np.ndarray, ss_y: np.ndarray | None = None) -> np.ndarray:
    """Residual sum of squares of every row of Y regressed on X.

    Y : (p, n), X : (n, k), k << n. Returns (p,).

    Uses an orthonormal basis Q of col(X) from a QR decomposition, so that
    RSS = ||y||^2 - ||Q^T y||^2. This costs one (p, n) x (n, k) matmul rather
    than forming and applying the (n, n) hat matrix -- an O(n/k) saving that is
    what makes the n=656 methylation cohort and its permutation nulls tractable.
    Rank deficiency is handled by dropping the null-space columns of R.

    `ss_y` (row sums of squares of Y) may be passed in to avoid recomputing it
    across many candidate designs on the same Y.
    """
    Q, R = np.linalg.qr(X)
    rank = int(np.sum(np.abs(np.diag(R)) > 1e-10 * max(1.0, np.abs(np.diag(R)).max())))
    Q = Q[:, :rank]
    if ss_y is None:
        ss_y = np.einsum("ij,ij->i", Y, Y)
    proj = Y @ Q                                    # (p, rank)
    return np.maximum(ss_y - np.einsum("ij,ij->i", proj, proj), 0.0)


def f_test(rss0: np.ndarray, rss1: np.ndarray, df0: int, df1: int, n: int):
    """Nested-model F statistic and its (nominal) p-value."""
    from scipy.stats import f as fdist
    num = (rss0 - rss1) / max(df1 - df0, 1)
    den = rss1 / max(n - df1, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        F = np.where(den > 0, num / den, np.nan)
    F = np.maximum(F, 0.0)
    return F, fdist.sf(F, df1 - df0, n - df1)


def natural_spline_basis(x: np.ndarray, df: int = 4) -> np.ndarray:
    """Natural cubic spline basis with `df` degrees of freedom (no intercept).

    Knots at equally spaced quantiles of x, boundary knots at the extremes --
    the standard `splines::ns` construction. Implemented directly to avoid a
    dependency and to keep the basis identical across permutations.
    """
    x = np.asarray(x, float)
    n_interior = max(df - 1, 0)
    qs = np.linspace(0, 1, n_interior + 2)[1:-1]
    knots = np.quantile(x, qs) if n_interior else np.array([])
    bk = np.array([x.min(), x.max()])
    all_k = np.concatenate([[bk[0]], knots, [bk[1]]])

    def d(k, j):
        num = np.clip(x - all_k[k], 0, None) ** 3 - np.clip(x - all_k[j], 0, None) ** 3
        return num / (all_k[j] - all_k[k])

    K = all_k.size
    cols = [x]
    for k in range(K - 2):
        cols.append(d(k, K - 1) - d(K - 2, K - 1))
    B = np.column_stack(cols)
    # Standardise columns for numerical conditioning; span is unchanged.
    B = (B - B.mean(0)) / np.where(B.std(0) > 0, B.std(0), 1.0)
    return B


# --------------------------------------------------------------- Q1 / Q2 tests
def linear_and_spline_tests(Y: np.ndarray, ages: np.ndarray, df: int = 4):
    """Q1 and Q2 for every molecule.

    Returns dict with F/p for (a) linear vs intercept-only and (b) spline vs
    linear, plus R^2 of each model.
    """
    n = ages.size
    X0 = np.ones((n, 1))
    X1 = np.column_stack([np.ones(n), ages])
    X2 = np.column_stack([np.ones(n), natural_spline_basis(ages, df)])

    ss_y = np.einsum("ij,ij->i", Y, Y)
    rss0 = _rss_all(Y, X0, ss_y)
    rss1 = _rss_all(Y, X1, ss_y)
    rss2 = _rss_all(Y, X2, ss_y)

    F_lin, p_lin = f_test(rss0, rss1, 1, 2, n)
    F_nl, p_nl = f_test(rss1, rss2, 2, X2.shape[1], n)
    with np.errstate(divide="ignore", invalid="ignore"):
        r2_lin = 1 - rss1 / rss0
        r2_spl = 1 - rss2 / rss0
    return dict(F_linear=F_lin, p_linear=p_lin, F_nonlinear=F_nl, p_nonlinear=p_nl,
                r2_linear=r2_lin, r2_spline=r2_spl, df_spline=X2.shape[1])


# ------------------------------------------------------------------- Q3 sup-F
def segmented_supF(Y: np.ndarray, ages: np.ndarray, psi_grid: np.ndarray,
                   min_per_side: int = 10):
    """Broken-stick regression: sup-F over candidate breakpoints.

    For each candidate breakpoint psi the model is
        y = b0 + b1*age + b2*(age - psi)_+ + e
    i.e. a continuous piecewise-linear fit with a slope change at psi. The test
    statistic is the maximum F for the added term over the grid; the maximiser is
    the breakpoint estimate.

    The distribution of sup-F is NOT F: the breakpoint is unidentified under the
    null (Davies 1987), so nominal p-values are anticonservative. Calibrate with
    `permutation_null_supF`.

    Returns (supF (p,), psi_hat (p,), F_all (p, n_psi)).
    """
    n = ages.size
    X1 = np.column_stack([np.ones(n), ages])
    ss_y = np.einsum("ij,ij->i", Y, Y)
    rss1 = _rss_all(Y, X1, ss_y)

    keep = [i for i, psi in enumerate(psi_grid)
            if (ages < psi).sum() >= min_per_side and (ages >= psi).sum() >= min_per_side]
    F_all = np.full((Y.shape[0], psi_grid.size), np.nan)
    for i in keep:
        psi = psi_grid[i]
        X2 = np.column_stack([np.ones(n), ages, np.clip(ages - psi, 0, None)])
        rss2 = _rss_all(Y, X2, ss_y)
        F_all[:, i], _ = f_test(rss1, rss2, 2, 3, n)

    supF = np.nanmax(F_all, axis=1)
    psi_hat = np.where(np.isfinite(supF),
                       psi_grid[np.nanargmax(np.nan_to_num(F_all, nan=-1), axis=1)],
                       np.nan)
    return supF, psi_hat, F_all


def permutation_null_supF(Y: np.ndarray, ages: np.ndarray, psi_grid: np.ndarray,
                          n_perm: int, rng, min_per_side: int = 10):
    """Permutation null for the three tests, sharing one set of permutations.

    Returns dict of (n_perm, p) arrays for supF / F_linear / F_nonlinear, plus
    the (n_perm, p) matrix of null breakpoint locations. Permuting age labels is
    equivalent to permuting the columns of Y, which leaves every design matrix
    and pseudo-inverse untouched -- the reason this is affordable.
    """
    out = dict(supF=[], psi=[], F_linear=[], F_nonlinear=[])
    for _ in range(n_perm):
        Yp = Y[:, rng.permutation(Y.shape[1])]
        t = linear_and_spline_tests(Yp, ages)
        s, psi, _ = segmented_supF(Yp, ages, psi_grid, min_per_side)
        out["supF"].append(s)
        out["psi"].append(psi)
        out["F_linear"].append(t["F_linear"])
        out["F_nonlinear"].append(t["F_nonlinear"])
    return {k: np.asarray(v) for k, v in out.items()}


def perm_pvalues(obs: np.ndarray, null: np.ndarray) -> np.ndarray:
    """Per-molecule permutation p-values using the POOLED null across molecules.

    Pooling is valid because all molecules share the same design and, under the
    null of no age association, the same exchangeability structure; it buys
    n_perm * p effective null draws instead of n_perm, which is what makes
    resolution below 1/n_perm possible. Ties are handled by the usual +1 rule.

    Non-finite observed statistics (a variable with missing values yields a NaN
    RSS and hence a NaN F) return NaN, NOT a small p-value. `np.searchsorted`
    places NaN at the end of the sorted null, which would otherwise hand every
    unevaluable variable the most extreme possible p-value -- a bug that
    initially manifested as 12/51 spurious "transitions" in the one iPOP layer
    that has missing data.
    """
    obs = np.asarray(obs, float)
    flat = np.sort(null.ravel()[np.isfinite(null.ravel())])
    out = np.full(obs.shape, np.nan)
    ok = np.isfinite(obs)
    if flat.size and ok.any():
        idx = np.searchsorted(flat, obs[ok], side="left")
        out[ok] = (1.0 + (flat.size - idx)) / (1.0 + flat.size)
    return out
