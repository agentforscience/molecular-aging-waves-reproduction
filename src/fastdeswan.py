"""Vectorised DE-SWAN, plus the power-equalised variant proposed in this work.

Three window-defining schemes are provided so that the *power confound* named as
Artifact 3 by Carbonneau et al. can be isolated:

  `fixed_width`  -- the published scheme. A window of total width 2*half_width
                    centred on `mid`; young = [mid-h, mid), old = [mid, mid+h).
                    Group sizes vary along the age axis with the age density, so
                    statistical power varies along the age axis. This is the
                    scheme whose curve shape is under suspicion.
  `equal_n`      -- our correction. Take the k nearest subjects *below* mid and
                    the k nearest *above* it, so (n_young, n_old) == (k, k) at
                    every window centre by construction. Power due to sample size
                    is then constant along the age axis; window *width* varies
                    instead, and is reported.
  `quantile`     -- the DEswan R package default (window centres at age
                    quantiles). Included because Carbonneau et al. note it does
                    *not* fix the problem, since window length stays fixed.

Two test statistics are provided:
  `wilcoxon` -- Shen et al. 2024's "modified DE-SWAN" (Mann-Whitney U).
  `linear`   -- Lehallier et al. 2019's original formulation (OLS on a
                below/above indicator, equivalent to a two-sample t-test here).

Every function returns per-window group sizes alongside the p-values, because
those sizes are the mechanism under test and no published DE-SWAN curve reports
them.
"""
from __future__ import annotations

import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

DEFAULT_MIDPOINTS = np.arange(40.0, 65.0, 1.0)  # Shen et al. plot 40-65 in 1-yr steps


# --------------------------------------------------------------------- helpers
def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values, NaN-safe, applied along a 1-D array."""
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    if ok.sum():
        q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q


def by(p: np.ndarray) -> np.ndarray:
    """Benjamini-Yekutieli q-values (valid under arbitrary dependence)."""
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    if ok.sum():
        q[ok] = multipletests(p[ok], method="fdr_by")[1]
    return q


def mannwhitney_p(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Two-sided Mann-Whitney U p-values for every row of A vs the same row of B.

    A : (p, n1), B : (p, n2). Uses the asymptotic approximation with continuity
    and tie correction -- the same normal approximation R's `wilcox.test` falls
    back to when ties are present, which they always are on LOESS grid output at
    these group sizes. Vectorised over rows.
    """
    A = np.asarray(A, float)
    B = np.asarray(B, float)
    if A.shape[1] == 0 or B.shape[1] == 0:
        return np.full(A.shape[0], np.nan)
    res = stats.mannwhitneyu(A, B, alternative="two-sided", axis=1, method="asymptotic")
    return np.asarray(res.pvalue, float)


def _avg_ranks_and_tiesum(C: np.ndarray):
    """Row-wise average (mid-)ranks and the tie-correction sum, from one argsort.

    `scipy.stats.rankdata(..., axis=1)` is the readable way to do this but is the
    dominant cost inside permutation loops. The trick here: after sorting each
    row, every tie group occupies a contiguous block, so its average rank is
    simply (first_index + last_index)/2 + 1. Both endpoints are obtainable with
    running max / running min over the group-boundary indicators, which is a
    handful of vectorised passes rather than a per-row Python loop.

    Returns (ranks, tie_sum) where tie_sum[i] = sum over tie groups of (t^3 - t),
    the quantity entering the tie-corrected variance of U. Note that a group of
    size t contributes (t^3 - t)/t = t^2 - 1 per member, so the sum over groups
    equals a sum over positions -- which keeps it vectorised.
    """
    m = C.shape[1]
    order = np.argsort(C, axis=1, kind="stable")
    Cs = np.take_along_axis(C, order, axis=1)

    pos = np.arange(m)
    start = np.empty(Cs.shape, bool)
    start[:, 0] = True
    start[:, 1:] = Cs[:, 1:] != Cs[:, :-1]
    end = np.empty(Cs.shape, bool)
    end[:, -1] = True
    end[:, :-1] = Cs[:, :-1] != Cs[:, 1:]

    first = np.maximum.accumulate(np.where(start, pos, 0), axis=1)
    last = np.minimum.accumulate(np.where(end, pos, m - 1)[:, ::-1], axis=1)[:, ::-1]

    avg_sorted = (first + last) / 2.0 + 1.0
    ranks = np.empty_like(avg_sorted)
    np.put_along_axis(ranks, order, avg_sorted, axis=1)

    t = (last - first + 1).astype(float)          # tie-group size at each position
    tie_sum = (t * t - 1.0).sum(axis=1)
    return ranks, tie_sum


def mannwhitney_p_fast(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Mann-Whitney U p-values, row-wise -- a faster route to the same answer.

    Numerically identical to `mannwhitney_p` (two-sided asymptotic, continuity
    correction, tie-corrected variance); it only replaces `scipy.stats.rankdata`
    with the argsort-based mid-rank computation above, which is ~3-4x cheaper and
    is what makes the n=656 permutation nulls affordable.

    Ties are handled exactly, which matters: methylation beta values stored as
    float32 do contain exact ties, and an earlier tie-ignoring version was caught
    by `validate_fast_mwu` disagreeing with the reference by up to 0.9 in
    p-value. Every experiment that uses this path validates it against the
    reference on its own data before the permutation loop starts.
    """
    A = np.asarray(A, float)
    B = np.asarray(B, float)
    n1, n2 = A.shape[1], B.shape[1]
    if n1 == 0 or n2 == 0:
        return np.full(A.shape[0], np.nan)
    n = n1 + n2
    C = np.concatenate([A, B], axis=1)
    ranks, tie_sum = _avg_ranks_and_tiesum(C)

    U1 = ranks[:, :n1].sum(axis=1) - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    var = (n1 * n2 / 12.0) * ((n + 1.0) - tie_sum / (n * (n - 1.0)))
    sd = np.sqrt(np.maximum(var, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(sd > 0, (np.abs(U1 - mu) - 0.5) / sd, 0.0)
    return 2.0 * stats.norm.sf(np.maximum(z, 0.0))


def validate_fast_mwu(A: np.ndarray, B: np.ndarray, tol: float = 1e-9) -> float:
    """Max absolute discrepancy between the fast and reference MWU p-values."""
    d = float(np.nanmax(np.abs(mannwhitney_p_fast(A, B) - mannwhitney_p(A, B))))
    return d


def welch_p(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Two-sided Welch t-test p-values, row-wise. The `linear` DE-SWAN variant."""
    A, B = np.asarray(A, float), np.asarray(B, float)
    if A.shape[1] < 2 or B.shape[1] < 2:
        return np.full(A.shape[0], np.nan)
    return np.asarray(stats.ttest_ind(A, B, axis=1, equal_var=False).pvalue, float)


def cohens_d(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Row-wise Cohen's d (pooled SD). Effect size with a scale, per Carbonneau
    et al.'s recommendation to prefer effect sizes over counts of significance."""
    A, B = np.asarray(A, float), np.asarray(B, float)
    n1, n2 = A.shape[1], B.shape[1]
    if n1 < 2 or n2 < 2:
        return np.full(A.shape[0], np.nan)
    s1, s2 = A.var(axis=1, ddof=1), B.var(axis=1, ddof=1)
    sp = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    with np.errstate(divide="ignore", invalid="ignore"):
        return (A.mean(axis=1) - B.mean(axis=1)) / sp


# ------------------------------------------------------------- window builders
def windows_fixed_width(ages: np.ndarray, midpoints: np.ndarray, half_width: float = 10.0):
    """Published scheme. Yields (mid, young_mask, old_mask) per window centre."""
    ages = np.asarray(ages, float)
    for mid in midpoints:
        y = (ages >= mid - half_width) & (ages < mid)
        o = (ages >= mid) & (ages < mid + half_width)
        yield mid, y, o


def windows_equal_n(ages: np.ndarray, midpoints: np.ndarray, k: int):
    """Power-equalised scheme: the k nearest subjects on each side of `mid`.

    Group sizes are (k, k) at every window centre by construction, so the
    sample-size component of statistical power (Artifact 3a) is held constant
    along the age axis. Windows near the edges of the age range become wide
    rather than small; the realised width is reported by `run_deswan`.
    Window centres with fewer than k subjects on a side are skipped.
    """
    ages = np.asarray(ages, float)
    for mid in midpoints:
        below = np.where(ages < mid)[0]
        above = np.where(ages >= mid)[0]
        y = np.zeros(ages.size, bool)
        o = np.zeros(ages.size, bool)
        if below.size >= k and above.size >= k:
            y[below[np.argsort(mid - ages[below])[:k]]] = True
            o[above[np.argsort(ages[above] - mid)[:k]]] = True
        yield mid, y, o


def windows_quantile(ages: np.ndarray, n_centres: int, half_width: float = 10.0):
    """DEswan package default: window centres at equally spaced age quantiles,
    fixed window *length*. Included to show that this does not fix Artifact 3a."""
    ages = np.asarray(ages, float)
    qs = np.linspace(0.1, 0.9, n_centres)
    mids = np.quantile(ages, qs)
    yield from windows_fixed_width(ages, mids, half_width)


# ------------------------------------------------------------------ main entry
def run_deswan(Y: np.ndarray, ages: np.ndarray, midpoints=None, half_width: float = 10.0,
               test: str = "wilcoxon", scheme: str = "fixed_width", k: int | None = None,
               min_n: int = 3, with_effect: bool = False, fast: bool = False):
    """Run DE-SWAN over a set of window centres.

    Parameters
    ----------
    Y : (p, n) variables x samples. May be raw per-subject data (`ages` are then
        subject ages) or a LOESS grid matrix (`ages` are then grid points).
    ages : (n,) age coordinate of every column of Y.
    midpoints : window centres. Defaults to 40..64 in 1-year steps (Shen et al.).
    half_width : half the total window width, in years (published: 10 => 20-year
        window).
    test : "wilcoxon" (Shen) or "linear" (Lehallier).
    scheme : "fixed_width", "equal_n", or "quantile".
    k : subjects per side, required when scheme == "equal_n".
    with_effect : also return per-window Cohen's d.
    fast : use `mannwhitney_p_fast` (no tie handling). Permutation loops only;
        validate with `validate_fast_mwu` on the data at hand first.

    Returns
    -------
    dict with keys
        midpoints (m,), P (p, m), Q_bh (p, m) [BH within each window centre, as
        published], n_young (m,), n_old (m,), width (m,) realised window width in
        years, and optionally D (p, m) Cohen's d.
    """
    Y = np.asarray(Y, float)
    ages = np.asarray(ages, float)
    if midpoints is None:
        midpoints = DEFAULT_MIDPOINTS
    midpoints = np.asarray(midpoints, float)

    if scheme == "fixed_width":
        wins = list(windows_fixed_width(ages, midpoints, half_width))
    elif scheme == "equal_n":
        if k is None:
            raise ValueError("scheme='equal_n' requires k")
        wins = list(windows_equal_n(ages, midpoints, k))
    elif scheme == "quantile":
        wins = list(windows_quantile(ages, len(midpoints), half_width))
        midpoints = np.array([w[0] for w in wins])
    else:
        raise ValueError(f"unknown scheme {scheme!r}")

    m = len(wins)
    P = np.full((Y.shape[0], m), np.nan)
    D = np.full((Y.shape[0], m), np.nan) if with_effect else None
    n_young = np.zeros(m, int)
    n_old = np.zeros(m, int)
    width = np.full(m, np.nan)

    if test == "wilcoxon":
        testfn = mannwhitney_p_fast if fast else mannwhitney_p
    else:
        testfn = welch_p
    for j, (mid, ymask, omask) in enumerate(wins):
        n_young[j], n_old[j] = ymask.sum(), omask.sum()
        if ymask.any() and omask.any():
            sel = ages[ymask | omask]
            width[j] = sel.max() - sel.min()
        if n_young[j] < min_n or n_old[j] < min_n:
            continue
        A, B = Y[:, ymask], Y[:, omask]
        P[:, j] = testfn(A, B)
        if with_effect:
            D[:, j] = cohens_d(A, B)

    Q = np.column_stack([bh(P[:, j]) for j in range(m)])
    out = dict(midpoints=midpoints, P=P, Q_bh=Q, n_young=n_young, n_old=n_old, width=width)
    if with_effect:
        out["D"] = D
    return out


def count_significant(Q: np.ndarray, threshold: float = 0.05) -> np.ndarray:
    """Per-window count of significant variables -- the published DE-SWAN curve."""
    return np.nansum(Q < threshold, axis=0).astype(int)


def crest_age(counts: np.ndarray, midpoints: np.ndarray) -> float:
    """Age of the maximum of a DE-SWAN curve. This is the quantity the waves
    literature interprets, and the quantity for which no sampling distribution
    has ever been stated."""
    counts = np.asarray(counts, float)
    if not np.isfinite(counts).any():
        return np.nan
    return float(np.asarray(midpoints, float)[int(np.nanargmax(counts))])


def local_maxima(counts: np.ndarray, midpoints: np.ndarray, min_prominence: float = 0.0):
    """All local maxima of a DE-SWAN curve, for cohorts where more than one
    crest is claimed (Lehallier: 34/60/78; Shen: 44/60)."""
    from scipy.signal import find_peaks
    counts = np.nan_to_num(np.asarray(counts, float))
    idx, props = find_peaks(counts, prominence=min_prominence)
    return np.asarray(midpoints, float)[idx], counts[idx], props.get("prominences", np.array([]))
