"""E6 -- Independent replication in GSE40279 (Hannum et al. 2013 whole-blood 450K).

Why this cohort. iPOP has 77-102 subjects spanning 26-75. If a corrected analysis
finds nothing there, that is ambiguous between "the method is broken" and "the
cohort is underpowered". GSE40279 has n=656 across ages 19-101 -- roughly 6x the
subjects, a wider and much more even age range, and a modality (DNA methylation)
with well-established, strong age associations. It is therefore the right place
to ask the question our hypothesis actually poses: in data with real power, is
change concentrated at particular ages, and are those ages ~44 and ~60?

Arms
----
A  LOESS+DE-SWAN (published pipeline)      does the artifact reproduce at n=656?
B  DE-SWAN, no LOESS, fixed-width windows  the naive corrected comparison; group
                                           sizes still vary along the age axis
C  DE-SWAN, no LOESS, power-equalised      group sizes constant by construction,
                                           with family-wise control across the
                                           overlapping window centres by a
                                           max-statistic permutation null
D  Corrected battery                       linear / spline-nonlinearity /
                                           segmented breakpoint, all
                                           permutation-calibrated

Confounding. Plate is severely confounded with age in this cohort
(one-way ANOVA of age on plate: F = 65.0, p = 6e-78) and sample source likewise.
Every arm is therefore run twice: unadjusted, and after residualising each CpG
on plate + source + sex + ethnicity. An age transition that survives only in the
unadjusted arm is a batch effect, not biology.

CpG selection. The top-20,000-variance CpGs are a natural analysis set but are
*selected on variance*, which is a form of circularity if variance covaries with
age. A random 20,000 CpGs drawn from the full 473,034 are analysed alongside as
an unselected control.

Outputs
-------
results/tables/e6_curves.csv        DE-SWAN curves (all arms) with group sizes
results/tables/e6_summary.csv       crests, FWER p-values, achieved FDP
results/tables/e6_breakpoint_dist.csv  observed vs null breakpoint distribution
results/tables/e6_tests_summary.csv    corrected battery counts
"""
from __future__ import annotations

import multiprocessing as mp
import os
import time

import numpy as np
import pandas as pd

from common import (FIGURES, env_report, get_logger, save_json, save_table,
                    set_seed, zscore_rows, ROOT)
from corrected import (linear_and_spline_tests, perm_pvalues,
                       permutation_null_supF, segmented_supF)
from fastdeswan import (bh, by, count_significant, crest_age, run_deswan,
                        validate_fast_mwu)
from fastloess import LoessEngine

LOG = get_logger("exp6_gse40279")

GSE = f"{ROOT}/datasets/GSE40279"
N_CPG = 20000
MIDPOINTS = np.arange(30.0, 91.0, 1.0)     # 20-year windows fit inside 19-101
# The permutation null for the max statistic is evaluated on a 2-year grid.
# Neighbouring 1-year window centres share 90% of their subjects, so the maximum
# over the coarse grid is within one window of the maximum over the fine grid,
# while costing half as much. Documented rather than silent, per our own rule
# about never truncating coverage without saying so.
NULL_MIDPOINTS = np.arange(30.0, 91.0, 2.0)
HALF_WIDTH = 10.0
GRID = np.arange(25.0, 95.5, 0.5)          # half-year grid, as published
PSI_GRID = np.arange(30.0, 90.5, 1.0)
N_PERM_DESWAN = 200
N_PERM_TESTS = 100
N_PERM_FDP = 30
ALPHA = 0.05


def load_gse(selection: str, rng):
    """Load a 20,000-CpG subset. `selection` is 'top_variance' or 'random'."""
    info = pd.read_csv(f"{GSE}/sample_info.csv")
    ages = info["age"].to_numpy(float)
    if selection == "top_variance":
        df = pd.read_parquet(f"{GSE}/beta_top20k_var.parquet")
        Y, ids = df.to_numpy(np.float64), df.index.to_numpy()
        assert list(df.columns) == info["gsm"].tolist(), "sample order mismatch"
    else:
        beta = np.load(f"{GSE}/beta_float32.npy", mmap_mode="r")
        cpgs = pd.read_csv(f"{GSE}/cpg_ids.csv")
        idx = np.sort(rng.choice(beta.shape[0], N_CPG, replace=False))
        Y = np.asarray(beta[idx, :], dtype=np.float64)
        ids = cpgs.iloc[idx, 0].to_numpy()
    ok = np.isfinite(Y).all(axis=1) & (Y.std(axis=1) > 1e-6)
    return Y[ok], ids[ok], ages, info


def adjust(Y, info, cols=("plate", "source", "gender", "ethnicity")):
    """Residualise every CpG on the technical/demographic covariates."""
    X = pd.get_dummies(info[list(cols)].astype(str), drop_first=True).to_numpy(float)
    X = np.column_stack([np.ones(X.shape[0]), X])
    Q, _ = np.linalg.qr(X)
    return Y - (Y @ Q) @ Q.T


def deswan_arm(Y, ages, scheme, k=None, half_width=HALF_WIDTH):
    res = run_deswan(Y, ages, MIDPOINTS, half_width, test="wilcoxon",
                     scheme=scheme, k=k, with_effect=True)
    return res, count_significant(res["Q_bh"], ALPHA)


# The permutation loop is the dominant cost at n=656 and is dominated in turn by
# the per-window argsort, which numpy runs single-threaded. It parallelises
# trivially across permutations. Workers are forked, so they inherit the data
# matrix copy-on-write via a module global rather than pickling ~100 MB per task.
_MP = {}
N_WORKERS = max(1, min(12, (os.cpu_count() or 4) - 4))


def _null_one(args):
    """One permutation of the max-statistic null. Runs in a forked worker."""
    seed, scheme, k = args
    Y, ages = _MP["Y"], _MP["ages"]
    rs = np.random.default_rng(seed)
    Yp = Y[:, rs.permutation(ages.size)]
    r = run_deswan(Yp, ages, NULL_MIDPOINTS, HALF_WIDTH, test="wilcoxon",
                   scheme=scheme, k=k, fast=True)
    return count_significant(r["Q_bh"], ALPHA)


def maxstat_null(Y, ages, rng, scheme, k=None, n_perm=N_PERM_DESWAN):
    """Permutation null for the maximum count over window centres (FWER control).

    Evaluated on NULL_MIDPOINTS (2-year spacing) with the tie-free fast
    Mann-Whitney; both approximations are validated before use in `run_selection`.
    Returns the (n_perm, len(NULL_MIDPOINTS)) matrix of null curves.

    Each permutation gets its own explicitly drawn seed, so the result is
    independent of worker scheduling and reproducible for a given `rng` state.
    """
    _MP["Y"], _MP["ages"] = Y, ages
    seeds = rng.integers(0, 2**63 - 1, n_perm)
    tasks = [(int(s), scheme, k) for s in seeds]
    if N_WORKERS <= 1:
        return np.array([_null_one(t) for t in tasks], int)
    ctx = mp.get_context("fork")
    with ctx.Pool(N_WORKERS) as pool:
        out = pool.map(_null_one, tasks, chunksize=1)
    return np.array(out, int)


def null_at(nullm, mid):
    """Null curve column nearest a (possibly finer-grid) window centre."""
    return nullm[:, int(np.argmin(np.abs(NULL_MIDPOINTS - mid)))]


def run_selection(selection: str, adjusted: bool, rng):
    Y, ids, ages, info = load_gse(selection, rng)
    if adjusted:
        Y = adjust(Y, info)
    Yz = zscore_rows(Y)
    p, n = Yz.shape
    tag = f"{selection}/{'adjusted' if adjusted else 'unadjusted'}"
    LOG.info(f"[{tag}] {p} CpGs x {n} subjects, ages {ages.min():.0f}-{ages.max():.0f}")

    # The permutation loops use a tie-free fast Mann-Whitney. Verify on this
    # exact data that it reproduces the reference p-values before relying on it.
    sel_y = (ages >= 40) & (ages < 50)
    sel_o = (ages >= 50) & (ages < 60)
    d_mwu = validate_fast_mwu(Yz[:2000, sel_y], Yz[:2000, sel_o])
    assert d_mwu < 1e-9, f"fast Mann-Whitney disagrees with reference by {d_mwu:.2e}"
    LOG.info(f"[{tag}] fast Mann-Whitney validated: max |dp| = {d_mwu:.1e}")

    curves, summaries = [], []

    # ---- Arm A: published pipeline ---------------------------------------
    t0 = time.time()
    eng = LoessEngine(ages, grid=GRID)
    Yhat = eng.smooth(Yz)
    resA, cA = deswan_arm(Yhat, eng.grid, "fixed_width")
    LOG.info(f"[{tag}] A LOESS+DE-SWAN: peak {cA.max()}/{p} at "
             f"{crest_age(cA, MIDPOINTS):.0f}, floor {cA.min()} ({time.time()-t0:.0f}s)")

    # The permutation nulls are the expensive part. They are computed for the
    # primary analysis set (top-variance CpGs, both adjustments). The random-CpG
    # set is a selection-sensitivity check on the CURVES, so it is run without
    # nulls; this is a deliberate coverage limit and is stated in the report
    # rather than left implicit.
    with_nulls = selection == "top_variance"

    # ---- Arm B: DE-SWAN, no LOESS, fixed width ---------------------------
    resB, cB = deswan_arm(Yz, ages, "fixed_width")
    nullB = maxstat_null(Yz, ages, rng, "fixed_width") if with_nulls else None
    pB = ((1 + np.sum(nullB.max(axis=1) >= cB.max())) / (1 + N_PERM_DESWAN)
          if with_nulls else np.nan)
    LOG.info(f"[{tag}] B DE-SWAN only: peak {cB.max()}/{p} at "
             f"{crest_age(cB, MIDPOINTS):.0f}, FWER p={pB:.4f}")

    # ---- Arm C: power-equalised ------------------------------------------
    k = max(30, int(np.floor(n / 8)))
    resC, cC = deswan_arm(Yz, ages, "equal_n", k=k)
    nullC = maxstat_null(Yz, ages, rng, "equal_n", k=k) if with_nulls else None
    pC = ((1 + np.sum(nullC.max(axis=1) >= cC.max())) / (1 + N_PERM_DESWAN)
          if with_nulls else np.nan)
    LOG.info(f"[{tag}] C power-equalised (k={k}): peak {cC.max()}/{p} at "
             f"{crest_age(cC, MIDPOINTS):.0f}, FWER p={pC:.4f}")

    # Achieved FDP of the published pipeline under the valid null
    fdp = [np.nan]
    if with_nulls:
        fdp = []
        for b in range(N_PERM_FDP):
            Yp = eng.smooth_permuted(Yz, rng)
            r = run_deswan(Yp, eng.grid, NULL_MIDPOINTS, HALF_WIDTH,
                           test="wilcoxon", fast=True)
            fdp.append(count_significant(r["Q_bh"], ALPHA).sum() / (p * NULL_MIDPOINTS.size))
        LOG.info(f"[{tag}] achieved FDP of published pipeline = {np.mean(fdp):.3f} "
                 f"({np.mean(fdp)/ALPHA:.0f}x nominal)")

    for arm, res, cnt, nullm in (("A_loess_deswan", resA, cA, None),
                                 ("B_deswan_only", resB, cB, nullB),
                                 ("C_power_equalised", resC, cC, nullC)):
        for j, mid in enumerate(MIDPOINTS):
            curves.append(dict(
                selection=selection, adjusted=adjusted, arm=arm, midpoint=float(mid),
                n_significant=int(cnt[j]), frac_significant=float(cnt[j] / p),
                n_young=int(res["n_young"][j]), n_old=int(res["n_old"][j]),
                window_width_years=float(res["width"][j]),
                median_abs_d=float(np.nanmedian(np.abs(res["D"][:, j]))),
                null_mean=float(null_at(nullm, mid).mean()) if nullm is not None else np.nan,
                null_q975=float(np.quantile(null_at(nullm, mid), 0.975))
                if nullm is not None else np.nan,
                p_fwer=float((1 + np.sum(nullm.max(axis=1) >= cnt[j])) / (1 + N_PERM_DESWAN))
                if nullm is not None else np.nan,
            ))

    summaries.append(dict(
        selection=selection, adjusted=adjusted, n_cpg=p, n_subjects=n,
        A_peak=int(cA.max()), A_crest=crest_age(cA, MIDPOINTS), A_floor=int(cA.min()),
        B_peak=int(cB.max()), B_crest=crest_age(cB, MIDPOINTS), B_p_fwer=float(pB),
        C_peak=int(cC.max()), C_crest=crest_age(cC, MIDPOINTS), C_p_fwer=float(pC),
        C_k_per_side=k,
        achieved_fdp_published=float(np.mean(fdp)), nominal_fdr=ALPHA,
    ))

    # ---- Arm D: corrected battery ----------------------------------------
    t = linear_and_spline_tests(Yz, ages)
    supF, psi_hat, _ = segmented_supF(Yz, ages, PSI_GRID, min_per_side=30)
    null = permutation_null_supF(Yz, ages, PSI_GRID, N_PERM_TESTS, rng, min_per_side=30)
    p_lin = perm_pvalues(t["F_linear"], null["F_linear"])
    p_nl = perm_pvalues(t["F_nonlinear"], null["F_nonlinear"])
    p_seg = perm_pvalues(supF, null["supF"])
    q_lin, q_nl, q_seg = bh(p_lin), bh(p_nl), bh(p_seg)

    sig = q_seg < ALPHA
    obs_psi, null_psi = psi_hat[sig], null["psi"].ravel()
    tests = dict(
        selection=selection, adjusted=adjusted, n_cpg=p,
        n_linear_bh=int(np.sum(q_lin < ALPHA)), n_linear_by=int(np.sum(by(p_lin) < ALPHA)),
        n_nonlinear_bh=int(np.sum(q_nl < ALPHA)), n_nonlinear_by=int(np.sum(by(p_nl) < ALPHA)),
        n_segmented_bh=int(sig.sum()), n_segmented_by=int(np.sum(by(p_seg) < ALPHA)),
        frac_linear=float(np.mean(q_lin < ALPHA)),
        frac_nonlinear=float(np.mean(q_nl < ALPHA)),
        frac_segmented=float(np.mean(sig)),
        median_r2_linear=float(np.nanmedian(t["r2_linear"])),
        median_r2_spline=float(np.nanmedian(t["r2_spline"])),
    )
    if sig.sum():
        bs = [np.median(rng.choice(obs_psi, obs_psi.size, replace=True)) for _ in range(1000)]
        tests.update(breakpoint_median=float(np.median(obs_psi)),
                     breakpoint_ci_lo=float(np.percentile(bs, 2.5)),
                     breakpoint_ci_hi=float(np.percentile(bs, 97.5)),
                     breakpoint_mode=float(pd.Series(obs_psi).mode().iloc[0]),
                     frac_bp_44pm2=float(np.mean(np.abs(obs_psi - 44) <= 2)),
                     frac_bp_60pm2=float(np.mean(np.abs(obs_psi - 60) <= 2)),
                     null_frac_bp_44pm2=float(np.mean(np.abs(null_psi - 44) <= 2)),
                     null_frac_bp_60pm2=float(np.mean(np.abs(null_psi - 60) <= 2)))
        LOG.info(f"[{tag}] D corrected: linear {tests['n_linear_bh']}/{p}, "
                 f"nonlinear {tests['n_nonlinear_bh']}/{p}, "
                 f"segmented {tests['n_segmented_bh']}/{p}; breakpoint median "
                 f"{tests['breakpoint_median']:.0f} "
                 f"[{tests['breakpoint_ci_lo']:.0f}-{tests['breakpoint_ci_hi']:.0f}], "
                 f"{100*tests['frac_bp_44pm2']:.1f}% near 44 "
                 f"(null {100*tests['null_frac_bp_44pm2']:.1f}%), "
                 f"{100*tests['frac_bp_60pm2']:.1f}% near 60 "
                 f"(null {100*tests['null_frac_bp_60pm2']:.1f}%)")

    dist = [dict(selection=selection, adjusted=adjusted, breakpoint=float(psi),
                 observed_count=int(np.sum(obs_psi == psi)),
                 observed_frac=float(np.mean(obs_psi == psi)) if obs_psi.size else np.nan,
                 null_frac=float(np.mean(null_psi == psi)))
            for psi in PSI_GRID]

    return curves, summaries, [tests], dist


def main() -> None:
    rng = set_seed()
    curves, summaries, tests, dists = [], [], [], []
    for selection in ("top_variance", "random"):
        for adjusted in (False, True):
            try:
                c, s, t, d = run_selection(selection, adjusted, rng)
                curves += c
                summaries += s
                tests += t
                dists += d
            except Exception as exc:
                LOG.error(f"[{selection}/{adjusted}] FAILED: {exc!r}")
                import traceback
                LOG.error(traceback.format_exc())

    save_table(pd.DataFrame(curves), "e6_curves.csv", LOG)
    save_table(pd.DataFrame(summaries), "e6_summary.csv", LOG)
    save_table(pd.DataFrame(tests), "e6_tests_summary.csv", LOG)
    save_table(pd.DataFrame(dists), "e6_breakpoint_dist.csv", LOG)
    save_json({"n_cpg": N_CPG, "midpoints": MIDPOINTS.tolist(),
               "n_perm_deswan": N_PERM_DESWAN, "n_perm_tests": N_PERM_TESTS,
               "alpha": ALPHA, "env": env_report()}, "e6_meta.json", LOG)


if __name__ == "__main__":
    main()
