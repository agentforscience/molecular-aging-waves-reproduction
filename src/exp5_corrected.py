"""E5 -- Corrected inference on iPOP: given the pipeline is broken, does any
age-localised change survive, and where?

This is the only experiment that can REFUTE our hypothesis, and it is set up so
that it can. Five analyses, each answering a claim the literature conflates:

 5.1  Is there any age association at all?     linear F-test, BH and BY across
                                               molecules, permutation-calibrated.
                                               ("aging is real" -- weak claim)
 5.2  Is the age association nonlinear?         natural cubic spline (df=4) vs
                                               linear, nested F, same control.
                                               ("aging is nonlinear" -- middle claim)
 5.3  Is there a discrete transition, and at    segmented regression, sup-F over
      what age?                                 breakpoints 35-70, calibrated
                                               against the permutation null
                                               because the breakpoint is
                                               unidentified under H0 (Davies 1987).
                                               Breakpoint locations reported with
                                               bootstrap CIs.
                                               ("there are transition ages" and
                                                "they are 44 and 60" -- strong claims)
 5.4  Power-equalised DE-SWAN                   equal-n windows (group sizes
                                               constant by construction) on raw
                                               per-subject data, with family-wise
                                               control ACROSS the 25 overlapping
                                               window centres via a max-statistic
                                               permutation null.
 5.5  How uncertain is the published crest?     nonparametric bootstrap CI for the
                                               crest location of the published
                                               LOESS+DE-SWAN pipeline. The waves
                                               literature reports this as a point
                                               estimate with no uncertainty at all.

Outputs
-------
results/tables/e5_tests_summary.csv    per layer: counts surviving each test
results/tables/e5_breakpoints.csv      per-molecule breakpoints (significant only)
results/tables/e5_breakpoint_dist.csv  observed vs null breakpoint distribution
results/tables/e5_powerequal_deswan.csv  power-equalised curve + FWER p-values
results/tables/e5_crest_bootstrap.csv  bootstrap CI of the published crest
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from common import (complete_rows, env_report, get_logger, load_layer, save_json, save_table,
                    set_seed, zscore_rows)
from corrected import (linear_and_spline_tests, perm_pvalues,
                       permutation_null_supF, segmented_supF)
from fastdeswan import (DEFAULT_MIDPOINTS, bh, by, count_significant, crest_age,
                        run_deswan)
from fastloess import DEFAULT_GRID, LoessEngine

LOG = get_logger("exp5_corrected")

LAYERS = ["plasma_transcriptome", "plasma_proteomics", "plasma_metabolomics_metabolite",
          "plasma_lipidomics", "plasma_cytokine", "clinical_test"]
PSI_GRID = np.arange(35.0, 70.5, 1.0)
N_PERM = 200          # shared across the three tests
N_PERM_DESWAN = 500   # for the power-equalised max-statistic null
N_BOOT = 400
ALPHA = 0.05


def bootstrap_crest(Y, ages, rng, n_boot=N_BOOT):
    """Nonparametric bootstrap CI for the crest age of the published pipeline.

    Subjects are resampled with replacement and the entire pipeline (span
    selection, LOESS, DE-SWAN, BH) re-run. This is the uncertainty that the
    published point estimate of "44" omits.
    """
    crests = []
    Yz = zscore_rows(Y)
    for _ in range(n_boot):
        idx = rng.integers(0, ages.size, ages.size)
        a = ages[idx]
        if a.max() - a.min() < 40:
            continue
        try:
            eng = LoessEngine(a, grid=DEFAULT_GRID)
            M = eng.smooth(Yz[:, idx])
            res = run_deswan(M, eng.grid, DEFAULT_MIDPOINTS, 10.0, test="wilcoxon")
            cnt = count_significant(res["Q_bh"], ALPHA)
            crests.append(crest_age(cnt, DEFAULT_MIDPOINTS) if cnt.max() > 0 else np.nan)
        except Exception:
            continue
    return np.array(crests, float)


def power_equalised_deswan(Yz, ages, rng, n_perm=N_PERM_DESWAN):
    """Equal-n DE-SWAN on raw data with max-statistic FWER control.

    Group sizes are (k, k) at every window centre, so the sample-size component
    of power (Artifact 3a) cannot shape the curve. The maximum count over the 25
    overlapping window centres is calibrated against its permutation null,
    giving one FWER-controlled test of "is there ANY window with more change
    than chance allows".
    """
    n = ages.size
    k = max(10, int(np.floor(n / 6)))
    obs = run_deswan(Yz, ages, DEFAULT_MIDPOINTS, test="wilcoxon",
                     scheme="equal_n", k=k, with_effect=True)
    cnt = count_significant(obs["Q_bh"], ALPHA)

    null_max = np.zeros(n_perm, int)
    null_curves = np.zeros((n_perm, DEFAULT_MIDPOINTS.size), int)
    for b in range(n_perm):
        Yp = Yz[:, rng.permutation(n)]
        r = run_deswan(Yp, ages, DEFAULT_MIDPOINTS, test="wilcoxon", scheme="equal_n", k=k)
        c = count_significant(r["Q_bh"], ALPHA)
        null_curves[b] = c
        null_max[b] = c.max()
    p_fwer = (1 + np.sum(null_max >= cnt.max())) / (1 + n_perm)
    # Per-window p-values controlled family-wise by the max-statistic null
    p_win = np.array([(1 + np.sum(null_max >= cnt[j])) / (1 + n_perm)
                      for j in range(cnt.size)])
    return obs, cnt, p_fwer, p_win, null_curves, k


def run_layer(layer, rng):
    Y, ages, info, var_ids = load_layer(layer)
    # The least-squares battery propagates NaN, so incomplete variables are
    # excluded explicitly (they would otherwise yield NaN statistics).
    comp = complete_rows(Y)
    if not comp.all():
        LOG.info(f"[{layer}] excluding {int((~comp).sum())} variables with missing "
                 f"values from the corrected battery")
        Y, var_ids = Y[comp], var_ids[comp]
    Yz = zscore_rows(Y)
    p, n = Yz.shape
    LOG.info(f"[{layer}] {p} variables x {n} subjects -- corrected battery")

    # ---- 5.1 / 5.2 -------------------------------------------------------
    t = linear_and_spline_tests(Yz, ages)
    # ---- 5.3 -------------------------------------------------------------
    supF, psi_hat, _ = segmented_supF(Yz, ages, PSI_GRID)

    null = permutation_null_supF(Yz, ages, PSI_GRID, N_PERM, rng)
    p_lin_perm = perm_pvalues(t["F_linear"], null["F_linear"])
    p_nl_perm = perm_pvalues(t["F_nonlinear"], null["F_nonlinear"])
    p_seg_perm = perm_pvalues(supF, null["supF"])

    q_lin, q_nl, q_seg = bh(p_lin_perm), bh(p_nl_perm), bh(p_seg_perm)
    qby_lin, qby_nl, qby_seg = by(p_lin_perm), by(p_nl_perm), by(p_seg_perm)

    summary = dict(
        layer=layer, n_variables=p, n_subjects=n, n_perm=N_PERM,
        n_linear_bh=int(np.sum(q_lin < ALPHA)), n_linear_by=int(np.sum(qby_lin < ALPHA)),
        frac_linear_bh=float(np.mean(q_lin < ALPHA)),
        n_nonlinear_bh=int(np.sum(q_nl < ALPHA)), n_nonlinear_by=int(np.sum(qby_nl < ALPHA)),
        frac_nonlinear_bh=float(np.mean(q_nl < ALPHA)),
        n_segmented_bh=int(np.sum(q_seg < ALPHA)), n_segmented_by=int(np.sum(qby_seg < ALPHA)),
        frac_segmented_bh=float(np.mean(q_seg < ALPHA)),
        median_r2_linear=float(np.nanmedian(t["r2_linear"])),
        median_r2_spline=float(np.nanmedian(t["r2_spline"])),
    )

    # ---- breakpoint distribution, observed vs null -----------------------
    sig = q_seg < ALPHA
    bp_rows = []
    if sig.sum():
        for vid, ps, s, q in zip(var_ids[sig], psi_hat[sig], supF[sig], q_seg[sig]):
            bp_rows.append(dict(layer=layer, variable=str(vid), breakpoint=float(ps),
                                supF=float(s), q_value=float(q)))
    obs_psi = psi_hat[sig] if sig.sum() else np.array([])
    null_psi = null["psi"].ravel()
    dist_rows = []
    for psi in PSI_GRID:
        dist_rows.append(dict(
            layer=layer, breakpoint=float(psi),
            observed_count=int(np.sum(obs_psi == psi)),
            observed_frac=float(np.mean(obs_psi == psi)) if obs_psi.size else np.nan,
            null_frac=float(np.mean(null_psi == psi)),
        ))
    summary["n_breakpoints_significant"] = int(sig.sum())
    if sig.sum():
        summary["breakpoint_median"] = float(np.median(obs_psi))
        # Bootstrap CI for the modal/median breakpoint across molecules
        bs = [np.median(rng.choice(obs_psi, obs_psi.size, replace=True)) for _ in range(1000)]
        summary["breakpoint_median_ci_lo"] = float(np.percentile(bs, 2.5))
        summary["breakpoint_median_ci_hi"] = float(np.percentile(bs, 97.5))
        summary["frac_bp_within2_of_44"] = float(np.mean(np.abs(obs_psi - 44) <= 2))
        summary["frac_bp_within2_of_60"] = float(np.mean(np.abs(obs_psi - 60) <= 2))
        null_44 = float(np.mean(np.abs(null_psi - 44) <= 2))
        null_60 = float(np.mean(np.abs(null_psi - 60) <= 2))
        summary["null_frac_bp_within2_of_44"] = null_44
        summary["null_frac_bp_within2_of_60"] = null_60
    else:
        for kk in ("breakpoint_median", "breakpoint_median_ci_lo", "breakpoint_median_ci_hi",
                   "frac_bp_within2_of_44", "frac_bp_within2_of_60",
                   "null_frac_bp_within2_of_44", "null_frac_bp_within2_of_60"):
            summary[kk] = np.nan

    LOG.info(f"[{layer}] linear assoc: {summary['n_linear_bh']}/{p} (BH) | "
             f"nonlinear: {summary['n_nonlinear_bh']}/{p} | "
             f"segmented: {summary['n_segmented_bh']}/{p}")
    if sig.sum():
        LOG.info(f"[{layer}] breakpoints: median {summary['breakpoint_median']:.1f} "
                 f"[95% CI {summary['breakpoint_median_ci_lo']:.1f}-"
                 f"{summary['breakpoint_median_ci_hi']:.1f}], "
                 f"{100*summary['frac_bp_within2_of_44']:.0f}% within 2y of 44 "
                 f"(null {100*summary['null_frac_bp_within2_of_44']:.0f}%)")

    # ---- 5.4 power-equalised DE-SWAN --------------------------------------
    obs, cnt, p_fwer, p_win, null_curves, k = power_equalised_deswan(Yz, ages, rng)
    pe_rows = []
    for j, mid in enumerate(DEFAULT_MIDPOINTS):
        pe_rows.append(dict(layer=layer, midpoint=float(mid), k_per_side=k,
                            n_significant=int(cnt[j]),
                            n_young=int(obs["n_young"][j]), n_old=int(obs["n_old"][j]),
                            window_width_years=float(obs["width"][j]),
                            median_abs_d=float(np.nanmedian(np.abs(obs["D"][:, j]))),
                            null_mean=float(null_curves[:, j].mean()),
                            null_q975=float(np.quantile(null_curves[:, j], 0.975)),
                            p_fwer_window=float(p_win[j])))
    summary["powerequal_peak"] = int(cnt.max())
    summary["powerequal_crest"] = crest_age(cnt, DEFAULT_MIDPOINTS) if cnt.max() > 0 else np.nan
    summary["powerequal_p_fwer"] = float(p_fwer)
    summary["powerequal_k"] = k
    LOG.info(f"[{layer}] power-equalised DE-SWAN (k={k}/side): peak {cnt.max()}/{p} "
             f"at {summary['powerequal_crest']}, FWER p = {p_fwer:.3f}")

    # ---- 5.5 bootstrap CI of the published crest --------------------------
    bc = bootstrap_crest(Y, ages, rng)
    boot_row = dict(layer=layer, n_boot=int(np.isfinite(bc).sum()),
                    crest_point_estimate=float(crest_age(
                        count_significant(run_deswan(
                            LoessEngine(ages, grid=DEFAULT_GRID).smooth(Yz),
                            DEFAULT_GRID, DEFAULT_MIDPOINTS, 10.0)["Q_bh"], ALPHA),
                        DEFAULT_MIDPOINTS)),
                    boot_median=float(np.nanmedian(bc)),
                    ci_lo=float(np.nanpercentile(bc, 2.5)),
                    ci_hi=float(np.nanpercentile(bc, 97.5)),
                    ci_width_years=float(np.nanpercentile(bc, 97.5) - np.nanpercentile(bc, 2.5)),
                    frac_in_44pm2=float(np.nanmean(np.abs(bc - 44) <= 2)),
                    frac_in_60pm2=float(np.nanmean(np.abs(bc - 60) <= 2)))
    LOG.info(f"[{layer}] published crest {boot_row['crest_point_estimate']:.0f}, "
             f"bootstrap 95% CI [{boot_row['ci_lo']:.0f}, {boot_row['ci_hi']:.0f}] "
             f"({boot_row['ci_width_years']:.0f} years wide)")

    return summary, bp_rows, dist_rows, pe_rows, boot_row


def main() -> None:
    rng = set_seed()
    summaries, bps, dists, pes, boots = [], [], [], [], []
    for layer in LAYERS:
        try:
            s, b, d, pe, bt = run_layer(layer, rng)
            summaries.append(s)
            bps += b
            dists += d
            pes += pe
            boots.append(bt)
        except Exception as exc:
            LOG.error(f"[{layer}] FAILED: {exc!r}")
            import traceback
            LOG.error(traceback.format_exc())

    save_table(pd.DataFrame(summaries), "e5_tests_summary.csv", LOG)
    save_table(pd.DataFrame(bps), "e5_breakpoints.csv", LOG)
    save_table(pd.DataFrame(dists), "e5_breakpoint_dist.csv", LOG)
    save_table(pd.DataFrame(pes), "e5_powerequal_deswan.csv", LOG)
    save_table(pd.DataFrame(boots), "e5_crest_bootstrap.csv", LOG)
    save_json({"psi_grid": PSI_GRID.tolist(), "n_perm": N_PERM,
               "n_perm_deswan": N_PERM_DESWAN, "n_boot": N_BOOT, "alpha": ALPHA,
               "env": env_report()}, "e5_meta.json", LOG)


if __name__ == "__main__":
    main()
