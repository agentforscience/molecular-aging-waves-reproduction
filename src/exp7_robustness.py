"""E7 -- Robustness of the one positive finding, and outlier diagnostics.

E5 returned a single positive result: in the iPOP plasma cytokine panel, 46/66
analytes showed spline nonlinearity and 35/66 a significant breakpoint, all
placed at age 41. Before that is reported as a surviving transition it has to be
interrogated, because it has the signature of the very failure mode this project
is about. Manual inspection of the strongest hit (IL23, sup-F = 48.5) showed the
entire effect coming from ONE subject, aged 25.9, at z = +9.4; 83% of cytokine
analytes have some |z| > 4. Raw cytokine concentrations are strongly
right-skewed and were not log-transformed in our reconstruction.

This script therefore:

  7.1  Quantifies outlier burden per layer (max |z|, fraction of variables with
       an extreme point, excess kurtosis).
  7.2  Re-runs the corrected battery under a rank-based inverse-normal transform
       (RINT), which is monotone -- so it preserves any genuine monotone or
       ordered relationship with age -- while bounding the influence of any one
       observation. If the cytokine transition is real, it survives; if it is one
       subject, it does not.
  7.3  Checks whether estimated breakpoints pile up at the EDGE of the admissible
       breakpoint grid, which is the fingerprint of an outlier near the age range
       boundary rather than of a transition in the interior.
  7.4  Leave-one-subject-out influence: refit and record how often the finding
       disappears when a single subject is dropped.

Outputs
-------
results/tables/e7_outlier_diagnostics.csv
results/tables/e7_rint_battery.csv
results/tables/e7_influence.csv
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from common import (complete_rows, env_report, get_logger, load_layer, save_json,
                    save_table,
                    set_seed, zscore_rows)
from corrected import (linear_and_spline_tests, perm_pvalues,
                       permutation_null_supF, segmented_supF)
from fastdeswan import bh

LOG = get_logger("exp7_robustness")

LAYERS = ["plasma_transcriptome", "plasma_proteomics", "plasma_metabolomics_metabolite",
          "plasma_lipidomics", "plasma_cytokine", "clinical_test"]
PSI_GRID = np.arange(35.0, 70.5, 1.0)
N_PERM = 200
ALPHA = 0.05


def rint(Y: np.ndarray) -> np.ndarray:
    """Rank-based inverse-normal transform, applied per variable.

    y -> Phi^-1((rank - 3/8) / (n + 1/4))  (Blom scores)

    Monotone, so any genuinely ordered relationship with age is preserved, but
    the influence of an extreme observation is bounded by the normal quantile of
    the top rank rather than by its raw magnitude. This is the standard
    robustness check for a finding that might be one data point.
    """
    n = Y.shape[1]
    r = stats.rankdata(Y, axis=1, method="average")
    return stats.norm.ppf((r - 0.375) / (n + 0.25))


def diagnostics(Y: np.ndarray, layer: str) -> dict:
    Yz = zscore_rows(Y)
    mx = np.abs(Yz).max(axis=1)
    kurt = stats.kurtosis(Yz, axis=1, fisher=True, bias=False)
    return dict(layer=layer, n_variables=Yz.shape[0], n_subjects=Yz.shape[1],
                median_max_abs_z=float(np.median(mx)),
                max_abs_z=float(mx.max()),
                frac_vars_with_z_gt_4=float(np.mean(mx > 4)),
                frac_vars_with_z_gt_6=float(np.mean(mx > 6)),
                median_excess_kurtosis=float(np.nanmedian(kurt)),
                frac_vars_kurtosis_gt_5=float(np.nanmean(kurt > 5)))


def battery(Y, ages, rng, label, layer):
    """Run the E5 corrected battery on an arbitrary transform of the data."""
    t = linear_and_spline_tests(Y, ages)
    supF, psi_hat, _ = segmented_supF(Y, ages, PSI_GRID)
    null = permutation_null_supF(Y, ages, PSI_GRID, N_PERM, rng)
    q_lin = bh(perm_pvalues(t["F_linear"], null["F_linear"]))
    q_nl = bh(perm_pvalues(t["F_nonlinear"], null["F_nonlinear"]))
    q_seg = bh(perm_pvalues(supF, null["supF"]))
    sig = q_seg < ALPHA
    row = dict(layer=layer, transform=label, n_variables=Y.shape[0],
               n_linear=int((q_lin < ALPHA).sum()),
               n_nonlinear=int((q_nl < ALPHA).sum()),
               n_segmented=int(sig.sum()))
    if sig.sum():
        bp = psi_hat[sig]
        row.update(breakpoint_median=float(np.median(bp)),
                   breakpoint_mode=float(pd.Series(bp).mode().iloc[0]),
                   # A pile-up at the first or last admissible breakpoint is the
                   # fingerprint of an edge outlier, not an interior transition.
                   frac_at_grid_edge=float(np.mean((bp <= PSI_GRID[1]) |
                                                   (bp >= PSI_GRID[-2]))),
                   frac_within2_of_44=float(np.mean(np.abs(bp - 44) <= 2)),
                   frac_within2_of_60=float(np.mean(np.abs(bp - 60) <= 2)))
    else:
        row.update(breakpoint_median=np.nan, breakpoint_mode=np.nan,
                   frac_at_grid_edge=np.nan, frac_within2_of_44=np.nan,
                   frac_within2_of_60=np.nan)
    return row, sig, psi_hat


def influence(Y, ages, sig_idx, rng, layer, n_drop=None):
    """Leave-one-subject-out: how many findings survive dropping one subject?

    Reports, for each subject, the number of previously-significant variables
    that remain significant when that subject is removed. A finding that
    collapses when one particular subject is dropped is that subject.
    """
    if sig_idx.sum() == 0:
        return []
    Ysig = Y[sig_idx]
    n = ages.size
    base_supF, _, _ = segmented_supF(Ysig, ages, PSI_GRID)
    rows = []
    for i in range(n if n_drop is None else min(n, n_drop)):
        keep = np.ones(n, bool)
        keep[i] = False
        s, _, _ = segmented_supF(Ysig[:, keep], ages[keep], PSI_GRID)
        rows.append(dict(layer=layer, dropped_subject_index=int(i),
                         dropped_subject_age=float(ages[i]),
                         median_supF_full=float(np.median(base_supF)),
                         median_supF_dropped=float(np.median(s)),
                         supF_retention=float(np.median(s) / max(np.median(base_supF), 1e-9))))
    return rows


def main() -> None:
    rng = set_seed()
    diag, bat, infl = [], [], []

    for layer in LAYERS:
        Y, ages, info, _ = load_layer(layer)
        comp = complete_rows(Y)          # see common.complete_rows: NaN propagates
        if not comp.all():
            LOG.info(f"[{layer}] excluding {int((~comp).sum())} incomplete variables")
            Y = Y[comp]
        diag.append(diagnostics(Y, layer))
        LOG.info(f"[{layer}] max|z| {diag[-1]['max_abs_z']:.1f}, "
                 f"{100*diag[-1]['frac_vars_with_z_gt_4']:.0f}% of variables have a "
                 f"point beyond 4 SD, median excess kurtosis "
                 f"{diag[-1]['median_excess_kurtosis']:.1f}")

        Yz = zscore_rows(Y)
        r_raw, sig_raw, _ = battery(Yz, ages, rng, "zscore (as in E5)", layer)
        r_rint, sig_rint, psi_rint = battery(rint(Y), ages, rng, "rank-inverse-normal", layer)
        bat += [r_raw, r_rint]
        LOG.info(f"[{layer}] segmented transitions: {r_raw['n_segmented']} (z-score) -> "
                 f"{r_rint['n_segmented']} (RINT); nonlinear "
                 f"{r_raw['n_nonlinear']} -> {r_rint['n_nonlinear']}")
        if r_raw["n_segmented"]:
            LOG.info(f"[{layer}]   z-score breakpoints: mode "
                     f"{r_raw['breakpoint_mode']:.0f}, "
                     f"{100*r_raw['frac_at_grid_edge']:.0f}% at the grid edge")

        # Influence analysis only where something was found, and only on the
        # layers small enough for a full leave-one-out sweep.
        if sig_raw.sum() and Yz.shape[0] <= 2000:
            infl += influence(Yz, ages, sig_raw, rng, layer)

    save_table(pd.DataFrame(diag), "e7_outlier_diagnostics.csv", LOG)
    save_table(pd.DataFrame(bat), "e7_rint_battery.csv", LOG)
    if infl:
        df = pd.DataFrame(infl)
        save_table(df, "e7_influence.csv", LOG)
        for layer, s in df.groupby("layer"):
            w = s.loc[s.supF_retention.idxmin()]
            LOG.info(f"[{layer}] most influential subject: age {w.dropped_subject_age:.0f}; "
                     f"dropping it alone leaves {100*w.supF_retention:.0f}% of the "
                     f"median sup-F statistic")
    save_json({"psi_grid": PSI_GRID.tolist(), "n_perm": N_PERM, "alpha": ALPHA,
               "env": env_report()}, "e7_meta.json", LOG)


if __name__ == "__main__":
    main()
