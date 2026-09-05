"""E2 -- Pipeline ablation: which step manufactures the crest?

"The pipeline is broken" is not actionable until the responsible component is
named. Two independent factors are conflated in LOESS+DE-SWAN:

  (i)  SMOOTHING     -- each fitted value is a weighted average of overlapping
       neighbourhoods, so fitted values are strongly positively dependent and
       nearly monotone within a window. A rank test on them separates almost
       perfectly.
  (ii) DENSIFICATION -- the fit is evaluated on a 99-point half-year grid rather
       than at the 77-102 observed ages, and the grid is *uniform* in age while
       the subjects are not.

These are separable, and nobody has separated them. We therefore run:

  raw              raw per-subject values                     (neither)
  interp_only      linear interpolation of raw values onto     (densification only)
                   the grid -- no local averaging
  loess_at_ages    LOESS evaluated at the observed ages        (smoothing only)
  loess_grid       LOESS evaluated on the grid  = PUBLISHED    (both)

plus one-factor-at-a-time ablations of the remaining published choices:
test statistic (Wilcoxon vs Lehallier's linear model), window scheme
(fixed width vs quantile centres vs our power-equalised equal-n), bucket width,
q-threshold, LOESS span, and confounder adjustment.

Outputs
-------
results/tables/e2_components.csv   the 4-way smoothing/densification decomposition
results/tables/e2_onefactor.csv    one-factor-at-a-time ablations
results/tables/e2_width_q_grid.csv the published robustness grid (width x q)
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from common import (FAST_LAYERS, env_report, get_logger, load_layer, regress_out,
                    save_json, save_table, set_seed, zscore_rows)
from fastdeswan import DEFAULT_MIDPOINTS, count_significant, crest_age, run_deswan
from fastloess import DEFAULT_GRID, LoessEngine

LOG = get_logger("exp2_ablation")
LAYERS = ["plasma_transcriptome"] + FAST_LAYERS


def _summarise(res, cnt, p, **tags) -> dict:
    """Common summary row for one DE-SWAN curve."""
    mids = res["midpoints"]
    return dict(
        **tags, n_variables=p,
        crest_age=crest_age(cnt, mids) if cnt.max() > 0 else np.nan,
        peak_n=int(cnt.max()), peak_frac=float(cnt.max() / p),
        min_n=int(cnt.min()), min_frac=float(cnt.min() / p),
        total_significant=int(cnt.sum()),
        n_windows_with_any=int((cnt > 0).sum()),
        # A curve that is uniformly near-saturated carries no localisation
        # information at all; this ratio distinguishes "a peak" from "a plateau".
        peak_to_floor=float(cnt.max() / cnt.min()) if cnt.min() > 0 else np.inf,
        n_young_min=int(res["n_young"].min()), n_young_max=int(res["n_young"].max()),
        n_old_min=int(res["n_old"].min()), n_old_max=int(res["n_old"].max()),
        median_width=float(np.nanmedian(res["width"])),
    )


def build_representations(Yz, ages, eng):
    """The four smoothing x densification cells."""
    grid = eng.grid
    # Linear interpolation onto the grid: densification with no local averaging.
    order = np.argsort(ages)
    xs = ages[order]
    interp = np.vstack([np.interp(grid, xs, row[order]) for row in Yz])
    # LOESS evaluated at the observed ages: smoothing with no densification.
    eng_at_ages = LoessEngine(ages, grid=np.sort(ages))
    spans, _ = eng_at_ages.select_spans(Yz)
    loess_ages = eng_at_ages.smooth(Yz, spans=spans)
    return {
        "raw":           (Yz, ages, "neither"),
        "interp_only":   (interp, grid, "densification only"),
        "loess_at_ages": (loess_ages, np.sort(ages), "smoothing only"),
        "loess_grid":    (eng.smooth(Yz), grid, "both (PUBLISHED)"),
    }


def main() -> None:
    set_seed()
    comp_rows, one_rows, grid_rows = [], [], []

    for layer in LAYERS:
        Y, ages, info, _ = load_layer(layer)
        Yz = zscore_rows(Y)
        p = Yz.shape[0]
        eng = LoessEngine(ages, grid=DEFAULT_GRID)
        LOG.info(f"[{layer}] {p} variables x {Yz.shape[1]} subjects")

        # -------- A. smoothing x densification decomposition ------------------
        reps = build_representations(Yz, ages, eng)
        for name, (M, coord, label) in reps.items():
            res = run_deswan(M, coord, DEFAULT_MIDPOINTS, 10.0, test="wilcoxon")
            cnt = count_significant(res["Q_bh"], 0.05)
            comp_rows.append(_summarise(res, cnt, p, layer=layer,
                                        representation=name, component=label))
            LOG.info(f"  [{layer}] {name:14s} ({label:20s}) peak {cnt.max():6d}/{p} "
                     f"floor {cnt.min():6d} crest "
                     f"{comp_rows[-1]['crest_age']}")

        Yhat = reps["loess_grid"][0]

        # -------- B. one-factor-at-a-time ablations ---------------------------
        # B1 test statistic (on the published representation)
        for test in ("wilcoxon", "linear"):
            res = run_deswan(Yhat, eng.grid, DEFAULT_MIDPOINTS, 10.0, test=test)
            cnt = count_significant(res["Q_bh"], 0.05)
            one_rows.append(_summarise(res, cnt, p, layer=layer, factor="test",
                                       level=test, arm="loess_grid"))
            # and on raw data, to show the choice of test is not what matters
            res = run_deswan(Yz, ages, DEFAULT_MIDPOINTS, 10.0, test=test)
            cnt = count_significant(res["Q_bh"], 0.05)
            one_rows.append(_summarise(res, cnt, p, layer=layer, factor="test",
                                       level=test, arm="raw"))

        # B2 window scheme
        k = max(8, int(np.floor(np.median([
            ((ages >= m - 10) & (ages < m)).sum() for m in DEFAULT_MIDPOINTS]))))
        for scheme, kw in (("fixed_width", {}), ("quantile", {}), ("equal_n", dict(k=k))):
            for arm, M, coord in (("loess_grid", Yhat, eng.grid), ("raw", Yz, ages)):
                kk = kw.copy()
                if scheme == "equal_n" and arm == "loess_grid":
                    kk["k"] = 20      # grid points per side under the published width
                res = run_deswan(M, coord, DEFAULT_MIDPOINTS, 10.0,
                                 test="wilcoxon", scheme=scheme, **kk)
                cnt = count_significant(res["Q_bh"], 0.05)
                one_rows.append(_summarise(res, cnt, p, layer=layer, factor="scheme",
                                           level=scheme, arm=arm))

        # B3 LOESS span held fixed instead of chosen by CV
        for span in (0.3, 0.4, 0.5, 0.6):
            M = eng.smooth(Yz, spans=np.full(p, span))
            res = run_deswan(M, eng.grid, DEFAULT_MIDPOINTS, 10.0, test="wilcoxon")
            cnt = count_significant(res["Q_bh"], 0.05)
            one_rows.append(_summarise(res, cnt, p, layer=layer, factor="span",
                                       level=str(span), arm="loess_grid"))

        # B4 confounder adjustment (sex / ethnicity / IRIS), flagged in STATE.md
        # as an untested assumption of our reconstruction
        Yr, used = regress_out(Y, info)
        if used:
            Yrz = zscore_rows(Yr)
            for arm, M, coord in (("loess_grid", eng.smooth(Yrz), eng.grid),
                                  ("raw", Yrz, ages)):
                res = run_deswan(M, coord, DEFAULT_MIDPOINTS, 10.0, test="wilcoxon")
                cnt = count_significant(res["Q_bh"], 0.05)
                one_rows.append(_summarise(res, cnt, p, layer=layer,
                                           factor="confounder_adjusted",
                                           level="+".join(used), arm=arm))

        # -------- C. the published robustness grid: bucket width x q ---------
        for bucket in (15, 20, 25, 30):
            for arm, M, coord in (("loess_grid", Yhat, eng.grid), ("raw", Yz, ages)):
                res = run_deswan(M, coord, DEFAULT_MIDPOINTS, bucket / 2.0, test="wilcoxon")
                for q_thr in (1e-4, 1e-3, 1e-2, 0.05):
                    cnt = count_significant(res["Q_bh"], q_thr)
                    grid_rows.append(_summarise(res, cnt, p, layer=layer, arm=arm,
                                                bucket=bucket, q_threshold=q_thr))

    save_table(pd.DataFrame(comp_rows), "e2_components.csv", LOG)
    save_table(pd.DataFrame(one_rows), "e2_onefactor.csv", LOG)
    save_table(pd.DataFrame(grid_rows), "e2_width_q_grid.csv", LOG)
    save_json({"layers": LAYERS, "env": env_report()}, "e2_meta.json", LOG)


if __name__ == "__main__":
    main()
