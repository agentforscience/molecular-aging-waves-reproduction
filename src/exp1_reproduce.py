"""E1 -- Reproduce the published LOESS+DE-SWAN pipeline across all 10 iPOP layers,
and run DE-SWAN on the same raw per-subject data without LOESS.

This is the foundation every later experiment compares against. Two arms per layer:

  A. LOESS+DE-SWAN (published)   z-score -> per-variable CV-span LOESS onto the
     half-year grid 26..75 -> Mann-Whitney between window halves -> BH within
     each window centre -> count q < 0.05. Windows: 20 years wide, centres 40..64.
  B. DE-SWAN only (no LOESS)     the identical test applied to the raw
     per-subject values.

Both arms report per-window (n_young, n_old) -- the quantity Carbonneau et al.
identify as the mechanism behind Artifact 3, and which no published DE-SWAN
curve has ever been shown alongside.

Outputs
-------
results/tables/e1_curves.csv     per layer x arm x window: counts, group sizes
results/tables/e1_summary.csv    per layer x arm: crest age, peak count, totals
results/e1_meta.json             environment + configuration
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd

from common import (LAYERS, env_report, get_logger, load_layer, save_json,
                    save_table, set_seed, zscore_rows)
from fastdeswan import DEFAULT_MIDPOINTS, count_significant, crest_age, run_deswan
from fastloess import DEFAULT_GRID, LoessEngine

LOG = get_logger("exp1_reproduce")
Q_THRESHOLD = 0.05
HALF_WIDTH = 10.0  # 20-year window, as published


def run_layer(layer: str) -> tuple[pd.DataFrame, list[dict]]:
    t0 = time.time()
    Y, ages, info, var_ids = load_layer(layer)
    Yz = zscore_rows(Y)
    p, n = Yz.shape
    LOG.info(f"[{layer}] {p} variables x {n} subjects, log1p={info.attrs['log1p']}")

    rows, summaries = [], []

    # ---- Arm A: the published pipeline -------------------------------------
    eng = LoessEngine(ages, grid=DEFAULT_GRID)
    spans, _ = eng.select_spans(Yz)
    Yhat = eng.smooth(Yz, spans=spans)
    resA = run_deswan(Yhat, eng.grid, DEFAULT_MIDPOINTS, HALF_WIDTH,
                      test="wilcoxon", scheme="fixed_width", with_effect=True)
    cA = count_significant(resA["Q_bh"], Q_THRESHOLD)

    # ---- Arm B: same test, raw per-subject data, no smoothing --------------
    resB = run_deswan(Yz, ages, DEFAULT_MIDPOINTS, HALF_WIDTH,
                      test="wilcoxon", scheme="fixed_width", with_effect=True)
    cB = count_significant(resB["Q_bh"], Q_THRESHOLD)

    for arm, res, cnt in (("loess_deswan", resA, cA), ("deswan_only", resB, cB)):
        for j, mid in enumerate(res["midpoints"]):
            rows.append(dict(
                layer=layer, arm=arm, midpoint=float(mid),
                n_significant=int(cnt[j]),
                frac_significant=float(cnt[j] / p),
                n_young=int(res["n_young"][j]), n_old=int(res["n_old"][j]),
                min_p=float(np.nanmin(res["P"][:, j])) if np.isfinite(res["P"][:, j]).any() else np.nan,
                median_abs_d=float(np.nanmedian(np.abs(res["D"][:, j]))),
                max_abs_d=float(np.nanmax(np.abs(res["D"][:, j]))) if np.isfinite(res["D"][:, j]).any() else np.nan,
            ))
        summaries.append(dict(
            layer=layer, arm=arm, n_variables=p, n_subjects=n,
            crest_age=crest_age(cnt, res["midpoints"]),
            peak_n_significant=int(np.max(cnt)),
            peak_frac=float(np.max(cnt) / p),
            min_n_significant=int(np.min(cnt)),
            total_significant_tests=int(np.sum(cnt)),
            n_windows_with_any=int(np.sum(cnt > 0)),
            median_span=float(np.median(spans)) if arm == "loess_deswan" else np.nan,
        ))
    LOG.info(f"[{layer}] LOESS+DE-SWAN peak {cA.max()}/{p} at age "
             f"{crest_age(cA, resA['midpoints']):.0f} | DE-SWAN-only peak {cB.max()}/{p} "
             f"(total {cB.sum()} significant tests over 25 windows) | {time.time()-t0:.1f}s")
    return pd.DataFrame(rows), summaries


def main() -> None:
    set_seed()
    all_rows, all_sum = [], []
    for layer in LAYERS:
        try:
            df, s = run_layer(layer)
            all_rows.append(df)
            all_sum.extend(s)
        except Exception as exc:                       # keep going; record the failure
            LOG.error(f"[{layer}] FAILED: {exc!r}")
            all_sum.append(dict(layer=layer, arm="ERROR", error=repr(exc)))

    curves = pd.concat(all_rows, ignore_index=True)
    summary = pd.DataFrame(all_sum)
    save_table(curves, "e1_curves.csv", LOG)
    save_table(summary, "e1_summary.csv", LOG)
    save_json({"config": dict(q_threshold=Q_THRESHOLD, half_width=HALF_WIDTH,
                              midpoints=DEFAULT_MIDPOINTS.tolist(),
                              grid=[float(DEFAULT_GRID[0]), float(DEFAULT_GRID[-1]),
                                    len(DEFAULT_GRID)]),
               "env": env_report()}, "e1_meta.json", LOG)

    LOG.info("\n" + summary.to_string(index=False))


if __name__ == "__main__":
    main()
