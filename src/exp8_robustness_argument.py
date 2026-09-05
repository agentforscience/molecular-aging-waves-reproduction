"""E8 -- Does the published robustness analysis distinguish signal from artifact?

Shen et al. support the ~44 crest by showing it is stable across bucket widths
{15, 20, 25, 30} and q-thresholds {1e-4, 1e-3, 1e-2, 0.05}. Our E2 reproduces
that stability: on the real iPOP transcriptome the crest sits at 43-46 in all 16
cells of the grid.

That argument is only evidence for biology if an ARTIFACTUAL crest would be
unstable across the same grid. Nobody has checked. This experiment does: it runs
the identical 16-cell robustness grid on data generated with no age relationship
at all, at the real iPOP ages, and asks whether the artifact is any less "robust"
than the published result.

If the null crest is equally stable, then robustness to bucket width and
q-threshold carries no evidential weight for or against a biological transition
-- because both hypotheses predict it. That would make the published robustness
analysis uninformative rather than merely insufficient.

Also included: the same grid on the real data and on the null, side by side, and
a direct comparison of the crest-location spread.

Outputs
-------
results/tables/e8_robustness_grid.csv
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from common import (env_report, get_logger, load_layer, save_json, save_table,
                    set_seed, zscore_rows)
from fastdeswan import DEFAULT_MIDPOINTS, count_significant, crest_age, run_deswan
from fastloess import DEFAULT_GRID, LoessEngine

LOG = get_logger("exp8_robustness_argument")

BUCKETS = (15, 20, 25, 30)
Q_THRESHOLDS = (1e-4, 1e-3, 1e-2, 0.05)
N_SIM = 100
N_MOL = 1000
LAYERS = ["plasma_transcriptome", "plasma_proteomics"]


def grid_crests(Yhat, grid):
    """Crest age in each cell of the published bucket x q-threshold grid."""
    out = {}
    for bucket in BUCKETS:
        res = run_deswan(Yhat, grid, DEFAULT_MIDPOINTS, bucket / 2.0, test="wilcoxon")
        for q in Q_THRESHOLDS:
            cnt = count_significant(res["Q_bh"], q)
            out[(bucket, q)] = crest_age(cnt, DEFAULT_MIDPOINTS) if cnt.max() > 0 else np.nan
    return out


def main() -> None:
    rng = set_seed()
    rows = []

    for layer in LAYERS:
        Y, ages, _, _ = load_layer(layer)
        Yz = zscore_rows(Y)
        eng = LoessEngine(ages, grid=DEFAULT_GRID)

        # ---- the real data, as in the published robustness analysis --------
        real = grid_crests(eng.smooth(Yz), eng.grid)
        rv = np.array([v for v in real.values() if np.isfinite(v)])
        LOG.info(f"[{layer}] REAL data: crest {np.nanmedian(rv):.0f} across the 16-cell "
                 f"grid, range {rv.min():.0f}-{rv.max():.0f}, sd {rv.std():.1f}")

        # ---- the same grid on data with NO age relationship ----------------
        null_cells = {k: [] for k in real}
        for b in range(N_SIM):
            Yn = rng.standard_normal((N_MOL, ages.size))
            g = grid_crests(eng.smooth(Yn), eng.grid)
            for k, v in g.items():
                null_cells[k].append(v)

        for k in real:
            nv = np.array(null_cells[k], float)
            rows.append(dict(
                layer=layer, bucket=k[0], q_threshold=k[1],
                real_crest=real[k],
                null_crest_median=float(np.nanmedian(nv)),
                null_crest_q025=float(np.nanpercentile(nv, 2.5)),
                null_crest_q975=float(np.nanpercentile(nv, 97.5)),
                null_frac_within2_of_real=float(np.nanmean(np.abs(nv - real[k]) <= 2)),
            ))

        # Spread of the crest across the grid: real vs one null realisation.
        per_sim_spread = []
        for b in range(N_SIM):
            v = np.array([null_cells[k][b] for k in real], float)
            v = v[np.isfinite(v)]
            if v.size:
                per_sim_spread.append(v.max() - v.min())
        real_spread = float(rv.max() - rv.min())
        p_spread = float(np.mean(np.array(per_sim_spread) <= real_spread))
        LOG.info(f"[{layer}] NULL data: crest spread across the same 16-cell grid is "
                 f"{np.median(per_sim_spread):.0f} years (median over {N_SIM} "
                 f"simulations) vs {real_spread:.0f} years for the real data; "
                 f"fraction of null simulations at least as 'robust' as the real "
                 f"data: {p_spread:.2f}")
        rows.append(dict(layer=layer, bucket=-1, q_threshold=-1.0,
                         real_crest=np.nan,
                         null_crest_median=float(np.median(per_sim_spread)),
                         null_crest_q025=np.nan, null_crest_q975=np.nan,
                         null_frac_within2_of_real=p_spread))

    save_table(pd.DataFrame(rows), "e8_robustness_grid.csv", LOG)
    save_json({"buckets": list(BUCKETS), "q_thresholds": list(Q_THRESHOLDS),
               "n_sim": N_SIM, "n_molecules": N_MOL, "env": env_report()},
              "e8_meta.json", LOG)


if __name__ == "__main__":
    main()
