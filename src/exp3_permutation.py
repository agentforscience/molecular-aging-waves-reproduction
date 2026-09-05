"""E3 -- The two permutation nulls, run on the same data in the same code.

Shen et al. 2024 report a permutation test under which their peaks disappear.
Carbonneau et al. 2026 report one under which they do not. The difference is a
single ordering step, and it decides whether the published crests have any
inferential support:

  Algorithm 1 (VALID, Carbonneau et al.)
      shuffle subject ages  ->  LOESS  ->  DE-SWAN
      Subject-level records are exchangeable under the null "no age
      relationship", so this null has the correct distribution. Every
      permutation goes through the full published pipeline, span selection
      included.

  Algorithm 2 (INVALID, Shen et al. Suppl. Fig. 4c)
      LOESS on the real data  ->  shuffle the age labels of the fitted grid
      points  ->  DE-SWAN
      LOESS output is *not* exchangeable: it is smooth by construction, so
      neighbouring grid values are near-duplicates. Shuffling after smoothing
      destroys exactly the structure that the test exploits, producing a null in
      which nothing is ever significant -- and therefore a null that the real
      data beats trivially, whether or not any signal exists.

We compute, for each layer and each null:
  * the null distribution of the DE-SWAN curve, window by window;
  * the null distribution of the *maximum* count over the 25 window centres,
    which gives family-wise control across the overlapping windows (the
    multiplicity dimension the published analysis never addressed);
  * the null distribution of the crest *location* -- where artifactual peaks
    land when there is provably no signal;
  * the ACHIEVED false discovery proportion of the published procedure, i.e.
    the fraction of tests it calls significant at nominal FDR 0.05 when every
    null hypothesis is true by construction.

Outputs
-------
results/tables/e3_null_curves.csv     mean/quantiles of the null curve per window
results/tables/e3_crest_null.csv      null crest locations (both algorithms)
results/tables/e3_pvalues.csv         permutation p-values + achieved FDR
results/e3_null_raw.npz               raw per-permutation curves
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd

from common import (env_report, get_logger, load_layer, save_json, save_table,
                    set_seed, zscore_rows, RESULTS)
from fastdeswan import DEFAULT_MIDPOINTS, count_significant, crest_age, run_deswan
from fastloess import DEFAULT_GRID, LoessEngine

LOG = get_logger("exp3_permutation")

# Permutation counts chosen so every layer finishes in the session; the
# transcriptome is 8,556 variables and the dominant cost.
N_PERM = {"plasma_transcriptome": 300}
N_PERM_DEFAULT = 1000
LAYERS = ["plasma_transcriptome", "plasma_proteomics", "plasma_cytokine",
          "clinical_test", "plasma_metabolomics_metabolite", "plasma_lipidomics"]
Q_THR = 0.05


def deswan_counts(M, coord, q_thr=Q_THR):
    res = run_deswan(M, coord, DEFAULT_MIDPOINTS, 10.0, test="wilcoxon")
    return count_significant(res["Q_bh"], q_thr), res


def run_layer(layer: str, rng) -> tuple[list, list, list, dict]:
    Y, ages, info, _ = load_layer(layer)
    Yz = zscore_rows(Y)
    p, n = Yz.shape
    n_perm = N_PERM.get(layer, N_PERM_DEFAULT)
    eng = LoessEngine(ages, grid=DEFAULT_GRID)

    # ------------------------------------------------ observed (published arm)
    Yhat = eng.smooth(Yz)
    obs_cnt, obs_res = deswan_counts(Yhat, eng.grid)
    obs_crest = crest_age(obs_cnt, DEFAULT_MIDPOINTS)
    LOG.info(f"[{layer}] observed: peak {obs_cnt.max()}/{p} at age {obs_crest}, "
             f"floor {obs_cnt.min()} | running {n_perm} permutations x 2 algorithms")

    null1 = np.zeros((n_perm, DEFAULT_MIDPOINTS.size), int)   # valid
    null2 = np.zeros((n_perm, DEFAULT_MIDPOINTS.size), int)   # invalid
    fdp1 = np.zeros(n_perm)                                   # achieved FDP, valid null

    t0 = time.time()
    Yhat_sorted_grid = Yhat            # smoothed real data, grid-ordered columns
    for b in range(n_perm):
        # --- Algorithm 1: permute ages, then smooth, then test ---------------
        Yp = eng.smooth_permuted(Yz, rng, reselect_spans=True)
        c1, _ = deswan_counts(Yp, eng.grid)
        null1[b] = c1
        # Under this null every hypothesis is true, so every rejection is false:
        # the realised false discovery proportion is the rejection rate itself.
        fdp1[b] = c1.sum() / (p * DEFAULT_MIDPOINTS.size)

        # --- Algorithm 2: smooth first, then permute the grid labels ---------
        gperm = rng.permutation(eng.grid.size)
        c2, _ = deswan_counts(Yhat_sorted_grid[:, gperm], eng.grid)
        null2[b] = c2

        if b == 0:
            LOG.info(f"[{layer}]   {time.time()-t0:.1f}s per permutation pair "
                     f"-> ~{(time.time()-t0)*n_perm/60:.1f} min total")

    LOG.info(f"[{layer}] permutations done in {(time.time()-t0)/60:.1f} min")

    # ------------------------------------------------------------- statistics
    max1, max2 = null1.max(axis=1), null2.max(axis=1)
    obs_max = obs_cnt.max()
    # Max-statistic permutation p-value = FWER-controlled test of "is there any
    # window where more molecules change than chance allows?"
    p_max_valid = (1 + np.sum(max1 >= obs_max)) / (1 + n_perm)
    p_max_invalid = (1 + np.sum(max2 >= obs_max)) / (1 + n_perm)

    crest1 = np.array([crest_age(c, DEFAULT_MIDPOINTS) if c.max() > 0 else np.nan
                       for c in null1])
    crest2 = np.array([crest_age(c, DEFAULT_MIDPOINTS) if c.max() > 0 else np.nan
                       for c in null2])

    curve_rows = []
    for j, mid in enumerate(DEFAULT_MIDPOINTS):
        pw_valid = (1 + np.sum(null1[:, j] >= obs_cnt[j])) / (1 + n_perm)
        pw_invalid = (1 + np.sum(null2[:, j] >= obs_cnt[j])) / (1 + n_perm)
        curve_rows.append(dict(
            layer=layer, midpoint=float(mid), n_variables=p,
            observed=int(obs_cnt[j]),
            null_valid_mean=float(null1[:, j].mean()),
            null_valid_q025=float(np.quantile(null1[:, j], 0.025)),
            null_valid_q975=float(np.quantile(null1[:, j], 0.975)),
            null_invalid_mean=float(null2[:, j].mean()),
            null_invalid_q975=float(np.quantile(null2[:, j], 0.975)),
            p_window_valid=float(pw_valid), p_window_invalid=float(pw_invalid),
            n_young=int(obs_res["n_young"][j]), n_old=int(obs_res["n_old"][j]),
        ))

    crest_rows = [dict(layer=layer, algorithm="valid_shuffle_before_loess",
                       permutation=int(b), crest=float(c), peak=int(m))
                  for b, (c, m) in enumerate(zip(crest1, max1))]
    crest_rows += [dict(layer=layer, algorithm="invalid_shuffle_after_loess",
                        permutation=int(b), crest=float(c), peak=int(m))
                   for b, (c, m) in enumerate(zip(crest2, max2))]

    pval_row = dict(
        layer=layer, n_variables=p, n_subjects=n, n_perm=n_perm,
        observed_peak=int(obs_max), observed_crest=obs_crest,
        observed_floor=int(obs_cnt.min()),
        null_valid_peak_mean=float(max1.mean()),
        null_valid_peak_q975=float(np.quantile(max1, 0.975)),
        null_invalid_peak_mean=float(max2.mean()),
        null_invalid_peak_q975=float(np.quantile(max2, 0.975)),
        p_max_valid=float(p_max_valid), p_max_invalid=float(p_max_invalid),
        achieved_fdp_mean=float(fdp1.mean()), achieved_fdp_max=float(fdp1.max()),
        nominal_fdr=Q_THR,
        null_crest_mode=float(pd.Series(crest1).mode().iloc[0]) if np.isfinite(crest1).any() else np.nan,
        null_crest_frac_40_48=float(np.mean((crest1 >= 40) & (crest1 <= 48))),
        null_crest_frac_56_64=float(np.mean((crest1 >= 56) & (crest1 <= 64))),
    )
    LOG.info(f"[{layer}] p(max | VALID null)   = {p_max_valid:.3f}   "
             f"[null peak mean {max1.mean():.0f} vs observed {obs_max}]")
    LOG.info(f"[{layer}] p(max | INVALID null) = {p_max_invalid:.3f}   "
             f"[null peak mean {max2.mean():.1f}]")
    LOG.info(f"[{layer}] achieved FDP at nominal FDR 0.05 = {fdp1.mean():.3f}  "
             f"({fdp1.mean()/Q_THR:.0f}x nominal)")

    raw = dict(null_valid=null1, null_invalid=null2, observed=obs_cnt, fdp=fdp1)
    return curve_rows, crest_rows, [pval_row], raw


def main() -> None:
    rng = set_seed()
    curves, crests, pvals, raws = [], [], [], {}
    for layer in LAYERS:
        try:
            c, cr, pv, raw = run_layer(layer, rng)
            curves += c
            crests += cr
            pvals += pv
            for k, v in raw.items():
                raws[f"{layer}__{k}"] = v
        except Exception as exc:
            LOG.error(f"[{layer}] FAILED: {exc!r}")

    save_table(pd.DataFrame(curves), "e3_null_curves.csv", LOG)
    save_table(pd.DataFrame(crests), "e3_crest_null.csv", LOG)
    save_table(pd.DataFrame(pvals), "e3_pvalues.csv", LOG)
    np.savez_compressed(f"{RESULTS}/e3_null_raw.npz", **raws)
    save_json({"n_perm": {**{"default": N_PERM_DEFAULT}, **N_PERM},
               "env": env_report()}, "e3_meta.json", LOG)
    LOG.info("\n" + pd.DataFrame(pvals).to_string(index=False))


if __name__ == "__main__":
    main()
