"""E4 -- Null simulations and age-distribution resampling: is the crest a
property of the biology, or of the estimator and the sampling design?

Four sub-experiments, each of which isolates a mechanism:

4a  PURE NULL AT THE REAL AGE DISTRIBUTION. Generate 1,000 molecules i.i.d.
    N(0,1), independent of age, at the observed iPOP ages. Run the *complete*
    published pipeline. If crests appear at ~44 and ~60 in data that contain no
    age relationship whatsoever, the crest ages are a property of the pipeline
    plus the age distribution, not of aging. Repeated over many simulations to
    get the null crest *distribution* rather than one anecdote.

4b  CURVE SHAPE VS AGE DISTRIBUTION. The same null generator under four age
    distributions (real iPOP, uniform, normal, bimodal), holding n fixed. Any
    difference in curve shape is attributable to the age distribution alone.

4c  RESAMPLING THE REAL DATA. Take the real iPOP transcriptome/proteome and
    resample subjects to (i) a flattened, near-uniform age distribution and
    (ii) a control that preserves the observed distribution -- both at the same
    n, so only the age distribution differs. If the published crest moves under
    (i) but not (ii), the crest location is set by the sampling design.

4d  DE-SWAN ON LINEAR DATA (Artifact 3, no LOESS involved). Generate strictly
    linear age trends, run DE-SWAN *without* smoothing, and show the resulting
    curve alongside the per-window group sizes. Variants: homoskedastic,
    heteroskedastic (increasing/decreasing variance), and outlier-contaminated.
    A monotone curve in the counts would mean the method is well behaved; peaks
    mean a linear signal has been rendered as a "wave".

Outputs
-------
results/tables/e4a_null_crests.csv    null crest locations under the real ages
results/tables/e4a_null_curves.csv    mean null curve per window
results/tables/e4b_age_dists.csv      curve shape by age distribution
results/tables/e4c_resampling.csv     crest movement under resampling
results/tables/e4d_linear.csv         DE-SWAN on linear data, with group sizes
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from common import (env_report, get_logger, load_layer, save_json, save_table,
                    set_seed, zscore_rows)
from fastdeswan import DEFAULT_MIDPOINTS, count_significant, crest_age, run_deswan
from fastloess import DEFAULT_GRID, LoessEngine

LOG = get_logger("exp4_nullsim")

N_MOL = 1000        # molecules per null simulation (Carbonneau et al. use 1,000)
N_SIM_A = 300       # simulations for the null crest distribution
N_SIM_B = 100       # simulations per age distribution
N_RESAMPLE = 200    # resamples per resampling scheme
Q_THR = 0.05


def pipeline_counts(Y, ages, eng=None, use_loess=True):
    """Run the published pipeline (or DE-SWAN alone) and return the count curve."""
    if use_loess:
        eng = eng or LoessEngine(ages, grid=DEFAULT_GRID)
        M, coord = eng.smooth(zscore_rows(Y)), eng.grid
    else:
        M, coord = zscore_rows(Y), ages
    res = run_deswan(M, coord, DEFAULT_MIDPOINTS, 10.0, test="wilcoxon")
    return count_significant(res["Q_bh"], Q_THR), res


# --------------------------------------------------------------------- 4a
def exp_4a(rng, ages_real):
    """Pure null (no age relationship) at the real iPOP age distribution."""
    eng = LoessEngine(ages_real, grid=DEFAULT_GRID)
    crests, curves, peaks = [], np.zeros((N_SIM_A, DEFAULT_MIDPOINTS.size)), []
    for b in range(N_SIM_A):
        Y = rng.standard_normal((N_MOL, ages_real.size))
        cnt, _ = pipeline_counts(Y, ages_real, eng)
        curves[b] = cnt
        peaks.append(cnt.max())
        crests.append(crest_age(cnt, DEFAULT_MIDPOINTS) if cnt.max() > 0 else np.nan)
    crests = np.array(crests, float)

    LOG.info(f"[4a] pure null, real iPOP ages: mean peak "
             f"{np.mean(peaks):.0f}/{N_MOL} molecules significant at FDR {Q_THR}")
    LOG.info(f"[4a] null crest distribution: median {np.nanmedian(crests):.1f}, "
             f"{100*np.mean((crests>=40)&(crests<=48)):.0f}% in 40-48, "
             f"{100*np.mean((crests>=56)&(crests<=64)):.0f}% in 56-64")

    crest_df = pd.DataFrame(dict(simulation=np.arange(N_SIM_A), crest=crests,
                                 peak=peaks))
    curve_df = pd.DataFrame(dict(
        midpoint=DEFAULT_MIDPOINTS,
        mean_significant=curves.mean(axis=0),
        q025=np.quantile(curves, 0.025, axis=0),
        q975=np.quantile(curves, 0.975, axis=0),
        frac_crest_here=[np.mean(crests == m) for m in DEFAULT_MIDPOINTS],
        n_molecules=N_MOL,
    ))
    return crest_df, curve_df


# --------------------------------------------------------------------- 4b
def age_distributions(rng, ages_real, n):
    lo, hi = 26.0, 75.0
    return {
        "ipop_real": ages_real.copy(),
        "uniform": rng.uniform(lo, hi, n),
        "normal": np.clip(rng.normal(50, 10, n), lo, hi),
        "bimodal": np.clip(np.where(rng.random(n) < 0.5,
                                    rng.normal(37, 5, n), rng.normal(64, 5, n)), lo, hi),
    }


def exp_4b(rng, ages_real):
    """Same null generator, four age distributions, n held fixed."""
    n = ages_real.size
    rows = []
    for name, ages in age_distributions(rng, ages_real, n).items():
        eng = LoessEngine(ages, grid=DEFAULT_GRID)
        curves = np.zeros((N_SIM_B, DEFAULT_MIDPOINTS.size))
        crests = []
        for b in range(N_SIM_B):
            Y = rng.standard_normal((N_MOL, n))
            cnt, res = pipeline_counts(Y, ages, eng)
            curves[b] = cnt
            crests.append(crest_age(cnt, DEFAULT_MIDPOINTS) if cnt.max() > 0 else np.nan)
        # group sizes of the *unsmoothed* windows, i.e. the sampling design
        _, res_raw = pipeline_counts(np.zeros((1, n)) + rng.standard_normal((1, n)),
                                     ages, use_loess=False)
        for j, mid in enumerate(DEFAULT_MIDPOINTS):
            rows.append(dict(age_dist=name, midpoint=float(mid),
                             mean_significant=float(curves[:, j].mean()),
                             q025=float(np.quantile(curves[:, j], 0.025)),
                             q975=float(np.quantile(curves[:, j], 0.975)),
                             n_young_raw=int(res_raw["n_young"][j]),
                             n_old_raw=int(res_raw["n_old"][j]),
                             crest_frac=float(np.mean(np.array(crests) == mid)),
                             n_molecules=N_MOL, n_subjects=n))
        LOG.info(f"[4b] {name:10s}: null crest median "
                 f"{np.nanmedian(crests):.1f}, mean peak {curves.max(axis=1).mean():.0f}")
    return pd.DataFrame(rows)


# --------------------------------------------------------------------- 4c
def resample_indices(rng, ages, n_out, target: str):
    """Draw n_out subjects without replacement, either flattening the age
    distribution or preserving it (control)."""
    if target == "preserve":
        return rng.choice(ages.size, n_out, replace=False)
    # Flatten: weight inversely to a kernel density estimate of the observed ages
    from scipy.stats import gaussian_kde
    dens = gaussian_kde(ages)(ages)
    w = 1.0 / np.maximum(dens, 1e-9)
    w /= w.sum()
    return rng.choice(ages.size, n_out, replace=False, p=w)


def exp_4c(rng, layers=("plasma_transcriptome", "plasma_proteomics")):
    rows = []
    for layer in layers:
        Y, ages, _, _ = load_layer(layer)
        n_out = int(0.7 * ages.size)
        for target in ("preserve", "flatten"):
            crests, peaks, cvs = [], [], []
            for b in range(N_RESAMPLE):
                idx = resample_indices(rng, ages, n_out, target)
                a = ages[idx]
                if a.max() - a.min() < 40:            # need coverage of 40-65
                    continue
                cnt, res = pipeline_counts(Y[:, idx], a)
                crests.append(crest_age(cnt, DEFAULT_MIDPOINTS) if cnt.max() > 0 else np.nan)
                peaks.append(cnt.max())
                # coefficient of variation of raw per-window group sizes: a
                # scalar summary of how uneven statistical power is along the axis
                _, rr = pipeline_counts(Y[:, idx], a, use_loess=False)
                ny = rr["n_young"].astype(float)
                cvs.append(ny.std() / max(ny.mean(), 1e-9))
            crests = np.array(crests, float)
            rows.append(dict(
                layer=layer, target=target, n_resamples=len(crests), n_subjects=n_out,
                crest_median=float(np.nanmedian(crests)),
                crest_iqr=float(np.nanpercentile(crests, 75) - np.nanpercentile(crests, 25)),
                crest_q025=float(np.nanpercentile(crests, 2.5)),
                crest_q975=float(np.nanpercentile(crests, 97.5)),
                crest_sd=float(np.nanstd(crests)),
                frac_crest_44pm2=float(np.mean(np.abs(crests - 44) <= 2)),
                frac_crest_60pm2=float(np.mean(np.abs(crests - 60) <= 2)),
                mean_peak=float(np.mean(peaks)),
                mean_groupsize_cv=float(np.mean(cvs)),
            ))
            LOG.info(f"[4c] {layer} / {target}: crest median "
                     f"{rows[-1]['crest_median']:.1f} "
                     f"[95% {rows[-1]['crest_q025']:.0f}-{rows[-1]['crest_q975']:.0f}], "
                     f"P(crest within 2y of 44) = {rows[-1]['frac_crest_44pm2']:.2f}")
    return pd.DataFrame(rows)


# --------------------------------------------------------------------- 4d
def exp_4d(rng, ages_real):
    """DE-SWAN (no LOESS) applied to strictly LINEAR age trends.

    Any structure in these curves is manufactured by age-varying statistical
    power, since the underlying truth is a straight line at every age.
    """
    n = ages_real.size
    z = (ages_real - ages_real.mean()) / ages_real.std()
    slope = 0.45          # effect size chosen to put counts in a detectable range
    rows = []
    scenarios = {
        "homoskedastic": lambda: np.ones(n),
        "variance_increasing": lambda: 0.5 + 1.5 * (z - z.min()) / (z.max() - z.min()),
        "variance_decreasing": lambda: 2.0 - 1.5 * (z - z.min()) / (z.max() - z.min()),
        "outlier_cluster_55": lambda: np.where(np.abs(ages_real - 55) < 3, 4.0, 1.0),
    }
    for name, sdfn in scenarios.items():
        sd = sdfn()
        curves = np.zeros((N_SIM_B, DEFAULT_MIDPOINTS.size))
        for b in range(N_SIM_B):
            Y = slope * z[None, :] + rng.standard_normal((N_MOL, n)) * sd[None, :]
            cnt, res = pipeline_counts(Y, ages_real, use_loess=False)
            curves[b] = cnt
        _, res = pipeline_counts(Y, ages_real, use_loess=False)
        mean_curve = curves.mean(axis=0)
        for j, mid in enumerate(DEFAULT_MIDPOINTS):
            rows.append(dict(scenario=name, midpoint=float(mid),
                             mean_significant=float(mean_curve[j]),
                             q025=float(np.quantile(curves[:, j], 0.025)),
                             q975=float(np.quantile(curves[:, j], 0.975)),
                             n_young=int(res["n_young"][j]), n_old=int(res["n_old"][j]),
                             n_molecules=N_MOL))
        pk = DEFAULT_MIDPOINTS[int(np.argmax(mean_curve))]
        # Correlation between the count curve and the harmonic-mean group size,
        # the standard sample-size proxy for two-sample power.
        hm = 2 / (1 / np.maximum(res["n_young"], 1) + 1 / np.maximum(res["n_old"], 1))
        r = np.corrcoef(mean_curve, hm)[0, 1]
        rows[-1]  # noqa
        LOG.info(f"[4d] {name:22s}: apparent crest at {pk:.0f} on strictly LINEAR data; "
                 f"corr(counts, harmonic-mean n) = {r:.2f}")
    return pd.DataFrame(rows)


def main() -> None:
    rng = set_seed()
    _, ages_real, _, _ = load_layer("plasma_transcriptome")

    crest_df, curve_df = exp_4a(rng, ages_real)
    save_table(crest_df, "e4a_null_crests.csv", LOG)
    save_table(curve_df, "e4a_null_curves.csv", LOG)

    save_table(exp_4b(rng, ages_real), "e4b_age_dists.csv", LOG)
    save_table(exp_4c(rng), "e4c_resampling.csv", LOG)
    save_table(exp_4d(rng, ages_real), "e4d_linear.csv", LOG)

    save_json({"n_molecules": N_MOL, "n_sim_a": N_SIM_A, "n_sim_b": N_SIM_B,
               "n_resample": N_RESAMPLE, "q_threshold": Q_THR,
               "env": env_report()}, "e4_meta.json", LOG)


if __name__ == "__main__":
    main()
