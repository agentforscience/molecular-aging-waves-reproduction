# Do the claimed "waves" of molecular aging survive rigorous statistical testing?

A reproduction and extension of the LOESS + DE-SWAN artifact analysis, testing whether the
widely-cited discrete aging transitions at **~44** and **~60** years are biological signal or
statistical artifacts of the estimator used to find them.

All ten omics layers of the iPOP cohort (Shen et al. 2024, *Nature Aging*) are re-analysed from
recovered raw data, plus an independent, six-times-larger methylation cohort (GSE40279, n = 656).

**Full write-up: [REPORT.md](REPORT.md).** Planning and direction budget: [planning.md](planning.md).
Literature synthesis: [literature_review.md](literature_review.md).

## Key findings

- **The published pipeline calls a majority of molecules significant at *every* age, in all 10
  omics layers.** LOESS + DE-SWAN declares 58–94 % of variables differentially expressed at its
  peak and never fewer than 12 % at its floor. The identical test on the same raw per-subject data
  without LOESS finds **8 significant tests out of 1,215,175** across all layers.
- **The crests are not significant under a valid permutation null.** Shuffling ages *before*
  smoothing (the exchangeable null) gives p = 0.37–0.66 in five of six layers — the observed peak
  is no larger than chance, and in four layers it is *below* the null mean. Shuffling *after*
  smoothing, as published, gives p ≤ 0.003 in every layer. The entire published significance rests
  on that ordering.
- **The procedure's achieved error rate is ~13× its nominal rate.** Measured on data where every
  null hypothesis is true by construction, the false discovery proportion is 0.60–0.69 against a
  nominal FDR of 0.05.
- **Both published transition ages are reproduced by data containing no age relationship at all.**
  1,000 pure-noise molecules at the real iPOP ages yield a DE-SWAN curve with an interior local
  maximum at **61** and, within the published 40–65 plotting window, an apparent crest in the low
  40s. Changing only the age distribution moves the crest to 57 (uniform), 62 (bimodal) or 64
  (normal).
- **The published robustness check cannot distinguish signal from artifact.** Across the same
  16-cell bucket-width × q-threshold grid used to support the ~44 crest, 79–100 % of pure-noise
  simulations produce a crest at least as stable as the real data.
- **DE-SWAN turns strictly linear data into waves even without LOESS**, with apparent crests at
  54–60 and counts correlating r = 0.60–0.82 with per-window sample size.
- **Under corrected inference, iPOP contains no discrete transitions at all.** Permutation-
  calibrated segmented regression, spline nonlinearity tests and a power-equalised DE-SWAN with
  family-wise control across window centres find no surviving transition in any layer (FWER
  p ≥ 0.25). The one apparent exception — 35/66 cytokines "transitioning" at 41 — is traced to a
  **single subject aged 25.9**: dropping that one subject leaves 0 significant analytes, while
  dropping a random subject leaves a median of 7.

## Reproducing

```bash
uv venv && source .venv/bin/activate && uv sync

python src/validate_fastloess.py        # validates the LOESS engine (must pass first)
python src/exp1_reproduce.py            # ~11 s   all 10 layers, published vs no-LOESS
python src/exp2_ablation.py             # ~40 s   smoothing/densification decomposition
python src/exp3_permutation.py          # ~25 min valid vs invalid permutation nulls
python src/exp4_nullsim.py              # ~13 min null simulations, age-distribution resampling
python src/exp5_corrected.py            # ~15 min corrected inference on iPOP
python src/exp6_gse40279.py             # ~90 min independent methylation cohort
python src/exp7_robustness.py           # ~2 min  outlier diagnostics, RINT battery
python src/exp7b_cytokine_followup.py   # ~1 min  the one positive finding, attacked
python src/exp8_robustness_argument.py  # ~1 min  is the published robustness check informative?

python src/analyse_results.py           # headline numbers + hypothesis verdicts
python src/make_figures.py              # all figures
```

Everything is seeded (`SEED = 42`, set in `src/common.py`). No GPU is required; the whole pipeline
is CPU-bound linear algebra and runs in roughly 2.5 hours on 32 cores.

## Layout

```
src/
  fastloess.py       vectorised R-compatible LOESS (5,000x the reference; exact to 1e-14)
  fastdeswan.py      DE-SWAN: fixed-width / quantile / power-equalised windows, both test statistics
  corrected.py       the replacement estimator: spline nonlinearity, segmented sup-F, permutation nulls
  common.py          data loading, seeding, logging, I/O
  exp1..exp8         the experiments, one file each
  analyse_results.py synthesis into headline numbers and verdicts
  make_figures.py    figures
  validate_fastloess.py  numerical validation against the reference and the authors' R output

results/tables/      every result as CSV
results/             headline_numbers.json, per-experiment metadata
figures/             fig1..fig7 PNG
logs/                per-experiment run logs
datasets/            iPOP (10 layers, recovered) + GSE40279; see datasets/README.md
code/                cloned reference repositories (Carbonneau et al., Shen et al., DEswan)
```

## What this adds to the existing critique

Carbonneau et al. (2026, bioRxiv — with the original waves authors as co-authors) established that
LOESS + DE-SWAN can manufacture waves. They analysed one cohort and mostly one modality, offered no
replacement estimator, and were explicit that they did not rule out real transitions detectable
with other cohorts and models. This work extends that in four directions: all ten iPOP layers
rather than one; a quantitative decomposition (achieved FDR, smoothing vs densification, crest
confidence intervals) rather than a qualitative demonstration; a corrected inferential battery
delivered rather than recommended; and an independent, well-powered cohort to separate "the method
is broken" from "the cohort is underpowered".
