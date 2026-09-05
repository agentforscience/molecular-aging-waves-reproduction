# Planning — Direction Budget

## Motivation & Novelty Assessment

### Why This Research Matters

A single methodological pipeline — LOESS smoothing followed by DE-SWAN sliding-window
differential expression — underwrites a claim that has reorganised how a field talks about
aging: that human molecular aging proceeds in **discrete waves**, notably at ages **~44 and
~60** (Shen et al. 2024, *Nature Aging*, ~400 citations in under two years) and **34/60/78**
(Lehallier et al. 2019, *Nat Med*, ~830 citations). Carbonneau et al. counted **>100 papers**
(2019–2026) inheriting LOESS, DE-SWAN, or both. A 2026 *Nature Reviews Genetics* review now
treats nonlinear aging as established, and "midlife critical window" framings are entering
clinical and public discourse. If the crests are artifacts of the estimator, then a large,
still-growing literature is built on a statistical mirage — and, worse, on a *reusable* mirage
that any group can regenerate from any cohort. The beneficiaries of settling this are aging
biologists deciding where to spend intervention studies, methodologists who need a defensible
replacement, and reviewers who currently have no calibrated benchmark to judge these curves
against.

### Gap in Existing Work

The June 2026 bioRxiv preprint of Carbonneau, Shutta, Miller, Snyder, Shen & Quackenbush
(with Shen and Snyder — first and last authors of the waves paper — as co-authors) already
establishes the *qualitative* core of our hypothesis: LOESS output is not exchangeable, feeding
it to a rank test is invalid, and DE-SWAN's significant-molecule count confounds effect size
with age-varying statistical power. Our literature review confirmed this is not an outside
attack but a jointly-endorsed correction. That materially changes what is worth doing here.
Four things the preprint explicitly does **not** do:

1. **It is a diagnosis without a replacement.** It recommends (permute before smoothing;
   characterise on nulls; prefer effect sizes to counts) but delivers no corrected estimator and
   no calibrated test for age-localised change.
2. **It never asks whether the biology survives.** It is explicit that it does "not rule out
   nonlinear molecular aging or age-associated transitions that may be detectable using other
   cohorts and statistical models." Our hypothesis is *stronger* than their conclusion, and the
   gap between the two is exactly the untested territory.
3. **One cohort, essentially one modality.** iPOP transcriptomics, with proteomics in a
   supplementary notebook. Nothing is known about the other 8 omics layers, or about what
   happens in a cohort with real power.
4. **No quantitative decomposition.** Nobody has partitioned a published crest into artifact
   vs. residual signal, reported the *achieved* FDR of the published procedure, or given a crest
   location with a confidence interval instead of a point estimate.

### Our Novel Contribution

We test the strong form of the hypothesis that the preprint left open, and we do it with a
falsifiable design that can lose. Three contributions:

- **A quantitative decomposition of the artifact** across all 10 iPOP omics layers (not one),
  including the achieved-vs-nominal FDR of the published procedure, the crest location's
  dependence on the age distribution alone, and both permutation nulls run side by side.
- **A corrected inferential battery** — the replacement the preprint calls for but does not
  supply: penalised-spline nonlinearity tests, permutation-calibrated segmented regression with
  breakpoint CIs, and a **power-equalised DE-SWAN** that holds per-window group sizes fixed by
  construction, with family-wise control across the overlapping window centres via a
  max-statistic permutation null.
- **An independent, well-powered replication** in GSE40279 (n=656, ages 19–101, ~6× iPOP and
  far more uniform in age), which is what separates *"the method is broken"* from *"iPOP is
  underpowered."*

We also keep three claims separate that the literature routinely conflates: *aging is
nonlinear*; *there exist discrete transition ages*; *those ages are 44 and 60*.

### Experiment Justification

- **E1 — Reproduce the published pipeline on all 10 layers, with per-window (n_young, n_old)
  reported.** Needed because everything downstream is a comparison against it, and because no
  published DE-SWAN curve has ever been shown alongside the group sizes that are the alleged
  mechanism. Establishes the artifact is present in *our* independent reconstruction of the data.
- **E2 — Pipeline ablation (LOESS on/off × grid on/off × test × bucket width × q-threshold).**
  Isolates *which* step manufactures the crest. Without this, "the pipeline is broken" is not
  actionable; with it, we can name the responsible component.
- **E3 — Valid (shuffle-before-LOESS) vs. invalid (shuffle-after-LOESS) permutation nulls.**
  This single ordering step is the entire disagreement between Shen et al.'s reported
  permutation check and Carbonneau et al.'s. It must be run on the same data, in the same code,
  to be conclusive rather than a he-said-she-said.
- **E4 — Null simulations with the real iPOP age distribution + age-distribution resampling.**
  Tests the mechanism directly: if crests appear at ~44 and ~60 in data generated with *no age
  relationship whatsoever*, and if the crest moves when only the age distribution is changed,
  the crest is a property of the estimator and the sampling design, not of biology.
- **E5 — Corrected inference on iPOP (splines, segmented regression, power-equalised DE-SWAN,
  two-dimensional multiplicity control).** Answers the question the preprint leaves open: given
  the pipeline is broken, does *any* age-localised change survive proper inference, and where?
  This is the only experiment that can *refute* our hypothesis.
- **E6 — The same corrected battery on GSE40279.** Distinguishes a broken method from a small
  cohort, and provides the independent test of "transitions at ages other than 44/60." A
  transition surviving here with a CI covering 44 or 60 would mean the method was wrong but the
  biological conclusion right — an outcome we commit in advance to reporting as such.

---

**Hypothesis.** The majority of reported discrete aging transitions in omics data are
statistical artifacts of the LOESS smoothing + DE-SWAN sliding-window pipeline rather than
biological signal. Corrected methods with proper multiple-testing control and cross-validated
model selection will find either no transitions or transitions at ages other than ~44/~60.

**Position after Phase 1 literature review.** A June 2026 bioRxiv preprint (Carbonneau,
Shutta, Miller, Snyder, Shen, Quackenbush — `papers/carbonneau2026_loess_deswan_artifactual_waves.pdf`)
already establishes the *qualitative* core of this hypothesis, and does so with the original
authors as co-authors. This **materially changes what is worth doing**: the novelty budget must
go to what that preprint explicitly leaves open, not to re-deriving what it settled. Its own
stated limits are the openings:

- It does *not* claim the waves are absent — only that the evidence does not support them
  ("do not rule out nonlinear molecular aging ... detectable using other cohorts and statistical models").
- It analyses **one cohort** (iPOP), mostly **transcriptomics** (proteomics in a supplementary notebook).
- It **diagnoses** but does not **replace**: it offers recommendations, not a corrected estimator
  or a calibrated test for age-localised change.
- It quantifies the artifact but never asks *what fraction* of the published crest is artifact
  vs. signal, nor where a properly-powered analysis would place a transition if one exists.

I verified the central claim end-to-end during Phase 1 before planning around it
(`code/deswan_py/smoke_test_proteomics.py`): on the reconstructed iPOP proteomics
cross-section, LOESS+DE-SWAN calls 174–245 of 302 proteins significant at FDR 0.05 in
*every* window with a maximum at age 44, while DE-SWAN on the same raw per-subject data
calls **1** test significant across all 25 windows. The artifact is real and reproducible here.

---

## Directions considered

Scored 1–5 on four axes. **Evi** = evidentiary support from the literature that the direction is
well-posed; **Rel** = relevance to the stated hypothesis; **Gain** = expected information gain
*given that Carbonneau et al. already exists*; **Feas** = implementation feasibility with the
resources actually in hand (no R, ~100 subjects, single machine).

| # | Direction | Evi | Rel | Gain | Feas | Total |
|---|-----------|:---:|:---:|:---:|:---:|:-----:|
| **D1** | **Reproduce + quantitatively decompose the artifact** — rerun the published LOESS+DE-SWAN across all 10 iPOP omics layers; ablate the pipeline (LOESS on/off, grid interpolation on/off, Wilcoxon vs linear test, bucket width, q-threshold); implement the valid *shuffle-before-LOESS* permutation null vs. the published *shuffle-after-LOESS* null; report FDR actually achieved vs. nominal. | 5 | 5 | 4 | 5 | **19** |
| **D2** | **Corrected inference: does *any* age-localised change survive?** — replace the pipeline with pre-registered alternatives on the same raw data: (a) GAM/penalised spline with cross-validated smoothness and a proper test of nonlinearity vs. linear; (b) segmented/changepoint regression with permutation-calibrated breakpoint CIs; (c) power-equalised DE-SWAN (equal *n* per half-window by design, not by quantile centres); all with BH/BY control and permutation-calibrated global nulls. Report transition ages with CIs, or a calibrated null result. | 4 | 5 | 5 | 4 | **18** |
| **D3** | **Independent, well-powered cohort test** — apply the same corrected battery to GSE40279 (Hannum whole-blood 450K, n=656, ages 19–101, 473k CpGs), a cohort ~6× larger and far more uniform in age than iPOP. Directly tests "transitions at different ages than ~44/~60", and separates *"the method is broken"* from *"iPOP is underpowered"*. | 4 | 5 | 5 | 4 | **18** |
| D4 | Sensitivity of the *published crest ages* to the age distribution alone — resample iPOP to uniform/skewed age distributions and track how far the crest moves. | 4 | 4 | 3 | 5 | 16 |
| D5 | Systematic audit of the >100 papers using LOESS/DE-SWAN identified by Carbonneau et al.'s Google Scholar sweep. | 3 | 3 | 3 | 1 | 10 |
| D6 | Re-analysis of Lehallier et al. 2019 individual-level plasma proteome (n=4,263). | 4 | 5 | 4 | 1 | 14 |
| D7 | New estimator development (uncertainty-propagating functional-data test for age-localised change). | 3 | 4 | 4 | 2 | 13 |
| D8 | Longitudinal (within-subject) modelling of the full iPOP time series rather than the per-subject cross-section. | 4 | 3 | 3 | 2 | 12 |
| D9 | Simulation-only replication of Carbonneau et al.'s null studies. | 5 | 3 | 1 | 5 | 14 |

## Selected: D1, D2, D3

Together they form one argument rather than three separate studies: **D1** establishes how much
of the published signal is pipeline artifact, **D2** asks what survives correct inference on the
same data, **D3** asks whether the answer changes when the data are actually well-powered. D1 is
the necessary foundation and is cheap; D2 and D3 are where the contribution beyond Carbonneau
et al. lies. D4 and D9 are folded in as *supporting analyses inside* D1 (resampling the age
distribution and null simulation are both required to interpret D1's ablations) rather than run
as separate directions — this keeps the budget at three without discarding the mechanism checks.

## Rejected, and why

- **D5 (literature audit)** — the bottleneck is per-paper data and code access, not analysis.
  Carbonneau et al. already declined this for the same reason. Low feasibility, and the result
  would be a bibliography, not evidence.
- **D6 (Lehallier re-analysis)** — the strongest possible replication target, but the
  individual-level SOMAscan data (INTERVAL/LonGenity/Stanford cohorts) is access-controlled and
  not in the supplementary release. Only per-protein summary statistics are available
  (`datasets/supplementary/lehallier2019_MOESM3.xlsx`), which cannot support a re-run of the
  sliding-window analysis. **Reopen only if individual-level access is obtained.**
- **D7 (new estimator)** — the right long-run answer and explicitly what Carbonneau et al. call
  for, but developing *and* validating a new functional-data method is a paper of its own. D2
  takes the tractable subset: apply existing, well-understood estimators correctly.
- **D8 (longitudinal modelling)** — the raw longitudinal objects *are* in hand (882–1,440 samples
  per layer, recovered from git history), so this is feasible in principle. Deferred because
  mean follow-up is ~1.7 years against a 50-year age range: within-subject slopes cannot resolve
  a transition at 44 vs. 60. It answers a different question than the hypothesis asks.
- **D9 (simulation-only)** — near-zero information gain as a standalone direction; the code and
  results already exist in `code/artifactual-waves-of-aging/simulations/`. Retained only as the
  calibration check inside D1.

## Falsifiable predictions

The hypothesis is stated strongly, and D2/D3 are set up so it can lose:

- **Supports the hypothesis:** D1 shows crest location tracks the age-density derivative and
  moves under resampling; the valid permutation null reproduces the crests; D2 finds no
  transition surviving calibrated testing, or finds one whose CI excludes 44 and 60; D3 finds no
  concentration of change at 44/60 in a 6×-larger cohort.
- **Refutes the hypothesis:** D2 finds a nonlinearity that survives cross-validated model
  selection *and* the shuffle-before-LOESS permutation null, with a breakpoint CI covering
  ~44 or ~60; **and** D3 independently localises change near the same ages. That outcome would
  mean the pipeline was wrong but the biological conclusion was right, and must be reported as such.

A likely intermediate outcome — real monotone age association with no defensible *discrete*
transition — refutes the "discrete waves" claim while sparing "nonlinear aging", and should be
reported as that distinction rather than collapsed into either extreme.

## Constraints carried into Phase 2

1. **No R.** LOESS is reimplemented in Python (`code/loess_py/rloess.py`) and validated against
   the published R output: median relative error 1.8e-3, median per-variable correlation 0.9989.
   Residual disagreement is R's `surface="interpolate"` kd-tree approximation plus CV span
   tie-breaking; document it, do not treat the Python fit as bit-identical.
2. **The published DE-SWAN input is post-LOESS.** Any comparison must state which matrix it used.
   Both are provided (`datasets/ipop_cross_section/` raw, `datasets/ipop_loess_R_csv/` R-LOESS).
3. **iPOP is small.** 77–102 subjects per layer; edge windows can fall to n=5 in one half. Report
   per-window group sizes alongside every DE-SWAN curve — the critique's whole point is that
   these drive the curve.
4. **Multiplicity is two-dimensional.** BH within a window (as published) does not control error
   across the 25 window centres, and windows overlap heavily. Family-wise control across the
   curve is required before any crest is called.
