# Literature Review

*Reproducing and Extending the LOESS-DE-SWAN Artifact Analysis: Do Claimed Waves of Molecular
Aging Survive Rigorous Statistical Testing?*

---

## 1. Research area overview

Since 2019 a body of work has claimed that human molecular aging is not gradual but proceeds in
**discrete waves** — ages at which unusually many biomolecules change at once. Two papers
established the claim and the method:

- **Lehallier et al. 2019** (*Nature Medicine*) introduced **DE-SWAN** (Differential Expression —
  Sliding WINdow ANalysis) on 2,925 plasma SOMAmers in 4,263 people, reporting crests at
  **34, 60 and 78**.
- **Shen et al. 2024** (*Nature Aging*) applied LOESS smoothing followed by a modified DE-SWAN to
  10 omics layers from 108 iPOP participants, reporting crests at **~44 and ~60**. This paper has
  ~400 citations in under two years and is the specific target of our hypothesis.

The method spread fast. Carbonneau et al. counted **>100 papers (2019–2026)** using LOESS,
DE-SWAN or both on aging data. The downstream corpus in `papers/` shows the pattern: undulating
brain aging, undulating metabolic aging, midlife "critical windows", "tipping points" in brain
aging, a microbiome "phase transition". A 2026 *Nature Reviews Genetics* review, *Embracing
non-linearity in human ageing*, treats nonlinear aging as established.

In **June 2026** that consensus was directly challenged.

## 2. The pivotal paper

### Carbonneau, Shutta, Miller, Snyder, Shen & Quackenbush (2026), bioRxiv — *LOESS and DE-SWAN can induce artifactual "waves" of molecular aging*

Filed on bioRxiv under the type **"contradictory results"**. Its author list is the single most
important fact about it: **Xiaotao Shen and Michael Snyder, first and last authors of the Nature
Aging waves paper, are co-authors** — they "clarified the original analytical workflow,
independently reviewed the reanalysis, contributed to the interpretation of affected conclusions".
This is not an outside attack that the original authors dispute; it is a jointly-endorsed
correction of the statistical basis of their own result.

**Three artifacts, each demonstrated on data constructed to contain no signal.**

**Artifact 1 — LOESS manufactures interpretable trajectories from noise.** The common
"smooth-cluster-interpret" workflow (LOESS each molecule → hierarchically cluster the smoothed
curves → interpret centroids → enrichment analysis) was run on 1,000 molecules drawn i.i.d. from
N(0,1), independent of age, under four age distributions including the real iPOP one. Raw data
give a flat dendrogram; after LOESS, clean clusters appear whose centroids resemble published
"characteristic aging trajectories". More clusters look *more* interpretable — that is overfitting,
not resolution.

**Artifact 2 — LOESS+DE-SWAN produces waves that no valid null supports.** Points on a fitted
LOESS curve are not independent observations: each is estimated from overlapping data and
constrained to be smooth. Feeding them to a test that assumes independence is invalid, and
especially so for a rank-based test, because a smooth curve is nearly monotone within a window so
ranks separate almost perfectly. On a single null gene, DE-SWAN gives p = 0.75 while LOESS+DE-SWAN
gives p = 1.8 × 10⁻⁵. Across 1,000 null molecules, LOESS+DE-SWAN declares a **majority**
significant at FDR 0.05 at *every* window centre; DE-SWAN alone correctly finds nearly none.
Crucially, on **null** data with the iPOP age distribution, two crests appear — one in the early
40s, one in the early 60s. On the real iPOP transcriptome, LOESS+DE-SWAN reproduces the published
peaks and DE-SWAN without LOESS finds **zero** transcripts at FDR 0.05.

**The permutation discrepancy.** Shen et al. reported a permutation test in which the peaks
disappeared. Carbonneau et al. ran one in which they did not. The difference is a single ordering
step:

| | procedure | valid? |
|---|---|---|
| Algorithm 1 (Carbonneau) | shuffle ages → **then** LOESS → DE-SWAN | ✅ |
| Algorithm 2 (Shen, Suppl. Fig. 4c) | LOESS → **then** shuffle ages → DE-SWAN | ❌ |

Permutation validity requires exchangeability under the null. Subject-level data are
exchangeable; LOESS-interpolated values are not, because LOESS imposes index-dependent structure.
Shuffling after smoothing destroys that structure and produces a null that is far too easy to
beat. Under the correct ordering the peaks survive permutation — i.e. they are indistinguishable
from artifact.

**Artifact 3 — DE-SWAN alone turns linear data into nonlinear curves.** A DE-SWAN curve is a count
of significant molecules, which confounds *how much is changing* with *the power to detect it*.
Power varies along the age axis through: **(a) window sample size** under non-uniform age
distributions — different age distributions give profoundly different DE-SWAN curves from
identically-generated linear data; **(b) intra-window sample density** — with the same linear
trend and the same 50/50 split, samples massed at the window boundaries give p = 2.3 × 10⁻⁵⁶ while
samples massed at the centre give p = 1.3 × 10⁻¹⁷; **(c) heteroskedasticity**, which is
biologically plausible and produces a different curve shape for each variance pattern;
**(d) outlier clusters**. They also note that the `DEswan` package's quantile-based window centres
do *not* fix (a), because window length stays fixed.

**What they do *not* claim.** They are explicit: this does "not rule out nonlinear molecular aging
or age-associated transitions that may be detectable using **other cohorts and statistical
models**", and it is "not a direct refutation" of individual downstream papers. **Our hypothesis
as stated is stronger than their conclusion**, which is precisely where the remaining work is.

**Their recommendations** — the design brief for Phase 2: critically evaluate LOESS on sparse
data and never analyse LOESS output without propagating its correlation structure; permute at the
level of exchangeability implied by the null; characterise any method on null simulations before
applying it; prefer effect sizes with uncertainty over counts of significant results.

**Gaps they name.** Current methods cannot exploit longitudinal data; between-molecule correlation
is ignored, so it is unknown whether crests reflect biology or just assay composition (a pathway
with many measured members contributes many correlated "discoveries"); and no replacement
estimator is offered.

## 3. What the original papers actually did

### Shen et al. 2024 — pipeline as specified in Methods

1. Cross-section: mean of each participant's **healthy** visits; age = mean age across visits.
2. Confounders (BMI, sex, IRIS, ethnicity) regressed out; **residuals** used downstream.
3. LOESS per molecule, span chosen by CV over {0.3, 0.4, 0.5, 0.6}, predicted on a **half-year
   grid from 26 to 75** (99 points).
4. Fuzzy c-means (`Mfuzz`) on smoothed curves → 11 clusters, membership > 0.5.
5. Modified DE-SWAN: 20-year window, **Wilcoxon** younger vs older half, **BH within each window
   centre**, 1-year steps, plotted 40–65. Robustness by bucket width {15,20,25,30} and q-threshold
   {1e-4, 1e-3, 1e-2, 0.05}.
6. Permutation: "phenotypes of the individuals are randomly permuted" — the invalid post-LOESS null.

Two admissions in the paper's own text deserve emphasis: the P values in Figs. 1d/1e are
**unadjusted**, and "data distribution was assumed to be normal, but this was not formally tested".

Independent of the smoothing issue, the multiplicity structure is incomplete: BH within a window
does not control error **across** the 25 heavily-overlapping window centres. The quantity being
interpreted — the *location of the maximum* of a curve of counts — has no stated sampling
distribution at all.

### Lehallier et al. 2019

Original DE-SWAN uses a **linear model** with type-II ANOVA, not Wilcoxon, and no LOESS
pre-smoothing. So Artifact 2 does not apply to it directly; **Artifact 3 does**. Its crests
(34/60/78) sit in a cohort assembled from several studies, so age distribution and cohort effects
are entangled by construction. Individual-level data are access-controlled.

## 4. Methods available for the corrected analysis

Everything below is implementable in the installed Python stack.

| Purpose | Method | Notes |
|---|---|---|
| Is the age relationship nonlinear at all? | Penalised spline / GAM with CV or REML smoothness selection | `statsmodels` GAM, or B-spline basis + ridge; test vs. linear by nested comparison |
| Where does a transition occur? | Segmented regression, breakpoint + CI | Muggeo 2003 (`segmented`); Davies 1987 test — nuisance parameter unidentified under H₀, so use permutation calibration |
| Multiple changepoints | PELT / binary segmentation with a penalty | Killick, Fearnhead & Eckley 2012; `ruptures` if added |
| Multiplicity within a window | Benjamini–Hochberg | Already used; correct as far as it goes |
| Multiplicity **across** windows and under correlation | Benjamini–Yekutieli, or permutation-based FWER over the whole curve | Windows overlap and molecules are correlated — BH within window is not enough |
| Valid null for the whole pipeline | Permute **ages before any smoothing**, rerun end to end | Carbonneau Algorithm 1 |
| Equalising power along the age axis | Fixed *n* per half-window rather than fixed width | Directly targets Artifact 3(a); note it does not fix 3(b) or 3(c) |
| Heteroskedasticity | Location–scale modelling, or variance-stabilising transform | Artifact 3(c) |

**Conceptual reference.** Kriegeskorte et al. 2009 (`papers/kriegeskorte2009_circular_analysis.pdf`)
is the canonical statement of the failure mode: selecting/transforming data with a procedure and
then testing on the same data with a test that ignores the transformation. LOESS+DE-SWAN is a
textbook instance.

**Classic references** (cited, not downloaded — standard and paywalled): Cleveland 1979;
Cleveland & Devlin 1988; Muggeo 2003, 2008; Wood 2011, 2017; Killick et al. 2012; Davies 1987;
Benjamini & Hochberg 1995; Benjamini & Yekutieli 2001.

## 5. Datasets in the literature

| Dataset | Used by | Status here |
|---|---|---|
| iPOP multi-omics (108 subjects, 10 layers) | Shen et al. 2024; Carbonneau et al. 2026 | ✅ recovered, all 10 layers — `datasets/ipop_cross_section/` |
| Lehallier plasma proteome (4,263 subjects) | Lehallier et al. 2019 | ❌ access-controlled; summary stats only |
| GSE40279 (656 subjects, 473k CpGs, ages 19–101) | Hannum et al. 2013 | ✅ downloaded — independent validation cohort |
| UK Biobank Olink plasma proteome | Several downstream papers | ❌ application required |

## 6. Gaps and opportunities

1. **One cohort, one modality.** Carbonneau et al. analysed iPOP transcriptomics (proteomics in
   supplementary). Nothing is known about whether the artifact reproduces across all 10 layers, or
   what happens in a cohort large enough to have real power.
2. **Diagnosis without replacement.** They recommend but do not deliver a corrected estimator.
   Nobody has asked: *given* the pipeline is broken, does any age-localised change survive proper
   inference, and where?
3. **No quantitative decomposition.** No one has partitioned the published crest into
   artifact vs. residual signal, or reported the *achieved* FDR of the published procedure.
4. **Power confound never controlled.** No published DE-SWAN analysis reports per-window group
   sizes alongside the curve, though these are the mechanism behind Artifact 3.
5. **Between-molecule correlation ignored.** Crest heights count correlated molecules as
   independent discoveries; a well-represented pathway inflates a crest by assay composition alone.
6. **Longitudinal data unused.** iPOP has ~47 samples per participant, discarded by averaging.

## 7. Recommendations for our experiments

**Datasets.** Primary: `datasets/ipop_cross_section/` (all 10 layers; transcriptome and proteome
first, since those are the published crest drivers). Independent: GSE40279. Do not use the DEswan
package's bundled `agingplasmaproteome` as data — it is a de-identified fixture.

**Baselines to reproduce before anything else.**
1. Published LOESS+DE-SWAN on iPOP — must show crests near 44 and 60.
2. DE-SWAN without LOESS on the same raw data — must show near-zero discoveries.
3. Shuffle-after-LOESS null (Shen) vs. shuffle-before-LOESS null (Carbonneau) — must show peaks
   disappearing under the former and persisting under the latter.

Baselines 1 and 2 are **already verified** here on proteomics
(`code/deswan_py/smoke_test_proteomics.py`): 174–245 of 302 proteins significant in every window
with a peak at **44** under LOESS+DE-SWAN, versus **1** significant test in total without LOESS.

**Metrics.** Report, for every configuration: (a) count of significant molecules per window with
per-window `(n_young, n_old)` alongside; (b) achieved FDR / type-I error on matched null data;
(c) crest location **with a confidence interval**, not a point estimate; (d) effect sizes with
uncertainty, per Carbonneau et al.'s recommendation; (e) sensitivity of crest location to bucket
width, q-threshold, span, and to resampling the age distribution.

**Methodological cautions.**
- Never test on LOESS output without propagating its correlation structure; if it is done for
  comparability with the published result, label it as the *reproduction* arm, not an analysis.
- Permute before smoothing. Always.
- Report per-window sample sizes with every DE-SWAN curve.
- Control multiplicity in **both** dimensions — across molecules and across window centres.
- Distinguish three distinct claims, which the literature routinely conflates: *aging is
  nonlinear* (weak, probably true), *there exist discrete transition ages* (the actual claim under
  test), and *those ages are 44 and 60* (the specific claim). Our results should address each
  separately.
- Preserve the possibility of refutation: if a transition survives calibrated testing with a CI
  covering 44 or 60, and replicates in GSE40279, report that the biological conclusion stood up
  despite the flawed method.
