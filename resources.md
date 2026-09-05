# Resources Catalog

## Summary

| Category | Count | Location |
|---|---:|---|
| Papers (11 PDF + 6 full-text XML) | 17 | `papers/` |
| Datasets | 4 | `datasets/` |
| Cloned repositories | 3 | `code/` |
| Modules written here | 2 | `code/loess_py/`, `code/deswan_py/` |

Details: `papers/README.md`, `datasets/README.md`, `code/README.md`.
Synthesis: `literature_review.md`. Direction budget: `planning.md`.

---

## Papers

| Title | Authors | Year | File | Key info |
|---|---|---|---|---|
| LOESS and DE-SWAN can induce artifactual "waves" of molecular aging | Carbonneau, Shutta, Miller, Snyder, Shen, Quackenbush | 2026 | `carbonneau2026_loess_deswan_artifactual_waves.pdf` | **The paper this project reproduces.** Original waves authors are co-authors. Code: QuackenbushLab/artifactual-waves-of-aging |
| Nonlinear dynamics of multi-omics profiles during human aging | Shen, Wang, Zhou, Zhou, Hornburg, Wu, Snyder | 2024 | `shen2024_nonlinear_multiomics_aging.pdf` | *Nature Aging*, ~400 cites. Source of the ~44/~60 claim. Full Methods extracted to `search/shen2024_fulltext.txt` |
| Undulating changes in human plasma proteome profiles across the lifespan | Lehallier et al. | 2019 | `lehallier2019_undulating_plasma_proteome.pdf` | *Nat Med*, ~830 cites. Origin of DE-SWAN; crests 34/60/78 |
| Plasma proteomics identify biomarkers and undulating changes of brain aging | Yin et al. | 2024 | `yin2024_plasma_proteomics_undulating_brain_aging.pdf` | Direct methodological descendant |
| Nonlinear DNA methylation trajectories in aging male mice | Sziráki et al. | 2024 | `sziraki2024_nonlinear_methylation_aging_mice.pdf` | Nonlinear claim, different organism/modality |
| Sex-specific nonlinear DNA methylation aging trajectories | Grolaux et al. | 2026 | `grolaux2026_sex_specific_nonlinear_methylation.fulltext.xml` | Closest analogue to our GSE40279 test |
| Brain aging shows nonlinear transitions, midlife "critical window" | Mujica-Parodi et al. | 2025 | `mujica2025_brain_aging_nonlinear_transitions_pnas.fulltext.xml` | *PNAS*; independent modality, midlife transition |
| Identifying tipping points during healthy brain aging | Wu et al. | 2025 | `wu2025_tipping_points_brain_aging.fulltext.xml` | snRNA-seq "tipping point" variant |
| Plasma proteomics reveals undulating changes in metabolic aging | Zhang et al. | 2025 | `zhang2025_plasma_proteomics_undulating_metabolic_aging.fulltext.xml` | Downstream application |
| Data-driven identification and classification of nonlinear aging patterns | Lai et al. | 2023 | `lai2023_nonlinear_aging_patterns_human_genomics.fulltext.xml` | Smooth-cluster-interpret workflow |
| Global metagenomic atlas of aging: microbiota phase transition | Fu et al. | 2026 | `fu2026_microbiota_phase_transition.pdf` | "Phase transition" claim |
| Organ aging signatures in plasma proteome track health and disease | Oh et al. | 2023 | `oh2023_organ_aging_plasma_proteome.pdf` | *Nature*; plasma proteome aging context |
| Embracing non-linearity in human ageing | — | 2026 | `natrevgenet2026_embracing_nonlinearity.pdf` | *Nat Rev Genet* review of the debate |
| Everything Everywhere All at Once: Unraveling the Waves of Aging | — | 2025 | `immunity2025_waves_of_aging.fulltext.xml` | *Immunity* commentary |
| Circular analysis in systems neuroscience | Kriegeskorte et al. | 2009 | `kriegeskorte2009_circular_analysis.pdf` | Canonical statement of the failure mode |
| Personal omics profiling reveals dynamic molecular…(iPOP) | Chen et al. | 2012 | `chen2012_ipop_personal_omics.fulltext.xml` | Origin of the iPOP cohort |
| Genome-wide methylation profiles reveal quantitative views of human aging rates | Hannum et al. | 2013 | `hannum2013_methylation_aging_rates.fulltext.xml` | Source publication for GSE40279 |

## Datasets

| Name | Source | Size | Task | Location | Notes |
|---|---|---|---|---|---|
| iPOP cross-section (10 omics) | `jaspershen-lab/ipop_aging` **git history** | 51–52,460 vars × 77–102 subjects | Reproduce/refute 44 & 60 crests | `datasets/ipop_cross_section/` | Transcriptome 8,556 × 97 — matches paper exactly |
| iPOP raw longitudinal | same | 630–1,440 samples/layer | Optional longitudinal extension | `datasets/ipop_raw/`, `datasets/ipop_csv/` | Retains visit structure, health status, covariates |
| iPOP published R LOESS output | same repo, working tree | 8 layers × 99 grid points | Ground truth for Python LOESS | `datasets/ipop_loess_R_csv/` | Validation target |
| Supplementary tables | Springer static content | 33 MB | Crest lists, cluster memberships, per-protein stats | `datasets/supplementary/` | Shen ×3, Lehallier ×3 |
| GSE40279 (Hannum 450K) | NCBI GEO | 473,034 × 656, ages 19–101 | Independent validation | `datasets/GSE40279/` | 1.24 GB source → float32 npy + 20k-CpG parquet |

**The iPOP recovery is the non-obvious part.** Shen et al.'s data-availability statement points to
portals that do not serve the analysed matrices, and the analysis repo's `ignore_large_files.sh`
excluded every file >5 MB — so the raw pre-LOESS objects are absent from `main`. They survive in
git history and are recovered by `datasets/recover_ipop_from_git.py` (verified byte-identical on
re-run). Without them, DE-SWAN-without-LOESS — the central comparison — is not runnable.

## Code repositories

| Name | URL | Purpose | Location | Notes |
|---|---|---|---|---|
| artifactual-waves-of-aging | github.com/QuackenbushLab/artifactual-waves-of-aging | The critique's own code | `code/artifactual-waves-of-aging/` | Python sims runnable as-is; ships **precomputed permutation results** (80 transcriptomic). R notebooks need R |
| ipop_aging | github.com/jaspershen-lab/ipop_aging | Original Shen et al. analysis | `code/ipop_aging/` | 2.9 GB. Contains crest molecule lists and the LOESS-interpolated matrix DE-SWAN was run on |
| DEswan | github.com/lehallib/DEswan | Reference DE-SWAN R package | `code/DEswan/` | Original linear/type-II-ANOVA test; quantile window-centre default |
| loess_py *(written here)* | — | R-compatible LOESS in Python | `code/loess_py/` | Validated: median rel. err 1.8e-3, median corr 0.9989 vs published R output |
| deswan_py *(written here)* | — | DE-SWAN in Python, both test variants | `code/deswan_py/` | Returns per-window group sizes; includes a passing end-to-end smoke test |

## Notes

### Search strategy

The paper-finder service was unavailable (`localhost:8000` not running), so literature search was
done manually against Europe PMC, Semantic Scholar (citation graph of both anchor papers — 396 +
836 citing papers screened by title regex), bioRxiv, and arXiv `stat.ME`. Screening the Shen 2024
citation graph is what surfaced the Carbonneau et al. critique, which does not appear under
obvious keyword queries.

### Selection criteria

Tier 1 = the three papers the hypothesis is directly about, deep-read in full. Tier 2 = downstream
papers inheriting the method, which define the scope of what is at stake. Tier 3 = methodological
and provenance context. Classic statistics references (Cleveland, Muggeo, Wood, Killick,
Benjamini–Hochberg/Yekutieli, Davies) are cited but not downloaded: paywalled, standard, and
already implemented in the installed stack.

### Challenges encountered and how they were resolved

| Problem | Resolution |
|---|---|
| paper-finder service down | Manual multi-source search incl. citation-graph screening |
| Shen et al. data not at the URLs given in the paper | Recovered from `ipop_aging` git history (`recover_ipop_from_git.py`) |
| `.rda` files are S4 `mass_dataset` objects; **no R in this environment** | Parsed with Python `rdata`; one over-strict assertion patched in `convert_ipop.py` |
| Pipeline requires R's `loess`; `statsmodels.lowess` is degree-1 and not equivalent | Reimplemented R's degree-2 tricube LOESS in Python and validated against the authors' own published R output |
| PMC blocked direct PDF scraping | Publisher PDFs via nature.com/PNAS; Europe PMC and NCBI eutils JATS XML otherwise |
| Semantic Scholar rate-limited (HTTP 429) | Switched to Europe PMC + arXiv for the methodological sweep |
| Lehallier individual-level data access-controlled | Documented as a hard limitation; direction D6 pruned in `planning.md` |

### Gaps and workarounds

- **Lehallier n=4,263 proteome is unavailable.** The 34/60/78 crests cannot be re-analysed
  directly. Workaround: GSE40279 (n=656) as an independent, well-powered cohort, plus the
  per-protein summary statistics in `lehallier2019_MOESM3.xlsx`.
- **No R.** Removed as a blocker by `loess_py` and `deswan_py`. `Mfuzz` fuzzy c-means has no
  direct Python equivalent; if cluster reproduction is needed, `skfuzzy.cmeans` is close but must
  be validated against `ipop_aging`'s saved cluster assignments (which are in the repo).
- **The published DE-SWAN input is post-LOESS.** Both matrices are provided so every comparison
  can state which it used.
- **`run_loess` is slow** (~1 min per 300 variables). Vectorise before permutation studies.

## Recommendations for experiment design

1. **Primary datasets** — `datasets/ipop_cross_section/` (transcriptome and proteomics first;
   these drive the published crests), then the other 8 layers for breadth. Independent validation:
   GSE40279.
2. **Baselines** — (i) published LOESS+DE-SWAN, (ii) DE-SWAN without LOESS, (iii) shuffle-after-LOESS
   vs shuffle-before-LOESS permutation nulls. (i) and (ii) are already verified here on proteomics.
3. **Metrics** — significant-molecule counts *with per-window group sizes*; achieved FDR on matched
   null data; crest location **with a CI**; effect sizes with uncertainty; sensitivity to bucket
   width, q-threshold, span, and age-distribution resampling.
4. **Code to reuse** — `code/deswan_py/deswan.py` and `code/loess_py/rloess.py` are the workhorses;
   `code/artifactual-waves-of-aging/simulations/DE-SWAN/` supplies ready null-simulation
   configurations for all of Artifact 3; its `permutation_rslts/*.csv` saves days of compute.
5. **Framing** — the hypothesis as written is stronger than what Carbonneau et al. concluded. Keep
   the three claims separate (*nonlinear*, *discrete transitions exist*, *they are at 44/60*), and
   keep the refutation path open: a transition that survives calibrated testing and replicates in
   GSE40279 would mean the method was wrong but the conclusion right, and must be reported that way.
