# Supplementary Methods

Supplementary detail for *Cell Type-Specific Remodelling of the Rat Hippocampus
by Parity and Age* — parameters, formulas, algorithm specifics, and
reproducibility notes that **complement rather than repeat** the main-text
Methods. Section numbering mirrors the manuscript (2.X ↔ S2.X); only sections
with additional detail beyond the main text are included.

**GEO accession:** [GSE329776](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329776)
**Code repository:** https://github.com/AJMcGovernLab/ReproductiveExperienceAndAgeDeconvolution

---

## Software and package versions

| Component | Version | Purpose |
|---|---|---|
| BCL2FASTQ | 2.20.0.422 | Demultiplexing of Illumina BCL output |
| Cutadapt | (latest at time of analysis) | Adapter and quality trimming |
| FastQC | 0.11.8 | Pre- and post-trim sequence QC |
| STAR | 2.6.1a_08-27 | Read alignment to *Rnor_6.0* (Ensembl 84) |
| HTSeq-count | 0.11.0 | Per-gene read quantification |
| SCDC | 0.0.0.9000 | Cell-type deconvolution |
| DESeq2 | 1.40.0 | Variance-stabilizing transform; binary parity DE |
| clusterProfiler | 4.8.0 | GSEA against GO / KEGG / Reactome / MSigDB |
| randomForest | 4.7-1.1 | Parity-classification feature selection |
| tidyverse | 2.0.0 | Data wrangling |
| ggplot2 | 3.4.2 | Visualization |

Original analysis under **R 4.3.0**; the public reproduction pipeline is verified
under **R 4.5.3**. Random seeds, environment manifests, and all driver scripts
are checked in to the project repository.

---

## §S2.3 Reference dataset cell counts and performance metrics

Per-reference cell counts in each sex configuration:

| Reference | Platform | Female cells | Male cells | Mixed cells |
|---|---|--:|--:|--:|
| Allen 10X Genomics whole-brain (2020) | 10X Chromium | 22,307 | 59,108 | 81,964 |
| Yao hippocampus 10X | 10X Chromium | 26,314 | 62,272 | 88,674 |
| Mouse-brain Smart-seq2 (2019, NeMO `dat-iye7gkp`) | Smart-seq2 | 3,693 | 1,849 | 5,869 |

Configurations were compared (Figure 1) on four size-normalized metrics:

- **Shannon entropy per 1,000 reference cells**: −Σ(pᵢ × log₂(pᵢ)) / 1,000, where pᵢ is the proportion of cell type *i*.
- **Relative cell-type detection rate**: proportion of maximum detectable cell types identified in each configuration.
- **Modified Simpson's diversity index**: standard Simpson's index normalized by log₁₀(n_ref_cells + 1) to control for reference size.
- **Detection efficiency**: cell types detected per 1,000 reference cells.

Female-matched configurations consistently outperformed mixed-sex and male-only
configurations on the 10X-based references; all deconvolution used the
**female-only** configuration.

## §S2.4 Cell-type harmonization detail

The harmonization computed pairwise Pearson correlations between all 349
cell-type signatures (extracted from the SCDC basis matrices), then performed
hierarchical clustering with average linkage on the distance matrix `1 − r` at a
**height cutoff of 0.2**. Per-reference reduction:

| Reference | Originals | Harmonized | Reduction |
|---|--:|--:|--:|
| mouse10x_2020 | 169 | 16 | ~10× |
| mouse_smartseq_2019 | 125 | 10 | ~12× |
| yao_hippo_10x | 55 | 1 | ~55× |
| **Total (sum across datasets)** | **349** | **27** | ~13× |

The **27** entries are the per-dataset-summed harmonized categories (Figure 2A);
collapsing categories that recur across references gives **23** unique harmonized
categories (the cell-class breakdown 12 neuronal + 6 glial + 3 vascular + 2
immune = 23 uses this deduplicated convention).

Two mapping files are used at different stages:
- [`cell_type_mapping_table.csv`](../outputs/02_fig2/cell_type_mapping_table.csv) — first-pass harmonization (Figure 2A / §2.4 counts).
- [`refined_cell_type_mapping.csv`](../outputs/02_fig2/refined_cell_type_mapping.csv) — second-pass refinement (additional h = 0.1 sub-clustering + biological renaming, e.g. `Sst_IN`, `CA3_dorsal`, `Oligodendrocyte`), used as the row identity for the harmonized three-way ANOVA.

Within- vs between-cluster correlation statistics (reported in §3.2) are on disk
in [`correlation_summary_stats.csv`](../outputs/02_fig2/correlation_summary_stats.csv).

## §S2.5 SCDC deconvolution parameters

- **Subject parameter**: animal ID, to account for biological replicates and improve estimation stability.
- **Ensemble weighting**: equal weights across the three reference datasets; sensitivity analyses with alternative weighting schemes produced concordant proportions (Pearson **r > 0.95**).
- **Convergence and reconstruction error** were verified on every sample.
- **Technical-replicate correlation** consistently exceeded **r = 0.85**.

Proportions estimated from harmonized vs. non-harmonized annotations were
compared to assess the impact of standardization on biological inference.

## §S2.6 Directional inference for cell-type proportions

For directional inference of significant main effects, **marginal log2 fold
change** estimates were computed as the unweighted mean of within-stratum log2
fold changes across the 4 (age × region) or (age × parity) cells of the design
— equivalent to estimated marginal means from the same factorial model. This
conditioned estimator matches the directional component of the three-way ANOVA
p-value, in contrast to a simple unconditioned pooled contrast (e.g. parous vs
nulliparous across all samples), which is biased by variance attributable to the
other two factors when the design is unbalanced. Per-stratum sign counts (out of
4 strata) accompany each directional log2 fold change in the Results as a
non-parametric robustness check (e.g. "4/4 parous-higher" indicates the parous
mean exceeds the nulliparous mean in every age × region subgroup).

Cross-dataset consensus: a finding was **high-confidence** when it showed the
same direction and p < 0.05 across **≥ 2 reference datasets**; cross-dataset
combination used **Fisher's combined probability test** with adjustment for
between-dataset correlation.

## §S2.7 Differential gene expression — additional detail

### §S2.7.1 Type III sums of squares
The gene-level three-way ANOVA (`expression ~ Age × Region × Parity`, `aov()`)
used **Type III sums of squares**, justified by the unbalanced factorial design;
Benjamini-Hochberg FDR < 0.05 was applied within each of the seven model terms.

### §S2.7.2 Enrichment weighting input (DESeq2)
For functional enrichment, weighted gene scores combined effect magnitude and
statistical confidence:

```
weight = 0.5 × |log2(fold-change)|  +  0.5 × (−log10(adjusted p-value))
```

### §S2.7.3 PCA driver-gene identification
PCA was performed on the variance-stabilized (DESeq2 VST) expression matrix; the
top 10 PCs explained > 60% of total variance. Each PC was tested for factor
associations with seven nested linear models (main effects and all interaction
combinations). PC5, PC6, and PC8 were significantly associated with binary
parity (Parity_Binary main-effect model: p = 0.048, 0.026, 0.013; R² = 0.066,
0.083, 0.101).

For parity-associated PCs, driver genes were identified by ranking genes by
absolute loading and quantifying each gene's contribution as the reduction in
model R² when its influence on the PC was removed: the score
`reduced = pc_scores − (gene_loading × gene_expression)` was substituted for the
original PC scores in the parity-binary model, and the R² loss (`r2_loss`) and
parity-coefficient p-value inflation (`log_p_degradation`) were combined into a
contribution score `log_p_degradation + 100 × r2_loss`; genes above the elbow of
that distribution were retained as drivers:

| PC | Driver genes |
|---|--:|
| PC5 | 63 |
| PC6 | 128 |
| PC8 | 168 |

### §S2.7.4 Random-forest parameters and reproducibility
The **multi-class** parity label (nulliparous, primiparous, biparous; 3 classes)
was the response variable — not the binary grouping used elsewhere — to allow the
classifier to leverage primiparous–biparous differences.

- **1,000 decision trees** (`ntree = 1000`); bootstrap sampling per tree; `mtry = floor(sqrt(n_features)) = 111`; `set.seed(12345)`.
- **Mean Decrease in Accuracy (MDA)** = drop in classification accuracy when each variable is permuted; **Mean Decrease in Gini (MDG)** = total node-impurity decrease averaged across trees.
- The number of parity-predictive genes (23) was determined from the **cross-validation classification accuracy across candidate gene-set sizes, which peaked at 23 genes** (5-fold cross-validation, 3 repeats; Supplementary Figure 5). The **top 23 genes by MDA** were retained; saved-model MDA range **1.83–3.04** (top ENSRNOG00000053712 = 3.036; bottom Vrk3 = 1.832).

**Reproducibility caveat.** The trained model objects are the canonical artefacts
for the published top-23 list — loading them reproduces it exactly on any
platform. End-to-end refit from the input matrix requires the original
R/randomForest versions (R 4.3.0 + randomForest 4.7-1.1): refit preserves the
binary-RF out-of-bag error to high precision (0.283), but multiclass top-23
membership shifts by ~17/23 across patch versions due to randomForest's
classification-mode RNG-consumption pattern. The pipeline
([`09_suppfig2to5_parity_genes/02_random_forest.R`](../scripts/09_suppfig2to5_parity_genes/02_random_forest.R))
defaults to loading the saved models; `REPRO_RF_REFIT=TRUE` retrains from scratch.

## §S2.8 Functional enrichment — weighting, GSEA, and meta-analysis

### §S2.8.1 Weighting schemes and databases
Four complementary weighting schemes integrated evidence from the DESeq2, PCA,
and random-forest analyses:

| Scheme | Composition |
|---|---|
| Association-weighted | 100% gene–cell-type association score |
| Expression-weighted | 60% association  +  40% \|log₂(fold-change)\| |
| Statistical-weighted | 60% association  +  40% −log₁₀(FDR) |
| Composite-weighted | 50% association  +  25% \|log₂(fold-change)\|  +  25% −log₁₀(FDR) |

**Seven pathway database collections** were queried: GO:BP, GO:MF, GO:CC, KEGG,
Reactome, and two MSigDB collections (Hallmark and Curated). For each cell type,
weighting scheme, and percentage-rank cutoff (top 5% / 10% / 25% by signed
weighted score, contributing meta-analysis weights 0.50 / 0.40 / 0.10), the
ranked gene list was submitted to GSEA against each database.

### §S2.8.2 GSEA implementation
GSEA used **clusterProfiler v4.8.0**: `gseGO()` for the three GO ontologies,
`gseKEGG()` for KEGG, `gsePathway()` (ReactomePA) for Reactome, and `GSEA()` for
the MSigDB collections.

- Minimum gene-set size **5**; maximum **500**.
- Default gene-set permutation (clusterProfiler fgsea backend); Benjamini-Hochberg adjustment.
- Pathway magnitude = **normalised enrichment score (NES)**, comparable across gene-set sizes; sign encodes the leading-edge position in the ranked input.
- **Directional interpretation**: because inputs are ranked by signed weighted score (positive = parity-upregulated), positive NES = upregulated, negative NES = downregulated in the relevant cell type.

### §S2.8.3 Meta-score and confidence tiers
For each pathway:

```
meta-score = consensus × evidence × max(0, NES_consistency)
```

- **consensus** = proportion of methods returning the pathway with FDR < 0.05.
- **evidence** = weighted mean −log₁₀(method-level FDR), method weights {association 0.30, expression 0.25, statistical 0.25, composite 0.20}.
- **NES consistency** = `1 − (SD / (|mean| + ε))` of per-method NES, clamped at ≥ 0.
- Method p-values combined via **Fisher's combined probability test** (`pchisq(−2 × Σ log(p), df = 2·n_methods)`); the resulting meta-FDR was used for filtering.

| Tier | Method count | Meta-score |
|---|---|---|
| Ultra-high | All 4 | > 7 |
| High | ≥ 3 | > 5 |
| Moderate | ≥ 2 | > 3 |
| Method-specific | 1 | — |

A pathway was classified **cell-type-specific** if enriched in ≤ 2 cell types
with a meta-score difference > 2 between enriched and non-enriched cell types.

---

## Cross-reference index

| Manuscript § | Driver script directory |
|---|---|
| 2.3 | [`scripts/01_fig1_reference_sex/`](../scripts/01_fig1_reference_sex/) |
| 2.4 | [`scripts/02_fig2_harmonization/`](../scripts/02_fig2_harmonization/) |
| 2.5 | [`scripts/01_fig1_reference_sex/01_run_deconvolution.R`](../scripts/01_fig1_reference_sex/01_run_deconvolution.R) (SCDC; outputs staged in `checkpoints/scdc_deconvolution/`) |
| 2.6 | [`scripts/03_fig3_harmonized_anova/`](../scripts/03_fig3_harmonized_anova/) (harmonized) and [`scripts/06_fig5ab_parity_proportions/`](../scripts/06_fig5ab_parity_proportions/) (per-dataset) |
| 2.7.1 | [`scripts/04_suppfig1_gene_anova/`](../scripts/04_suppfig1_gene_anova/) |
| 2.7.2 | [`scripts/09_suppfig2to5_parity_genes/01_deseq2_binary_parity.R`](../scripts/09_suppfig2to5_parity_genes/01_deseq2_binary_parity.R) |
| 2.7.3 | [`scripts/07_pca_preprocessing/`](../scripts/07_pca_preprocessing/) and [`scripts/09_suppfig2to5_parity_genes/03_pc_driver_analysis.R`](../scripts/09_suppfig2to5_parity_genes/03_pc_driver_analysis.R) |
| 2.7.4 | [`scripts/08_rf_preprocessing/`](../scripts/08_rf_preprocessing/) and [`scripts/09_suppfig2to5_parity_genes/02_random_forest.R`](../scripts/09_suppfig2to5_parity_genes/02_random_forest.R) |
| 2.8 | [`scripts/05_fig4_age_region_enrichment/`](../scripts/05_fig4_age_region_enrichment/) (age/region) and [`scripts/10_fig5_parity_enrichment/`](../scripts/10_fig5_parity_enrichment/) (parity); overlaps in [`scripts/11_fig5_pathway_overlaps/`](../scripts/11_fig5_pathway_overlaps/) |

**End-to-end reproduction:** `Rscript scripts/run_all_robust.R` from the
[Repository root](../../Repository/) regenerates every CSV, PDF, and PNG in
[`outputs/`](../outputs/) from the staged data and checkpoints. SCDC
deconvolution (§2.5) is checkpointed because SCDC is not on CRAN
(`remotes::install_github("meichendong/SCDC")`); its output `all_results_scdc.rds` is produced by
Section 01. See the [README](../README.md) for inputs and the run order.
