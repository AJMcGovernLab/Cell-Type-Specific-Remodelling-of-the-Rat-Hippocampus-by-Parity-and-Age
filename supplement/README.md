# Supplementary Data Bundle

Source data CSVs for every figure and supplementary figure in:

> **Cell Type-Specific Remodelling of the Rat Hippocampus by Parity and Age**
> McGovern, Duarte-Guterman, Galea

All files here are copies of the corresponding reproduced outputs in `../outputs/`.

## Main figures

### `Figure1_SourceData/`
Reference-sex deconvolution performance metrics.
- `normalized_metrics_table.csv` — 9 rows (3 datasets × 3 sex configurations). Columns: `reference`, `dataset`, `dataset_label`, `sex`, `technology`, `n_ref_cells`, `n_detected_celltypes`, `types_per_k_cells`, `entropy_per_k_cells`, `sparsity_per_k_cells`, `mean_prop_per_k_cells`, `detected_per_sample_per_k`, `diversity_index`. The entire Figure 1 is derived from this file.

### `Figure2_SourceData/`
Harmonization of 349 reference cell types into 23 biologically coherent categories.
- `refined_cell_type_mapping.csv` — 349 rows × 8 cols. The authoritative mapping from original cell-type names to harmonized `unified_name` labels (CA3_dorsal, Sst_IN, etc.). Underlies Figure 2a-c.
- `correlation_summary_stats.csv` — within/between-cluster Pearson correlation statistics cited in §3.2 (mean_within, sd_within, mean_between, sd_between, Wilcoxon W and p).

### `Figure3_SourceData/`
Three-way ANOVA on harmonized cell-type proportions.
- `harmonized_anova_results.csv` — 14 cell types × 7 effect-term p-values (age, parity, region, all 2-way interactions, 3-way). Source of every p-value cited in §3.3.1 (Sst_IN age 3.66×10⁻⁸, DG/CA1/CA2/CA3_dorsal region 6.51×10⁻²⁹, etc.).

### `Figure4_SourceData/`
Weighted GSEA meta-analysis for age/region effects on harmonized cell types.
- `cell_type_vulnerability_ranking.csv` — 17 rows. Source of Figure 4a-b. Columns: `effect`, `cell_type`, `mean_meta_score`, `high_conf_total`, `vulnerability_score`, `top_pathway`, `top_meta_score`. (Sst_IN age vulnerability = 8.192; top pathway "catalytic activity" meta-score 7.114.)
- `age_high_confidence_pathways.csv` — 37 rows. Per-cell-type high-confidence (Moderate+) pathways for the age effect.
- `region_high_confidence_pathways.csv` — 3 rows. Region effect high-confidence pathways.
- `interaction_high_confidence_pathways.csv` — 11 rows. Age × Region interaction effect pathways.

### `Figure5_SourceData/`
Parity effects on unharmonized cell-type proportions + cell-specific enrichment + pathway overlaps.
- `Figure5a_all_cell_types_parity_ANOVA.csv` — 27 rows. All individual-dataset cell types tested with three-way ANOVA (`fixed_individual_anova_results.csv`).
- `Figure5b_parity_responsive_4_cells.csv` — the 4 parity-significant cell types (356_CA3-do p=0.010, 358_CA3-do p=0.0086, 376_Astro p=0.022, 78_Sst_HPF p=0.029).
- `Figure5c_356_CA3_do_meta_analysis.csv` — CA3 dorsal (10x 2020) enriched pathway meta-analysis (50 rows).
- `Figure5c_358_CA3_do_meta_analysis.csv` — CA3 dorsal (Smart-seq 2019) meta-analysis (1,557 rows).
- `Figure5d_376_Astro_meta_analysis.csv` — Astrocyte (10x 2020) meta-analysis (4 rows).
- `Figure5e_78_Sst_HPF_meta_analysis.csv` — Sst_HPF (10x 2020) meta-analysis (77 rows).
- `Figure5cde_high_confidence_pathways_combined.csv` — Moderate+ tier pathways combined across all 4 cell types (source of §3.4.2 top values: GSK3B/NFE2L2 meta=7.88 NES=−3.55, Cyclin D meta=7.86 NES=−3.42).
- `Figure5fi_pathway_overlap_summary.csv` — 7-region Venn counts + percentages (419/110/22/11/25/19/13, source of Figure 5f-i).

## Supplementary figures

### `SupplementaryFigure1_SourceData/` — gene-level three-way ANOVA
- `threeway_anova_results_CORRECTED.csv` — 12,516 genes × 32 columns (p, FDR, F-stat, significance flags, R²).
- `significant_genes_{age,parity,region,age_region,age_parity}_CORRECTED.csv` — per-effect gene lists (age 8,163; region 6,086; parity 8; age×region 1,798; age×parity 7).
- `anova_significance_summary_CORRECTED.csv` — effect totals + percentages.

### `SupplementaryFigure2_SourceData/` — random-forest ∩ PCA gene overlap (Venn)
- `rf_vs_pc_parity_set_sizes.csv` — set sizes and overlaps between the random-forest and PCA (PC5/6/8) parity-gene sets.

### `SupplementaryFigure3_SourceData/` — DESeq2 binary-parity differential expression (design: ~ age + parity_binary)
- `significant_genes_parity.csv` — FDR < 0.05 genes (volcano highlights).
- `all_genes_parity_deseq2.csv` — full per-gene DESeq2 result (all genes; full volcano).

### `SupplementaryFigure4_SourceData/` — PCA driver genes + PC × model contribution
- `PC{5,6,8}_Parity_Binary_improved_drivers.csv` — driver gene lists (63 / 128 / 168 genes).
- `pc_model_driver_contribution.csv` — mean driver-gene contribution per PC × statistical model (heatmap data).

### `SupplementaryFigure5_SourceData/` — random-forest classification
- `top23_genes_with_names.csv` — ranked top-23 genes from the multi-class RF.
- `multiclass_feature_importance.csv` — full importance ranking across all genes.
- `cv_accuracy_by_geneset_20_40.csv` — cross-validation accuracy by gene-set size (peak at 23).

## Supplementary Data (non-figure pathway tables)

### `SupplementaryData_Pathways/`
Full pathway lists for each region of the 3-cell-type Venn (§3.4.3, Fig 5f-i). Cited counts in paper but full lists not shown.
- `All_three_cell_types_pathways.csv` — 419 pathways enriched in CA3 ∩ Astro ∩ Sst.
- `Astrocytes_and_Sst_HPF_pathways.csv` — 110 Astro ∩ Sst pairwise.
- `CA3_Dorsal_and_Astrocytes_pathways.csv` — 22 CA3 ∩ Astro.
- `CA3_Dorsal_and_Sst_HPF_pathways.csv` — 11 CA3 ∩ Sst.
- Three "only" lists — pathways unique to each single cell type.

## Reproduction chain

Each file here is a copy of the corresponding `../outputs/` file produced by the pipeline in
[../scripts/](../scripts/). Run `Rscript ../scripts/run_all_robust.R` to regenerate `../outputs/`;
see the top-level [README](../README.md) for the inputs and the figure → script map.
