# Parity Paper — Analysis Code

Analysis code for:

> **Cell Type-Specific Remodelling of the Rat Hippocampus by Parity and Age**
> McGovern AJ, Duarte-Guterman P, Galea LAM
> *Journal of Neuroendocrinology* (under review)

[![GEO](https://img.shields.io/badge/GEO-GSE329776-blue)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329776)

This repository contains **the scripts that build every reported figure and value**, plus the
produced results (`outputs/`) and figure source-data (`supplement/`). It does **not** bundle the
large inputs (bulk counts, single-cell references) or the computed intermediates — instead it tells
you exactly **where to download the inputs** and **which script produces each result**, so the full
analysis can be re-run. (A complete local archive with all data + intermediates is maintained
separately and deposited at Zenodo / GEO.)

---

## 1. Inputs — where to get them

Download these and place them as shown (paths are relative to the repo root; create `data/` yourself).

| Input | Source | Place at |
|---|---|---|
| **Bulk RNA-seq counts + metadata** (rat hippocampus, 60 samples) | NCBI GEO **[GSE329776](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329776)** (processed `normalized_filtered_counts.csv` + `sample_metadata.csv`) | `data/bulk_rnaseq/` |
| **Gene metadata** (Ensembl rat ID → symbol, Rnor_6.0) | `biomaRt::getBM()` against Ensembl rat, or the copy on the GEO/Zenodo record | `data/gene_metadata_with_rat_names.csv` |
| **Reference 1 — Allen 10x 2020** (Mouse Whole Cortex & Hippocampus) | [Allen Brain Map portal](https://celltypes.brain-map.org/) | `data/references/mouse10x_2020/` |
| **Reference 2 — Smart-seq2 (Tasic 2018)** | GEO **GSE115746** | `data/references/mouse_smartseq_2019/` |
| **Reference 3 — Yao et al. 2021** (Whole Cortex + HPF atlas) | NeMO archive (`mou_whole_cortex_subsample_10x`) | `data/references/yao_hippo_10x/` |

`sample_metadata.csv` is the manually-curated experimental design (animal age / region / parity) and
is the authoritative record, also deposited with the GEO submission.

---

## 2. How to run

**Setup** (once):
```sh
Rscript scripts/install_packages.R          # R deps (R 4.5.3 verified)
# SCDC is not on CRAN — install it explicitly:
Rscript -e 'remotes::install_github("renozao/xbioc"); remotes::install_github("meichendong/SCDC")'
pip install -r scripts/requirements.txt     # Python deps (reference filtering)
```

**Step A — references → deconvolution (Section 01).** Filter the three references by sex/region and
run SCDC to produce `all_results_scdc.rds` (the entry point for everything downstream):
```sh
# edit the input paths at the top of scripts/01_fig1_reference_sex/00_filter_*.py to your downloads
python scripts/01_fig1_reference_sex/00_filter_mouse10x_2020.py
python scripts/01_fig1_reference_sex/00_filter_smartseq_2019.py
python scripts/01_fig1_reference_sex/00_filter_yao_hippo.py
python scripts/01_fig1_reference_sex/00_validate_references.py
Rscript scripts/01_fig1_reference_sex/01_run_deconvolution.R   # -> all_results_scdc.rds
Rscript scripts/01_fig1_reference_sex/02_compute_metrics.R
Rscript scripts/01_fig1_reference_sex/03_generate_figure1.R    # Fig 1
```

**Step B — everything downstream of `all_results_scdc.rds` (Sections 02–09):**
```sh
Rscript scripts/run_all_robust.R            # all sections
Rscript scripts/run_all_robust.R 07         # one section (prefix match)
```
`run_all_robust.R` is the single entry point; it runs each section's canonical scripts in order and
prints a per-script OK/FAIL summary. Outputs land in `outputs/<section>/`. Set `REPO_ROOT=/path` to
run from outside the repo.

---

## 3. What builds each reported figure

Every reported panel, the script that produces it, and the output file:

| Figure | Script | Output |
|---|---|---|
| **Fig 1** reference/sex performance | `01_fig1_reference_sex/03_generate_figure1.R` | `Figure1_Sex_Reference_Performance.png` |
| **Fig 2a** per-dataset composition | `02_fig2_harmonization/06_regenerate_panel_a.R` | `Figure2_PanelA_per_dataset_sum.png` |
| **Fig 2b** correlation heatmap | `02_fig2_harmonization/01_harmonize_cell_types.R` | `correlation_heatmap.pdf` |
| **Fig 2c** merged composition | `02_fig2_harmonization/07_merge_composition.R` | `Figure2_PanelC_merge_composition.png` |
| **Fig 2d–e** UMAP | `02_fig2_harmonization/05_create_umap.R` | `final_harmonized_umap_by_{class,type}.pdf` |
| **Fig 3** age/region/interaction volcanos | `03_fig3_harmonized_anova/04_plot_log2fc_volcanos.R` | `{age,region}_effects_volcano.png`, `interaction_effects_plot.png` |
| **Fig 4a** pathway distribution / **4b** vulnerability | `05_fig4_age_region_enrichment/07_overview_panels_AB.R` | `Figure4_Panel{A,B}_*.png` |
| **Fig 4c–e** enrichment dotplots | `05_fig4_age_region_enrichment/06_dotplot_figure4.R` | `Figure4_combined_dotplots.png` |
| **Fig 5a** parity heatmap | `06_fig5ab_parity_proportions/06_parity_focused_heatmap.R` | `Figure5a_parity_focused_heatmap.png` |
| **Fig 5b** parity volcano | `06_fig5ab_parity_proportions/03_plot_parity_volcano.R` | `parity_volcano_unified.png` |
| **Fig 5c–e** cell-specific volcanos | `10_fig5_parity_enrichment/03e_themed_volcanos.R` | `volcano_themed_5{c,d,e}_*.png` |
| **Fig 5f–h** two-sided dotplots | `10_fig5_parity_enrichment/04_two_sided_dotplots.R` | `{CA3_Astro,CA3_Sst,Astro_Sst}_two_sided_dotplot.png` |
| **Fig 5 (Venn)** pathway overlap | `11_fig5_pathway_overlaps/01_pathway_venn_diagram.R` + `03_percentage_venns.R` | `pathway_venn_diagram.png` |
| **Fig 5i** shared-pathway dotplot | `11_fig5_pathway_overlaps/02_shared_pathways_dotplot.R` | `top20_shared_pathways_dotplot.png` |
| **Supp 1** gene-ANOVA volcanos (Age/Region/Parity + Age×Region + Age×Parity) | `04_suppfig1_gene_anova/02b_fc_volcano_plots.R` | `supp_fig1_{age,region,parity}_volcano.png`, `supp_fig1_age_{region,parity}_volcano.png` |
| **Supp 2** RF ∩ PCA Venn | `09_suppfig2to5_parity_genes/04_rf_pc_venn.R` | `supp_fig2_rf_pc_venn.png` |
| **Supp 3** DESeq2 binary-parity volcano | `09_suppfig2to5_parity_genes/07_binary_deseq_volcano.R` | `supp_fig3_binary_deseq_volcano.png` |
| **Supp 4** PCA drivers + PC×model contribution | `09_suppfig2to5_parity_genes/06_pc_drivers_plot.R` + `11_pc_model_contribution_heatmap.R` | `supp_fig4_pc_drivers.png`, `supp_fig4_pc_model_contribution.png` |
| **Supp 5** RF top-23 + CV accuracy + performance | `09_suppfig2to5_parity_genes/05_rf_top23_plot.R`, `08_supp_fig5_rf_panels.R`, `09_rf_binary_vs_multiclass.R` | `supp_fig5_rf_top23.png`, `supp_fig5_rf_cv_accuracy.png`, `supp_fig5_rf_binary_vs_multiclass.png` |

> Rendered figures are not stored in the repo (they're in the paper); the scripts regenerate them
> (named by their current Supp number), and the underlying results data live in `outputs/` and `supplement/`.

---

## 4. Repository layout

```
├── README.md
├── scripts/
│   ├── config.R                 central path config (set REPO_ROOT or REPRO_CONFIG)
│   ├── run_all_robust.R         single entry point (Sections 02–09)
│   ├── install_packages.R, requirements.txt
│   ├── 00_utilities/            helpers.R + adapter.R (sourced by section scripts)
│   └── 01_fig1_reference_sex/ … 11_fig5_pathway_overlaps/   one folder per figure group
├── outputs/                     results data per section (CSV/RDS; figures not stored — see §3)
├── supplement/                  paper-ready source data per figure
└── docs/Supplementary_Methods.md   canonical, parameter-complete methods
```

**Not included here** (by design): `data/` (download per §1), `checkpoints/` (regenerated by running
the pipeline), and the **rendered figures** (they're in the paper). These — plus the figures — live in
the complete local / Zenodo archive.

---

## 5. Reproduction status

Given the inputs above, the analysis regenerates the reported values. Verified totals:

- **Supp 1** gene ANOVA: 12,516 tested; 8,163 age / 6,086 region / 8 parity / 1,798 age×region / 7 age×parity significant (FDR < 0.05)
- **Fig 5a-b** parity-responsive cell types: CA3 (↑), Sst (↑), Astrocytes (↓)
- **Supp 4** PCA drivers: PC5 = 63, PC6 = 128, PC8 = 168 genes
- **Supp 5** RF: binary 71.7% / multi-class 58.3% OOB accuracy; top-23 genes (CV-accuracy peak at 23)
- **Fig 5f-i** pathway overlaps: 419 / 110 / 22 / 11 (of 1,600)

The published random-forest top-23 list loads from the saved model object for version-stable
reproduction; set `REPRO_RF_REFIT=TRUE` to refit from scratch.

## 6. Dependencies

- **R 4.5.3** (verified); Python 3 for reference filtering.
- SCDC (GitHub), DESeq2, clusterProfiler, ReactomePA, msigdbr, randomForest, tidyverse, ggplot2,
  ggrepel, patchwork, ComplexHeatmap, VennDiagram, viridis. See `scripts/install_packages.R`.

## License & citation

- Code: MIT (`LICENSE`). Prose/docs: CC-BY-4.0.
- Cite: McGovern AJ, Duarte-Guterman P, Galea LAM. *Cell Type-Specific Remodelling of the Rat
  Hippocampus by Parity and Age.* Journal of Neuroendocrinology (in press). Raw data: GEO **GSE329776**.
