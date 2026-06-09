# run_all_robust.R — end-to-end reproducer (canonical scripts only).
#
# Runs each section's pipeline scripts in dependency order via system2(Rscript)
# from the section's own directory, so each behaves as if you had cd'd into the
# folder and run `Rscript <script>`.
#
# PREREQUISITE — Section 01 (reference filtering + SCDC deconvolution) is NOT run
# here: it needs the three public scRNA-seq references downloaded (see README
# "Inputs") and the SCDC package (not on CRAN). Run Section 01 first to produce
# `all_results_scdc.rds` and the reference data, then this script reproduces
# Sections 02-09 (all reported main + supplementary figures).
#
# Usage:  Rscript scripts/run_all_robust.R          # sections 02-09
#         Rscript scripts/run_all_robust.R 07        # only sections starting "07"

source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

Rscript <- file.path(R.home("bin"),
                     if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
Sys.setenv(REPRO_CONFIG = file.path(REPO_ROOT, "scripts", "config.R"))

# Canonical script per section (the producer of each reported figure/data).
section_scripts <- list(
  "02_fig2_harmonization" = c(
    "01_harmonize_cell_types.R", "02_refine_mapping.R", "03_count_201_to_23.R",
    "04_within_between_correlations.R", "05_create_umap.R",
    "06_regenerate_panel_a.R", "07_merge_composition.R"),            # Fig 2
  "03_fig3_harmonized_anova" = c(
    "01_three_way_anova.R", "03_compute_log2fc.R",
    "04_plot_log2fc_volcanos.R", "05_marginal_age_region_log2fc.R"), # Fig 3
  "04_suppfig1_gene_anova" = c(
    "01_gene_level_anova.R", "02b_fc_volcano_plots.R"),              # Supp Fig 1
  "05_fig4_age_region_enrichment" = c(
    "01_prepare_gene_sets.R", "02_run_weighted_enrichment.R", "03_meta_analysis.R",
    "04_create_factor_summaries.R", "06_dotplot_figure4.R",
    "07_overview_panels_AB.R"),                                      # Fig 4
  "06_fig5ab_parity_proportions" = c(
    "01_fixed_individual_anova.R", "02_compute_log2fc.R", "03_plot_parity_volcano.R",
    "05_marginal_parity_log2fc.R", "06_parity_focused_heatmap.R"),   # Fig 5a-b
  "07_pca_preprocessing" = c("01_pca_and_factor_association.R"),
  "08_rf_preprocessing"  = c("01_residualize_and_prep_ml.R"),
  "09_suppfig2to5_parity_genes" = c(
    "01_deseq2_binary_parity.R", "02_random_forest.R", "03_pc_driver_analysis.R",
    "04_rf_pc_venn.R", "05_rf_top23_plot.R", "06_pc_drivers_plot.R",
    "07_binary_deseq_volcano.R", "08_supp_fig5_rf_panels.R",
    "09_rf_binary_vs_multiclass.R", "11_pc_model_contribution_heatmap.R"), # Supp 2-5
  "10_fig5_parity_enrichment" = c(
    "01_parity_weighted_enrichment.R", "02_meta_analysis_integration.R",
    "02b_v7_adapter.R", "03e_themed_volcanos.R",
    "04_two_sided_dotplots.R", "05_cell_type_specificity.R"),        # Fig 5c-h
  "11_fig5_pathway_overlaps" = c(
    "01_pathway_overlap_summary.R", "01_pathway_venn_diagram.R",
    "02_shared_pathways_dotplot.R", "03_percentage_venns.R")         # Fig 5f-i
)

args <- commandArgs(trailingOnly = TRUE)
sections <- if (length(args)) {
  unique(unlist(lapply(args, function(a)
    grep(paste0("^", a), names(section_scripts), value = TRUE))))
} else names(section_scripts)

results <- list()
for (sec in sections) {
  sec_dir <- file.path(REPO_ROOT, "scripts", sec)
  message(sprintf("\n==== %s ====", sec))
  for (scr in section_scripts[[sec]]) {
    scr_path <- file.path(sec_dir, scr)
    if (!file.exists(scr_path)) { message(sprintf("  [skip] %s (not found)", scr)); next }
    message(sprintf("  -> %s", scr))
    t0 <- Sys.time(); log_out <- tempfile(fileext = ".log")
    rc <- system2(Rscript, args = shQuote(scr_path), stdout = log_out, stderr = log_out)
    results[[paste(sec, scr, sep = "/")]] <-
      list(rc = rc, dt = as.numeric(Sys.time() - t0, units = "secs"), log = log_out)
    if (rc != 0) {
      message(sprintf("     FAILED (rc=%d). Last 8 log lines:", rc))
      for (l in tryCatch(tail(readLines(log_out, warn = FALSE), 8),
                         error = function(e) character(0))) message(sprintf("     | %s", l))
    }
  }
}

cat("\n========== SUMMARY ==========\n")
n_ok <- sum(vapply(results, function(r) r$rc == 0, logical(1)))
for (key in names(results)) {
  r <- results[[key]]
  if (r$rc == 0) cat(sprintf("  OK    %-62s (%.1f s)\n", key, r$dt))
  else cat(sprintf("  FAIL  %-62s rc=%d (%.1f s) log=%s\n", key, r$rc, r$dt, r$log))
}
cat(sprintf("\n%d OK, %d FAIL out of %d scripts\n", n_ok, length(results) - n_ok, length(results)))
