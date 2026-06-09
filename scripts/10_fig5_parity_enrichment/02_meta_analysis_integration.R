# 02_meta_analysis_integration.R — canonical §3.4.2 pipeline.
#
# Per-author decision (2026-05-07): the manuscript's enrichment values were
# produced by an over-representation analysis (`enrichGO` / `enrichKEGG` via
# the v7 weighted meta-analysis script). We therefore use v7's staged output
# as the canonical source for §3.4.2 rather than the alternative curated
# GSEA pipeline.
#
# v7 columns:
#   weighted_mean_NES   — actually log2(GeneRatio / BgRatio); rename for honesty
#   weighted_meta_padj  — Fisher-combined FDR across methods
#   meta_group          — CA3_do_combined / 376_Astro / 78_Sst_HPF
#
# This script reads the staged v7 file and produces:
#   parity_high_confidence_pathways.csv  — Moderate+ tier filtered table
#   {CELL}_meta_analysis.csv             — per-cell pathway list (Fig 5c-e source)
#
# Direction tagging (positive vs negative log2 fold enrichment) is added by
# 02c_tag_pathway_direction.R using mean DESeq2 log2FC of pathway members.

suppressPackageStartupMessages({ library(tidyverse) })
source(Sys.getenv(
  "REPRO_CONFIG",
  "f:/Parity/Final/Repository/scripts/config.R"
))

V7_CKPT <- file.path(
  CHECKPOINT_DIR, "enrichment_parity",
  "functional_enrichment_v7", "weighted_meta_analysis_significant_v7.csv"
)

main <- function() {
  out_dir <- OUT$fig5ce
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(V7_CKPT)) {
    stop(sprintf("v7 staged output not found: %s", V7_CKPT))
  }

  v7 <- read.csv(V7_CKPT, stringsAsFactors = FALSE)
  message(sprintf("[meta] read %d v7 meta-rows", nrow(v7)))

  # Rename `weighted_mean_NES` to `log2_fold_enrichment` for honesty.
  # Keep the v7 column under its original name for legacy compatibility.
  v7 <- v7 %>%
    mutate(
      log2_fold_enrichment = weighted_mean_NES,  # honest name
      mean_nes             = weighted_mean_NES,  # legacy alias for downstream code
      meta_fdr             = weighted_meta_padj,
      meta_pvalue          = weighted_meta_pvalue,
      meta_score           = weighted_enhanced_score,
      n_methods            = n_studies,
      cell_type            = meta_group,
      effect               = "parity",
      regulation_pattern_set_membership = case_when(
        log2_fold_enrichment > 0 ~ "Over-represented",
        log2_fold_enrichment < 0 ~ "Under-represented",
        TRUE                     ~ "Mixed"
      ),
      confidence_tier = case_when(
        n_methods >= 4 & meta_score > 7 ~ "Ultra-High",
        n_methods >= 3 & meta_score > 5 ~ "High",
        n_methods >= 2 & meta_score > 3 ~ "Moderate",
        TRUE                            ~ "Method-Specific"
      )
    )

  # High-confidence (Moderate+) table — paper cites top entries.
  hc <- v7 %>%
    filter(confidence_tier %in% c("Ultra-High", "High", "Moderate")) %>%
    arrange(desc(meta_score))
  write.csv(hc, file.path(out_dir, "parity_high_confidence_pathways.csv"),
            row.names = FALSE)
  message(sprintf("[meta] wrote %d high-confidence pathways", nrow(hc)))

  # Per-cell pathway lists (Fig 5c-e source).
  for (ct in unique(v7$cell_type)) {
    slug <- gsub("[^A-Za-z0-9]", "_", ct)
    rows <- v7 %>% filter(cell_type == !!ct) %>% arrange(desc(meta_score))
    write.csv(rows,
              file.path(out_dir, sprintf("%s_meta_analysis.csv", slug)),
              row.names = FALSE)
    message(sprintf("[meta]   %-25s %d pathways", ct, nrow(rows)))
  }

  invisible(v7)
}

if (sys.nframe() == 0L) main()
