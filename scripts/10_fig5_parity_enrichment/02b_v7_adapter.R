# 02b_v7_adapter.R — emit a v7-shaped CSV from the curated meta-analysis
# so legacy visualization scripts (03/04/05 in this section, 01/02/03 in
# section 09) can read it without rewriting their column references.
#
# This is a *compatibility shim only*. The canonical output is
# parity_high_confidence_pathways.csv produced by 02_meta_analysis_integration.R.
# This adapter renames columns to match the v7 schema so legacy plotting
# scripts (which reference `weighted_mean_NES` and `weighted_meta_padj_corrected`)
# can stay analysis-logic-unchanged.

suppressPackageStartupMessages({ library(tidyverse) })
source(Sys.getenv(
  "REPRO_CONFIG",
  "f:/Parity/Final/Repository/scripts/config.R"
))

main <- function() {
  curated <- read.csv(
    file.path(OUT$fig5ce, "parity_high_confidence_pathways.csv"),
    stringsAsFactors = FALSE
  )

  v7_shape <- curated %>%
    rename(
      weighted_mean_NES          = mean_nes,
      weighted_meta_padj_corrected = meta_fdr,
      weighted_meta_pvalue       = meta_pvalue,
      weighted_enhanced_score    = consensus_score
    ) %>%
    mutate(
      meta_group = case_when(
        cell_type %in% c("358_CA3_do", "356_CA3_do") ~ "CA3_do_combined",
        cell_type == "376_Astro"                     ~ "376_Astro",
        cell_type == "78_Sst_HPF"                    ~ "78_Sst_HPF",
        TRUE                                         ~ cell_type
      ),
      percentage = "5pct",  # legacy column required by some scripts
      n_cell_types = ifelse(is.na(n_methods), 1L, n_methods),
      analysis_version = "curated_GSEA_via_v7_shim"
    )

  out_file <- file.path(
    OUT$fig5ce,
    "weighted_meta_analysis_significant_v7.csv"
  )
  write.csv(v7_shape, out_file, row.names = FALSE)
  message(sprintf("[v7-adapter] wrote %s (%d rows)",
                  out_file, nrow(v7_shape)))
  invisible(v7_shape)
}

if (sys.nframe() == 0L) main()
