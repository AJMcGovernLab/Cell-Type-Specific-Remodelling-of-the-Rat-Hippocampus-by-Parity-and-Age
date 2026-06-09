# 04_create_factor_summaries.R — Figure 4 source CSVs.
# Condensed from Paper_Level/6_Cell_Specific_Functional_Enrichment_Analysis/
# comprehensive_weighted_enrichment/publication_materials/create_factor_summaries.R
# (290 lines → ~90).
#
# Input: the paper's pre-computed meta-analysis table
#   Repository/checkpoints/enrichment_fig4/all_meta_analysis_results.csv
#   Repository/checkpoints/enrichment_fig4/cell_type_summary.csv

suppressPackageStartupMessages({ library(tidyverse) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

FIG4_CKPT <- file.path(CHECKPOINT_DIR, "enrichment_fig4")
HIGH_CONF <- c("Ultra-High", "High", "Moderate")
EFFECT_LABELS <- c(age = "age", region = "region", interaction = "age_region")

main <- function() {
  out_dir <- OUT$fig4
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  meta <- read.csv(file.path(FIG4_CKPT, "all_meta_analysis_results.csv"),
                   stringsAsFactors = FALSE)
  cells <- read.csv(file.path(FIG4_CKPT, "cell_type_summary.csv"),
                    stringsAsFactors = FALSE)
  message(sprintf("[fig4] %d meta-analysis rows, %d cell-type/effect combos",
                  nrow(meta), nrow(cells)))

  # Per-effect high-confidence pathway tables
  for (lab in names(EFFECT_LABELS)) {
    eff <- EFFECT_LABELS[[lab]]
    df <- meta |>
      filter(effect == eff, confidence_tier %in% HIGH_CONF) |>
      arrange(desc(meta_score)) |>
      select(cell_type, description, meta_score, confidence_tier,
             n_methods, database)
    write.csv(df,
              file.path(out_dir,
                        sprintf("%s_high_confidence_pathways.csv", lab)),
              row.names = FALSE)
    message(sprintf("[fig4] %-12s %d high-conf pathways", lab, nrow(df)))
  }

  # Vulnerability ranking — Fig 4a/b source
  cell_vuln <- cells |>
    mutate(high_conf_total     = n_ultra_high + n_high + n_moderate,
           vulnerability_score = mean_meta_score * (high_conf_total + 1)) |>
    arrange(desc(vulnerability_score)) |>
    select(effect, cell_type, mean_meta_score, high_conf_total,
           vulnerability_score, top_pathway, top_meta_score)
  write.csv(cell_vuln,
            file.path(out_dir, "cell_type_vulnerability_ranking.csv"),
            row.names = FALSE)

  # Cross-factor summary
  invisible(list(meta = meta, vuln = cell_vuln))
}

if (sys.nframe() == 0L) main()
