# 04_within_between_correlations.R — computes within/between cluster correlation statistics.
# Produces §3.2's "r = 0.78 ± 0.12 within, 0.23 ± 0.18 between, Wilcoxon p < 0.001"
# numbers, which were not written to disk in the original pipeline.

suppressPackageStartupMessages(library(tidyverse))
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

main <- function() {
  out_dir <- OUT$fig2

  mapping  <- readRDS(file.path(out_dir, "mapping_table.rds"))    # from 01_
  combined <- readRDS(file.path(out_dir, "combined_correlation_data.rds"))
  cor_mat  <- combined$cor_matrix
  names    <- combined$cell_type_names

  cluster_of <- setNames(mapping$cluster_id[match(names, mapping$original_name)],
                         names)

  pairs <- expand.grid(i = seq_along(names), j = seq_along(names))
  pairs <- pairs[pairs$i < pairs$j, ]
  pairs$r      <- cor_mat[cbind(pairs$i, pairs$j)]
  pairs$within <- cluster_of[names[pairs$i]] == cluster_of[names[pairs$j]]

  within  <- pairs$r[pairs$within]
  between <- pairs$r[!pairs$within]

  wt <- wilcox.test(within, between, alternative = "greater")

  stats <- data.frame(
    statistic        = c("mean_within", "sd_within",
                         "mean_between", "sd_between",
                         "wilcoxon_W", "wilcoxon_p_value",
                         "n_within_pairs", "n_between_pairs"),
    value            = c(mean(within),  sd(within),
                         mean(between), sd(between),
                         unname(wt$statistic), wt$p.value,
                         length(within), length(between))
  )
  write.csv(stats, file.path(out_dir, "correlation_summary_stats.csv"), row.names = FALSE)
  message(sprintf("[corr] within %.3f ± %.3f | between %.3f ± %.3f | Wilcoxon p = %.3e",
                  mean(within), sd(within), mean(between), sd(between), wt$p.value))
}

if (sys.nframe() == 0L) main()
