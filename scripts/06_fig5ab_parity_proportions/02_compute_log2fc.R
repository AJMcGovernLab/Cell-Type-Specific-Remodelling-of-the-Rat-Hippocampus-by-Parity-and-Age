# 02_compute_log2fc.R — REAL log2 fold changes for unharmonized per-dataset cell types.
#
# Replaces the fabricated `calculate_parity_fold_changes()` in the original
# 06_Parity_Effect_Volcano_Plots/scripts/01_create_parity_volcano_plots.R, which
# used runif() bins + grepl(cell_type) to assign signs.
#
# We compute, per (dataset × cell_type), the log2 fold change for three contrasts:
#   age:    log2(mean_Old / mean_Young)
#   region: log2(mean_Dorsal / mean_Ventral)
#   parity: log2(mean_Parous / mean_Nulliparous)   where Parous = Primiparous + Biparous
#
# Outputs: outputs/06_fig5ab_parity_proportions/per_dataset_log2fc.csv  (one row per cell × dataset)
#          outputs/06_fig5ab_parity_proportions/parity_responsive_log2fc.csv  (cells with parity p<0.05)

suppressPackageStartupMessages({ library(tidyverse) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

MIN_VAR     <- 1e-8
MIN_MEAN    <- 1e-4
MIN_NONZERO <- 30
PSEUDOCOUNT <- 1e-6

safe_log2fc <- function(num, den) log2((num + PSEUDOCOUNT) / (den + PSEUDOCOUNT))

per_cell <- function(ct, deconv, dataset_name, metadata) {
  d <- deconv %>%
    select(sample, proportion = all_of(ct)) %>%
    left_join(metadata, by = "sample") %>%
    mutate(parity_binary = ifelse(parity == "Nulliparous", "Nulliparous", "Parous"))

  data.frame(
    cell_type    = ct,
    dataset      = dataset_name,
    n_samples    = nrow(d),
    grand_mean   = mean(d$proportion),
    mean_Young   = mean(d$proportion[d$age    == "Young"]),
    mean_Old     = mean(d$proportion[d$age    == "Old"]),
    mean_Dorsal  = mean(d$proportion[d$region == "Dorsal"]),
    mean_Ventral = mean(d$proportion[d$region == "Ventral"]),
    mean_Nullip  = mean(d$proportion[d$parity_binary == "Nulliparous"]),
    mean_Parous  = mean(d$proportion[d$parity_binary == "Parous"])
  ) %>%
    mutate(
      log2fc_age        = safe_log2fc(mean_Old,    mean_Young),
      pct_change_age    = 100 * (mean_Old    - mean_Young)   / mean_Young,
      log2fc_region     = safe_log2fc(mean_Dorsal, mean_Ventral),
      pct_change_region = 100 * (mean_Dorsal - mean_Ventral) / mean_Ventral,
      log2fc_parity     = safe_log2fc(mean_Parous, mean_Nullip),
      pct_change_parity = 100 * (mean_Parous - mean_Nullip)  / mean_Nullip
    )
}

analyze_dataset <- function(res, dataset_name, metadata) {
  deconv <- if (!is.null(res$proportions))            res$proportions
            else if (!is.null(res$raw_result$prop.est.mvw)) {
              df <- as.data.frame(res$raw_result$prop.est.mvw); df$sample <- rownames(df); df
            } else return(NULL)

  cell_cols <- setdiff(names(deconv), c("sample","sample_id"))
  ok <- vapply(cell_cols, function(ct) {
    x <- deconv[[ct]]
    var(x) > MIN_VAR && mean(x) > MIN_MEAN && sum(x > 1e-6) >= MIN_NONZERO
  }, logical(1))
  cell_cols <- cell_cols[ok]

  message(sprintf("[%s] %d cell types pass filters", dataset_name, length(cell_cols)))
  rows <- lapply(cell_cols, per_cell, deconv = deconv,
                 dataset_name = dataset_name, metadata = metadata)
  bind_rows(rows)
}

main <- function() {
  out_dir  <- OUT$fig5ab
  metadata <- read.csv(SAMPLE_META)

  female_rds <- c(mouse10x_2020       = "mouse10x_2020_female_scdc.rds",
                  mouse_smartseq_2019 = "mouse_smartseq_2019_female_scdc.rds",
                  yao_hippo_10x       = "yao_hippo_10x_female_scdc.rds")

  combined <- bind_rows(lapply(names(female_rds), function(name) {
    res <- readRDS(file.path(SCDC_CHECKPOINTS, female_rds[[name]]))
    analyze_dataset(res, name, metadata)
  }))

  anova_path <- file.path(out_dir, "fixed_individual_anova_results.csv")
  if (file.exists(anova_path)) {
    pvals <- read.csv(anova_path)
    combined <- combined %>%
      left_join(pvals %>% select(cell_type, dataset,
                                 age_pval, parity_pval, region_pval,
                                 age_parity_pval, age_region_pval,
                                 parity_region_pval, three_way_pval),
                by = c("cell_type", "dataset"))
  }

  write.csv(combined, file.path(out_dir, "per_dataset_log2fc.csv"), row.names = FALSE)

  responsive <- combined %>%
    filter(!is.na(parity_pval) & parity_pval < CONFIG$fdr_threshold) %>%
    arrange(parity_pval)
  write.csv(responsive,
            file.path(out_dir, "parity_responsive_log2fc.csv"),
            row.names = FALSE)

  message(sprintf("[log2fc] %d rows; %d parity-responsive (p < %.2f)",
                  nrow(combined), nrow(responsive), CONFIG$fdr_threshold))
  invisible(combined)
}

if (sys.nframe() == 0L) main()
