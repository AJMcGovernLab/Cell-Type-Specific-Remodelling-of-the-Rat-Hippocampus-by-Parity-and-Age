# 03_compute_log2fc.R — REAL log2 fold changes for the 23 harmonized cell types.
#
# Replaces the fabricated `calculate_fold_changes()` in the original
# Manuscript_Figures/.../08_create_volcano_plots.R, which used runif() + grepl()
# to generate placeholder values. Here we compute log2(mean_groupA / mean_groupB)
# from the actual harmonized SCDC proportions for the three contrasts:
#   age:    log2(mean_Old / mean_Young)
#   region: log2(mean_Dorsal / mean_Ventral)
#   parity: log2(mean_Parous / mean_Nulliparous), where Parous = Primiparous + Biparous
#
# Outputs: outputs/03_fig3_harmonized_anova/harmonized_log2fc.csv
#          outputs/03_fig3_harmonized_anova/harmonized_group_means.csv

suppressPackageStartupMessages({ library(tidyverse) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

MIN_SAMPLES <- 60
MIN_VAR     <- 1e-8
PSEUDOCOUNT <- 1e-6   # avoid log2(0) when a group mean is zero

harmonize_one <- function(res, dataset_name, mapping, metadata) {
  ds_map <- mapping %>% filter(dataset == dataset_name)
  deconv <- if (!is.null(res$proportions))            res$proportions
            else if (!is.null(res$raw_result$prop.est.mvw)) {
              df <- as.data.frame(res$raw_result$prop.est.mvw); df$sample <- rownames(df); df
            } else return(NULL)

  cell_cols <- intersect(names(deconv), ds_map$original_name)
  if (length(cell_cols) == 0) return(NULL)

  deconv %>%
    select(sample, all_of(cell_cols)) %>%
    pivot_longer(-sample, names_to = "original_cell_type", values_to = "proportion") %>%
    left_join(ds_map %>% select(original_name, unified_name, cell_class),
              by = c("original_cell_type" = "original_name")) %>%
    group_by(sample, unified_name, cell_class) %>%
    summarise(proportion = sum(proportion), .groups = "drop") %>%
    left_join(metadata, by = "sample") %>%
    mutate(dataset = dataset_name)
}

safe_log2fc <- function(num, den) {
  log2((num + PSEUDOCOUNT) / (den + PSEUDOCOUNT))
}

main <- function() {
  out_dir   <- OUT$fig3
  mapping   <- read.csv(file.path(OUT$fig2, "refined_cell_type_mapping.csv"))
  metadata  <- read.csv(SAMPLE_META)

  female_rds <- c(mouse10x_2020       = "mouse10x_2020_female_scdc.rds",
                  mouse_smartseq_2019 = "mouse_smartseq_2019_female_scdc.rds",
                  yao_hippo_10x       = "yao_hippo_10x_female_scdc.rds")

  combined <- bind_rows(lapply(names(female_rds), function(name) {
    res <- readRDS(file.path(SCDC_CHECKPOINTS, female_rds[[name]]))
    harmonize_one(res, name, mapping, metadata)
  })) %>%
    mutate(parity_binary = ifelse(parity == "Nulliparous", "Nulliparous", "Parous"))

  candidates <- combined %>%
    group_by(unified_name) %>%
    summarise(n = n(), v = var(proportion), .groups = "drop") %>%
    filter(n >= MIN_SAMPLES, v > MIN_VAR) %>%
    pull(unified_name)

  message(sprintf("[log2fc] %d harmonized cell types pass filters", length(candidates)))

  group_means <- combined %>%
    filter(unified_name %in% candidates) %>%
    group_by(unified_name) %>%
    summarise(
      n_total       = n(),
      n_datasets    = n_distinct(dataset),
      grand_mean    = mean(proportion),
      mean_Young    = mean(proportion[age == "Young"]),
      mean_Old      = mean(proportion[age == "Old"]),
      mean_Dorsal   = mean(proportion[region == "Dorsal"]),
      mean_Ventral  = mean(proportion[region == "Ventral"]),
      mean_Nullip   = mean(proportion[parity_binary == "Nulliparous"]),
      mean_Parous   = mean(proportion[parity_binary == "Parous"]),
      n_Young       = sum(age    == "Young"),
      n_Old         = sum(age    == "Old"),
      n_Dorsal      = sum(region == "Dorsal"),
      n_Ventral     = sum(region == "Ventral"),
      n_Nullip      = sum(parity_binary == "Nulliparous"),
      n_Parous      = sum(parity_binary == "Parous"),
      .groups = "drop"
    )

  log2fc <- group_means %>%
    transmute(
      cell_type        = unified_name,
      n_total, n_datasets, grand_mean,
      log2fc_age       = safe_log2fc(mean_Old,    mean_Young),
      pct_change_age   = 100 * (mean_Old    - mean_Young)   / mean_Young,
      log2fc_region    = safe_log2fc(mean_Dorsal, mean_Ventral),
      pct_change_region= 100 * (mean_Dorsal - mean_Ventral) / mean_Ventral,
      log2fc_parity    = safe_log2fc(mean_Parous, mean_Nullip),
      pct_change_parity= 100 * (mean_Parous - mean_Nullip)  / mean_Nullip,
      mean_Young, mean_Old, mean_Dorsal, mean_Ventral, mean_Nullip, mean_Parous
    )

  anova_path <- file.path(out_dir, "harmonized_anova_results.csv")
  if (file.exists(anova_path)) {
    pvals <- read.csv(anova_path)
    log2fc <- log2fc %>%
      left_join(pvals %>% select(cell_type, age_pval, region_pval, parity_pval,
                                 age_region_pval, age_parity_pval,
                                 parity_region_pval, three_way_pval),
                by = "cell_type")
  }

  write.csv(log2fc, file.path(out_dir, "harmonized_log2fc.csv"), row.names = FALSE)
  write.csv(group_means, file.path(out_dir, "harmonized_group_means.csv"), row.names = FALSE)

  message(sprintf("[log2fc] wrote %s", file.path(out_dir, "harmonized_log2fc.csv")))
  message(sprintf("[log2fc] wrote %s", file.path(out_dir, "harmonized_group_means.csv")))

  invisible(log2fc)
}

if (sys.nframe() == 0L) main()
