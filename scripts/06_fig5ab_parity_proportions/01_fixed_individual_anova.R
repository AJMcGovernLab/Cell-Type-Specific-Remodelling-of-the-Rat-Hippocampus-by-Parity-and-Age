# 01_fixed_individual_anova.R — per-dataset three-way ANOVA on unharmonized cell types.
# Condensed from Set 1/0Full Analysis/Step5_NoHarmony_Effects/fixed_individual_anova.R (325 lines).

suppressPackageStartupMessages({ library(tidyverse); library(broom) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

MIN_VAR   <- 1e-8
MIN_MEAN  <- 1e-4
MIN_NONZERO <- 30

anova_one <- function(ct, deconv, dataset_name, metadata) {
  d <- deconv %>%
    select(sample, proportion = all_of(ct)) %>%
    left_join(metadata, by = "sample") %>%
    mutate(age    = factor(age,    levels = c("Young","Old")),
           parity = factor(parity, levels = c("Nulliparous","Primiparous","Biparous")),
           region = factor(region, levels = c("Dorsal","Ventral")))
  tryCatch({
    a <- tidy(aov(proportion ~ age * parity * region, data = d))
    data.frame(
      cell_type          = ct,
      dataset            = dataset_name,
      mean_proportion    = mean(d$proportion),
      age_pval           = a$p.value[a$term == "age"],
      parity_pval        = a$p.value[a$term == "parity"],
      region_pval        = a$p.value[a$term == "region"],
      age_parity_pval    = a$p.value[a$term == "age:parity"],
      age_region_pval    = a$p.value[a$term == "age:region"],
      parity_region_pval = a$p.value[a$term == "parity:region"],
      three_way_pval     = a$p.value[a$term == "age:parity:region"],
      n_samples          = nrow(d)
    )
  }, error = function(e) NULL)
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
  rows <- lapply(cell_cols, anova_one, deconv = deconv,
                 dataset_name = dataset_name, metadata = metadata)
  bind_rows(rows)
}

main <- function() {
  out_dir <- OUT$fig5ab
  metadata <- read.csv(SAMPLE_META)

  female_rds <- c(mouse10x_2020       = "mouse10x_2020_female_scdc.rds",
                  mouse_smartseq_2019 = "mouse_smartseq_2019_female_scdc.rds",
                  yao_hippo_10x       = "yao_hippo_10x_female_scdc.rds")

  combined <- bind_rows(lapply(names(female_rds), function(name) {
    res <- readRDS(file.path(SCDC_CHECKPOINTS, female_rds[[name]]))
    analyze_dataset(res, name, metadata)
  }))

  write.csv(combined,
            file.path(out_dir, "fixed_individual_anova_results.csv"),
            row.names = FALSE)

  responsive <- combined %>%
    filter(parity_pval < CONFIG$fdr_threshold) %>%
    arrange(parity_pval)
  write.csv(responsive,
            file.path(out_dir, "parity_responsive_cell_types.csv"),
            row.names = FALSE)

  message(sprintf("[parity-anova] %d rows; %d parity-responsive (p < %.2f)",
                  nrow(combined), nrow(responsive), CONFIG$fdr_threshold))
  invisible(combined)
}

if (sys.nframe() == 0L) main()
