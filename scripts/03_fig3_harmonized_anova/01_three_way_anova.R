# 01_three_way_anova.R — three-way ANOVA on the 23 harmonized cell-type proportions.
# Condensed from Set 1/0Full Analysis/Step3_Harmonization/Statistical_Results/harmonized_three_way_anova.R (~360 lines).

suppressPackageStartupMessages({
  library(tidyverse); library(broom); library(lme4); library(lmerTest)
})
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

MIN_SAMPLES <- 60
MIN_VAR     <- 1e-8

# ---- harmonize one dataset's proportions onto unified_name -----------------

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

# ---- per-cell-type model --------------------------------------------------

run_anova_one <- function(data, cell_type) {
  d <- data %>% filter(unified_name == cell_type)
  if (nrow(d) < 20 || var(d$proportion) < MIN_VAR) return(NULL)

  out <- list(cell_type = cell_type,
              n_datasets = n_distinct(d$dataset),
              total_samples = nrow(d))

  tryCatch({
    out$fixed_anova <- tidy(aov(proportion ~ age * parity * region, data = d))
    if (n_distinct(d$dataset) > 1) {
      m <- lmer(proportion ~ age * parity * region + (1 | dataset),
                data = d, control = lmerControl(optimizer = "bobyqa"))
      out$mixed_anova <- anova(m)
      out$model_type  <- "mixed"
    } else {
      out$model_type <- "fixed"
    }
    out$summary_stats <- d %>%
      group_by(age, parity, region) %>%
      summarise(mean_prop = mean(proportion),
                se_prop   = sd(proportion) / sqrt(n()),
                n = n(), .groups = "drop")
  }, error = function(e) message(sprintf("  [skip] %s: %s", cell_type, e$message)))
  out
}

# ---- extract one row of p-values per cell type ----------------------------

extract_effects <- function(anova_results) {
  lapply(names(anova_results), function(ct) {
    r <- anova_results[[ct]]
    if (r$model_type == "mixed" && !is.null(r$mixed_anova)) {
      a <- r$mixed_anova
      data.frame(cell_type = ct, n_datasets = r$n_datasets,
                 age_pval            = a["age",                "Pr(>F)"],
                 parity_pval         = a["parity",             "Pr(>F)"],
                 region_pval         = a["region",             "Pr(>F)"],
                 age_parity_pval     = a["age:parity",         "Pr(>F)"],
                 age_region_pval     = a["age:region",         "Pr(>F)"],
                 parity_region_pval  = a["parity:region",      "Pr(>F)"],
                 three_way_pval      = a["age:parity:region",  "Pr(>F)"])
    } else if (!is.null(r$fixed_anova)) {
      a <- r$fixed_anova
      data.frame(cell_type = ct, n_datasets = r$n_datasets,
                 age_pval            = a$p.value[a$term == "age"],
                 parity_pval         = a$p.value[a$term == "parity"],
                 region_pval         = a$p.value[a$term == "region"],
                 age_parity_pval     = a$p.value[a$term == "age:parity"],
                 age_region_pval     = a$p.value[a$term == "age:region"],
                 parity_region_pval  = a$p.value[a$term == "parity:region"],
                 three_way_pval      = a$p.value[a$term == "age:parity:region"])
    }
  }) %>% bind_rows()
}

# ---- main -----------------------------------------------------------------

main <- function() {
  out_dir   <- OUT$fig3
  mapping   <- read.csv(file.path(OUT$fig2, "refined_cell_type_mapping.csv"))
  metadata  <- read.csv(SAMPLE_META)

  female_rds <- c(mouse10x_2020       = "mouse10x_2020_female_scdc.rds",
                  mouse_smartseq_2019 = "mouse_smartseq_2019_female_scdc.rds",
                  yao_hippo_10x       = "yao_hippo_10x_female_scdc.rds")

  harmonized_list <- lapply(names(female_rds), function(name) {
    res <- readRDS(file.path(SCDC_CHECKPOINTS, female_rds[[name]]))
    harmonize_one(res, name, mapping, metadata)
  })
  combined <- bind_rows(harmonized_list) %>%
    mutate(
      age    = factor(age,    levels = c("Young", "Old")),
      parity = factor(parity, levels = c("Nulliparous", "Primiparous", "Biparous")),
      region = factor(region, levels = c("Dorsal",  "Ventral"))
    )

  candidates <- combined %>%
    group_by(unified_name) %>%
    summarise(n = n(), v = var(proportion), .groups = "drop") %>%
    filter(n >= MIN_SAMPLES, v > MIN_VAR) %>%
    pull(unified_name)

  message(sprintf("[anova] testing %d cell types", length(candidates)))
  anova_results <- setNames(
    lapply(candidates, function(ct) run_anova_one(combined, ct)),
    candidates
  )
  anova_results <- anova_results[!sapply(anova_results, is.null)]

  effects_df <- extract_effects(anova_results)
  write.csv(effects_df, file.path(out_dir, "harmonized_anova_results.csv"),
            row.names = FALSE)
  message(sprintf("[anova] wrote %s",
                  file.path(out_dir, "harmonized_anova_results.csv")))
  invisible(effects_df)
}

if (sys.nframe() == 0L) main()
