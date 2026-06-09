# 05_marginal_age_region_log2fc.R — properly-conditioned age and region log2FC.
#
# For each harmonized cell type, computes:
#   - age main effect (Old vs Young), controlling for parity and region
#   - region main effect (Dorsal vs Ventral), controlling for parity and age
# Three estimators, all should agree:
#   (a) EMM (predicted-marginal-mean ratio from additive lm)
#   (b) lm coefficient (additive model: prop ~ age + parity_binary + region)
#   (c) stratum-averaged log2FC: per-stratum log2FC, then averaged
# Plus per-stratum sign count (out of 4 strata).
#
# Output: outputs/03_fig3_harmonized_anova/marginal_age_region_log2fc.csv

suppressPackageStartupMessages({ library(tidyverse) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

PSEUDOCOUNT <- 1e-6
safe_log2fc <- function(num, den) log2((num + PSEUDOCOUNT) / (den + PSEUDOCOUNT))

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

per_cell_age_region <- function(d_cell) {
  d <- d_cell %>%
    mutate(parity_binary = factor(ifelse(parity == "Nulliparous", "Nulliparous", "Parous"),
                                  levels = c("Nulliparous","Parous")),
           age    = factor(age,    levels = c("Young","Old")),
           region = factor(region, levels = c("Dorsal","Ventral")))
  if (var(d$proportion) < 1e-10) return(NULL)

  fit <- tryCatch(lm(proportion ~ age + parity_binary + region, data = d),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  cf <- coef(fit)

  # ---- AGE: stratify by parity_binary × region (4 strata)
  age_strata <- d %>%
    group_by(parity_binary, region, age) %>%
    summarise(mu = mean(proportion), .groups = "drop") %>%
    pivot_wider(names_from = age, values_from = mu) %>%
    mutate(log2fc = safe_log2fc(Old, Young))
  age_stratum_mean <- mean(age_strata$log2fc, na.rm = TRUE)
  age_strata_old_higher <- sum(age_strata$Old > age_strata$Young, na.rm = TRUE)

  age_grid <- expand.grid(age = c("Young","Old"),
                          parity_binary = c("Nulliparous","Parous"),
                          region = c("Dorsal","Ventral"))
  age_grid$pred <- predict(fit, newdata = age_grid)
  age_emm <- age_grid %>% group_by(age) %>% summarise(mu = mean(pred), .groups = "drop")
  age_emm_y <- age_emm$mu[age_emm$age == "Young"]
  age_emm_o <- age_emm$mu[age_emm$age == "Old"]
  age_emm_log2fc <- safe_log2fc(age_emm_o, age_emm_y)
  age_pct <- 100 * (age_emm_o - age_emm_y) / age_emm_y
  age_lm_log2fc <- safe_log2fc(cf["(Intercept)"] + cf["ageOld"], cf["(Intercept)"])

  # ---- REGION: stratify by parity_binary × age (4 strata)
  region_strata <- d %>%
    group_by(parity_binary, age, region) %>%
    summarise(mu = mean(proportion), .groups = "drop") %>%
    pivot_wider(names_from = region, values_from = mu) %>%
    mutate(log2fc = safe_log2fc(Dorsal, Ventral))
  region_stratum_mean <- mean(region_strata$log2fc, na.rm = TRUE)
  region_strata_dorsal_higher <- sum(region_strata$Dorsal > region_strata$Ventral, na.rm = TRUE)

  region_grid <- age_grid
  region_emm <- region_grid %>% group_by(region) %>% summarise(mu = mean(pred), .groups = "drop")
  region_emm_d <- region_emm$mu[region_emm$region == "Dorsal"]
  region_emm_v <- region_emm$mu[region_emm$region == "Ventral"]
  region_emm_log2fc <- safe_log2fc(region_emm_d, region_emm_v)
  region_pct <- 100 * (region_emm_d - region_emm_v) / region_emm_v
  # region coefficient (`regionVentral`) is Ventral relative to Dorsal (the
  # factor reference level). Compute Dorsal-vs-Ventral log2FC accordingly.
  region_lm_log2fc <- safe_log2fc(cf["(Intercept)"],
                                  cf["(Intercept)"] + cf["regionVentral"])

  data.frame(
    age_emm_log2fc      = age_emm_log2fc,
    age_pct_change      = age_pct,
    age_lm_log2fc       = age_lm_log2fc,
    age_stratum_log2fc  = age_stratum_mean,
    age_strata_old_higher = age_strata_old_higher,   # of 4
    age_three_agree     = (sign(age_emm_log2fc) == sign(age_lm_log2fc) &
                           sign(age_emm_log2fc) == sign(age_stratum_mean)),
    region_emm_log2fc      = region_emm_log2fc,
    region_pct_change      = region_pct,
    region_lm_log2fc       = region_lm_log2fc,
    region_stratum_log2fc  = region_stratum_mean,
    region_strata_dorsal_higher = region_strata_dorsal_higher,  # of 4
    region_three_agree     = (sign(region_emm_log2fc) == sign(region_lm_log2fc) &
                              sign(region_emm_log2fc) == sign(region_stratum_mean))
  )
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
  }))

  # use cells with n >= 60 (matches the published ANOVA filter)
  candidates <- combined %>%
    group_by(unified_name) %>%
    summarise(n = n(), v = var(proportion), .groups = "drop") %>%
    filter(n >= 60, v > 1e-10) %>%
    pull(unified_name)

  results <- bind_rows(lapply(candidates, function(ct) {
    out <- per_cell_age_region(combined %>% filter(unified_name == ct))
    if (is.null(out)) return(NULL)
    cbind(cell_type = ct, out)
  }))

  # join the published ANOVA p-values
  anova_path <- file.path(out_dir, "harmonized_anova_results.csv")
  if (file.exists(anova_path)) {
    pvals <- read.csv(anova_path) %>%
      select(cell_type, age_pval, region_pval, age_region_pval)
    results <- results %>% left_join(pvals, by = "cell_type")
  }

  results <- results %>%
    arrange(age_pval) %>%
    mutate(
      age_significant    = age_pval    < 0.05,
      region_significant = region_pval < 0.05
    )

  write.csv(results, file.path(out_dir, "marginal_age_region_log2fc.csv"), row.names = FALSE)

  cat("\n========== AGE main effect (sig only, p<0.05) ==========\n\n")
  age_sig <- results %>% filter(age_significant)
  for (i in seq_len(nrow(age_sig))) {
    r <- age_sig[i, ]
    cat(sprintf("%s   age_pval = %.3g\n", r$cell_type, r$age_pval))
    cat(sprintf("  EMM log2FC = %+.3f  (%+.1f%%)   lm = %+.3f   stratum-avg = %+.3f   Old higher in %d/4 strata   agree=%s\n\n",
                r$age_emm_log2fc, r$age_pct_change, r$age_lm_log2fc,
                r$age_stratum_log2fc, r$age_strata_old_higher,
                ifelse(r$age_three_agree, "YES", "NO")))
  }

  cat("\n========== REGION main effect (sig only, p<0.05) ==========\n\n")
  region_sig <- results %>% filter(region_significant)
  for (i in seq_len(nrow(region_sig))) {
    r <- region_sig[i, ]
    cat(sprintf("%s   region_pval = %.3g\n", r$cell_type, r$region_pval))
    cat(sprintf("  EMM log2FC = %+.3f  (%+.1f%%)   lm = %+.3f   stratum-avg = %+.3f   Dorsal higher in %d/4 strata   agree=%s\n\n",
                r$region_emm_log2fc, r$region_pct_change, r$region_lm_log2fc,
                r$region_stratum_log2fc, r$region_strata_dorsal_higher,
                ifelse(r$region_three_agree, "YES", "NO")))
  }

  invisible(results)
}

if (sys.nframe() == 0L) main()
