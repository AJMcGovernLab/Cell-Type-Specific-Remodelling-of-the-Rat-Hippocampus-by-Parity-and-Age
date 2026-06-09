# 05_marginal_parity_log2fc.R — properly-conditioned parity fold change.
#
# The simple parous-vs-nullip contrast (in 02_compute_log2fc.R) pools across
# age and region, so the dominant variance from age and region masks the
# parity signal. The ANOVA p-value comes from a model that *adjusts* for
# age and region, so the right matching directional estimate is the
# parity main effect AFTER conditioning on age and region.
#
# Computed three ways for robustness:
#   (a) stratum-mean log2FC: within each (age × region) cell, compute
#       log2(mean_parous / mean_nullip), then average the 4 strata-level
#       log2FCs (this is the unweighted estimated marginal log2FC).
#   (b) lm-coefficient on the proportion scale: fit
#       proportion ~ age + region + parity_binary, take the parity_binary
#       coefficient, divide by intercept to get a relative effect, then log2.
#   (c) EMM ratio via base-R fitted means: predict at parous and nullip
#       holding age and region balanced, ratio them, log2.
#
# Outputs: outputs/06_fig5ab_parity_proportions/parity_marginal_log2fc.csv

suppressPackageStartupMessages({ library(tidyverse); library(broom) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

PSEUDOCOUNT <- 1e-6
safe_log2fc <- function(num, den) log2((num + PSEUDOCOUNT) / (den + PSEUDOCOUNT))

per_cell_marginal <- function(ct, deconv, dataset_name, metadata) {
  d <- deconv %>%
    select(sample, proportion = all_of(ct)) %>%
    left_join(metadata, by = "sample") %>%
    mutate(parity_binary = factor(ifelse(parity == "Nulliparous", "Nulliparous", "Parous"),
                                  levels = c("Nulliparous", "Parous")),
           age    = factor(age,    levels = c("Young","Old")),
           region = factor(region, levels = c("Dorsal","Ventral")))

  # (a) stratum log2FC, then average
  stratum <- d %>%
    group_by(age, region, parity_binary) %>%
    summarise(mu = mean(proportion), .groups = "drop") %>%
    pivot_wider(names_from = parity_binary, values_from = mu) %>%
    mutate(log2fc = safe_log2fc(Parous, Nulliparous))
  stratum_log2fc <- mean(stratum$log2fc)
  n_strata_parous_higher <- sum(stratum$Parous > stratum$Nulliparous)

  # (b) lm-coefficient (additive model)
  fit <- lm(proportion ~ age + region + parity_binary, data = d)
  cf  <- coef(fit)
  intercept <- cf["(Intercept)"]
  parity_eff <- cf["parity_binaryParous"]
  # baseline (Nullip, Young, Dorsal) = intercept; parous shift = + parity_eff
  lm_log2fc <- safe_log2fc(intercept + parity_eff, intercept)

  # (c) EMM via balanced predict
  grid <- expand.grid(age = c("Young","Old"), region = c("Dorsal","Ventral"),
                      parity_binary = c("Nulliparous","Parous"))
  grid$pred <- predict(fit, newdata = grid)
  emm <- grid %>% group_by(parity_binary) %>% summarise(mu = mean(pred), .groups = "drop")
  emm_nullip <- emm$mu[emm$parity_binary == "Nulliparous"]
  emm_parous <- emm$mu[emm$parity_binary == "Parous"]
  emm_log2fc <- safe_log2fc(emm_parous, emm_nullip)
  emm_pct    <- 100 * (emm_parous - emm_nullip) / emm_nullip

  # full model p-value (matches manuscript)
  full_fit <- aov(proportion ~ age * parity_binary * region, data = d)
  ap <- summary(full_fit)[[1]]
  parity_pval <- ap["parity_binary", "Pr(>F)"]

  data.frame(
    cell_type        = ct,
    dataset          = dataset_name,
    n_total          = nrow(d),
    parity_pval_3way = parity_pval,
    EMM_nullip       = emm_nullip,
    EMM_parous       = emm_parous,
    log2fc_emm       = emm_log2fc,
    pct_change_emm   = emm_pct,
    log2fc_lm_additive = lm_log2fc,
    log2fc_stratum_avg = stratum_log2fc,
    n_strata_parous_higher = n_strata_parous_higher,  # out of 4
    direction_consistent = ifelse(
      sign(emm_log2fc) == sign(lm_log2fc) & sign(emm_log2fc) == sign(stratum_log2fc),
      "ALL THREE AGREE",
      "DISAGREE"
    )
  )
}

main <- function() {
  out_dir  <- OUT$fig5ab
  metadata <- read.csv(SAMPLE_META)

  targets <- list(
    list(cell = "356_CA3-do",  rds = "mouse10x_2020_female_scdc.rds",       ds = "mouse10x_2020"),
    list(cell = "358_CA3-do",  rds = "mouse_smartseq_2019_female_scdc.rds", ds = "mouse_smartseq_2019"),
    list(cell = "376_Astro",   rds = "mouse10x_2020_female_scdc.rds",       ds = "mouse10x_2020"),
    list(cell = "78_Sst HPF",  rds = "mouse10x_2020_female_scdc.rds",       ds = "mouse10x_2020")
  )

  results <- bind_rows(lapply(targets, function(t) {
    res <- readRDS(file.path(SCDC_CHECKPOINTS, t$rds))
    deconv <- if (!is.null(res$proportions)) res$proportions
              else { df <- as.data.frame(res$raw_result$prop.est.mvw); df$sample <- rownames(df); df }
    per_cell_marginal(t$cell, deconv, t$ds, metadata)
  }))

  write.csv(results, file.path(out_dir, "parity_marginal_log2fc.csv"), row.names = FALSE)

  cat("\n========== MARGINAL PARITY LOG2FC (adjusted for age + region) ==========\n\n")
  for (i in seq_len(nrow(results))) {
    r <- results[i, ]
    cat(sprintf("%s  (%s)   parity_pval = %.4f\n", r$cell_type, r$dataset, r$parity_pval_3way))
    cat(sprintf("  EMM      nullip=%.5f  parous=%.5f   log2FC = %+.4f   %+.2f%%\n",
                r$EMM_nullip, r$EMM_parous, r$log2fc_emm, r$pct_change_emm))
    cat(sprintf("  lm beta  (additive model)               log2FC = %+.4f\n", r$log2fc_lm_additive))
    cat(sprintf("  stratum mean of within-stratum log2FCs  log2FC = %+.4f\n", r$log2fc_stratum_avg))
    cat(sprintf("  parous higher in %d of 4 age×region strata\n", r$n_strata_parous_higher))
    cat(sprintf("  three estimates: %s\n\n", r$direction_consistent))
  }

  invisible(results)
}

if (sys.nframe() == 0L) main()
