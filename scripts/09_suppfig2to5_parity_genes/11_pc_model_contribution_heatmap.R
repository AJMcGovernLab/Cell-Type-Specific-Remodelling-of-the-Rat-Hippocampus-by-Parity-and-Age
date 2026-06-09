# 11_pc_model_contribution_heatmap.R — Supp Fig 4 panel.
#
# Mean driver-gene contribution for PC1..PC8 x 7 statistical models, showing the
# specificity of PC5/PC6/PC8 to parity. Per gene, the contribution to a PC's
# association with a model is
#   contribution = log10(reduced_p / base_p) + 100 * (base_R2 - reduced_R2)
# where reduced_* is the model refit after removing the gene's loading-weighted
# contribution from the PC scores (same definition as 03_pc_driver_analysis.R,
# generalised to all 7 models and vectorised over genes). Cell value = mean
# contribution of driver genes (contribution > 0.1); cells whose overall model
# fit is significant (p < 0.05) are outlined. Deterministic (prcomp + lm), no
# randomness.

suppressPackageStartupMessages({ library(tidyverse); library(ggplot2); library(data.table) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

DRIVER_THRESH <- 0.1
PCS <- paste0("PC", 1:8)
MODELS <- list(
  "Age"             = "~ age",
  "Region"          = "~ region",
  "Parity"          = "~ parity_binary",
  "Age x Region"    = "~ age + region + age:region",
  "Age x Parity"    = "~ age + parity_binary + age:parity_binary",
  "Region x Parity" = "~ region + parity_binary + region:parity_binary",
  "Three-way"       = "~ age + region + parity_binary + age:region + age:parity_binary + region:parity_binary + age:region:parity_binary"
)

main <- function() {
  ckpt <- file.path(CHECKPOINT_DIR, "pca_drivers")
  expr <- read.csv(file.path(BULK_DIR, "clean_normalized_counts.csv"), check.names = FALSE)
  expr <- as.matrix(expr); rownames(expr) <- paste0("Gene_", seq_len(nrow(expr)))
  loadings <- read.csv(file.path(ckpt, "pca_loadings_top10.csv"), stringsAsFactors = FALSE)
  coords <- read.csv(file.path(ckpt, "binary_parity_pca_coordinates.csv"), stringsAsFactors = FALSE) %>%
    mutate(age = factor(age), region = factor(region), parity_binary = factor(parity_binary))

  expr <- expr[, match(coords$sample, colnames(expr))]   # align samples to coords
  n <- nrow(coords); G <- nrow(expr)

  # overall-fit R2 + F-test p for a response matrix Y (n x G) under design X
  fit_stats <- function(Y, X) {
    qrX <- qr(X); p <- ncol(X); dfm <- p - 1; dfr <- n - p
    res <- qr.resid(qrX, Y)
    sse <- colSums(res^2)
    Yc  <- sweep(Y, 2, colMeans(Y)); sst <- colSums(Yc^2)
    r2  <- 1 - sse / sst
    Fst <- (r2 / dfm) / ((1 - r2) / dfr)
    list(r2 = r2, p = pf(Fst, dfm, dfr, lower.tail = FALSE))
  }

  res_rows <- list()
  for (pc in PCS) {
    sc  <- coords[[pc]]                                   # n-vector PC scores
    ld  <- loadings[[pc]][match(rownames(expr), loadings$gene)]
    contrib <- t(expr * ld)                               # n x G : gene g column = expr_g * loading_g
    Yred <- sc - contrib                                  # n x G reduced PC scores (sc recycled per column)
    for (mn in names(MODELS)) {
      X <- model.matrix(as.formula(sub("~", paste(pc, "~"), MODELS[[mn]])), data = coords)
      base <- fit_stats(matrix(sc, ncol = 1), X)
      red  <- fit_stats(Yred, X)
      contribution <- log10(pmax(red$p, 1e-300) / max(base$p, 1e-300)) + 100 * (base$r2 - red$r2)
      drv <- contribution[contribution > DRIVER_THRESH]
      res_rows[[length(res_rows) + 1]] <- data.frame(
        PC = pc, Model = mn,
        mean_contrib = if (length(drv)) mean(drv) else 0,
        n_drivers = length(drv),
        model_p = base$p, sig = base$p < 0.05)
    }
  }
  d <- bind_rows(res_rows) %>%
    mutate(PC = factor(PC, levels = rev(PCS)), Model = factor(Model, levels = names(MODELS)))

  sig_d <- subset(d, sig)
  p_ht <- ggplot(d, aes(Model, PC)) +
    geom_tile(fill = "white", colour = "grey88", linewidth = 0.5) +
    geom_tile(data = sig_d, aes(fill = mean_contrib), colour = "black", linewidth = 0.9) +
    geom_text(data = sig_d, aes(label = sprintf("%.1f", mean_contrib)),
              size = 3.9, fontface = "bold",
              colour = ifelse(sig_d$mean_contrib > max(sig_d$mean_contrib) * 0.55, "white", "black")) +
    scale_fill_viridis_c(option = "C", direction = -1, name = "Mean driver\ncontribution") +
    labs(title = "Mean driver-gene contribution by PC and statistical model",
         subtitle = "Significant associations only (model p < 0.05)",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(size = 13, face = "bold", colour = "black"),
          plot.subtitle = element_text(size = 9, colour = "gray35"),
          axis.text.x = element_text(size = 11, face = "bold", colour = "black", angle = 30, hjust = 1),
          axis.text.y = element_text(size = 12, face = "bold", colour = "black"),
          panel.grid = element_blank(), legend.position = "right")

  ggsave(file.path(OUT$suppfig234, "supp_fig4_pc_model_contribution.png"), p_ht, width = 8.6, height = 5.4, dpi = 600)
  write.csv(d %>% select(PC, Model, mean_contrib, n_drivers, model_p),
            file.path(OUT$suppfig234, "pc_model_driver_contribution.csv"), row.names = FALSE)
  message("[supp4] wrote supp_fig4_pc_model_contribution.png")
  print(d %>% select(PC, Model, mean_contrib, n_drivers, model_p) %>% arrange(PC, Model), digits = 3)
}

if (sys.nframe() == 0L) main()
