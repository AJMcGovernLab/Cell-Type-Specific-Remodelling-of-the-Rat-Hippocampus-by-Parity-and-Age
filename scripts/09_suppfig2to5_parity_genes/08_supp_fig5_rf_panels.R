# 08_supp_fig5_rf_panels.R — extra Supp Fig 5 random-forest panels:
#   * cross-validation accuracy vs gene-set size (peak at 23 genes)
#   * classification performance (confusion matrix from the saved multi-class model)
#
# Writes PNGs to outputs/09_suppfig2to5_parity_genes/.

suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(scales)
})
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

# ---------------------------------------------------------------------------
# Supp 5 — RF classification performance + cross-validation accuracy
# ---------------------------------------------------------------------------
rf_panels <- function() {
  rf_dir <- file.path(OUT$suppfig234, "Random_Forest")

  ## --- cross-validation accuracy vs gene-set size (peak at 23) ----------
  ## From the saved CV results (5-fold CV x 3 repeats over top-k genes by
  ## random-forest importance). Re-plotted from saved data, not refit.
  cv <- read.csv(file.path(rf_dir, "cv_accuracy_by_geneset_20_40.csv"),
                 stringsAsFactors = FALSE)
  peak <- cv$n_genes[which.max(cv$mean_accuracy)]
  p_cv <- ggplot(cv, aes(n_genes, mean_accuracy)) +
    geom_ribbon(aes(ymin = mean_accuracy - se_accuracy, ymax = mean_accuracy + se_accuracy),
                fill = "gray85", alpha = 0.6) +
    geom_line(colour = "gray40") +
    geom_point(aes(colour = n_genes == peak), size = 2.3) +
    geom_vline(xintercept = peak, linetype = "dashed", colour = "#d73027") +
    annotate("text", x = peak + 0.6, y = max(cv$mean_accuracy),
             label = sprintf("peak = %d genes\n(%.1f%% accuracy)", peak, 100 * max(cv$mean_accuracy)),
             hjust = 0, vjust = 1, size = 3.3, colour = "#d73027") +
    scale_colour_manual(values = c(`TRUE` = "#d73027", `FALSE` = "#2171b5"), guide = "none") +
    scale_x_continuous(breaks = seq(20, 40, 2)) +
    labs(title = "Cross-validation accuracy by gene-set size",
         subtitle = "5-fold cross-validation (3 repeats); top genes by random-forest importance",
         x = "Number of genes", y = "Mean CV classification accuracy") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(size = 13, face = "bold", colour = "black"),
          plot.subtitle = element_text(size = 9.5, colour = "gray35"),
          panel.grid.minor = element_blank())
  ggsave(file.path(OUT$suppfig234, "supp_fig5_rf_cv_accuracy.png"),
         p_cv, width = 7, height = 5, dpi = 600)
  message(sprintf("[supp5] wrote supp_fig5_rf_cv_accuracy.png (peak at %d genes)", peak))

  ## --- classification performance (confusion matrix from saved model) ---
  rdata <- file.path(REPO_ROOT, "checkpoints", "rf_inputs", "multiclass_random_forest.RData")
  if (!file.exists(rdata)) { message("[supp5] RF RData not found; skipping performance panel"); return(invisible()) }
  e <- new.env(); load(rdata, envir = e)
  rf <- NULL
  for (n in ls(e)) if (inherits(e[[n]], "randomForest")) { rf <- e[[n]]; break }
  if (is.null(rf)) { message("[supp5] no randomForest object in RData; skipping"); return(invisible()) }

  cm <- as.data.frame(rf$confusion)
  classes <- rownames(rf$confusion)
  conf <- rf$confusion[, classes, drop = FALSE]
  oob  <- if (!is.null(rf$err.rate)) rf$err.rate[nrow(rf$err.rate), "OOB"] else NA

  long <- as.data.frame(as.table(as.matrix(conf)))
  names(long) <- c("Actual", "Predicted", "n")
  long$Actual    <- factor(long$Actual, levels = rev(classes))
  long$Predicted <- factor(long$Predicted, levels = classes)

  p_conf <- ggplot(long, aes(Predicted, Actual, fill = n)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = n), size = 5, fontface = "bold") +
    scale_fill_gradient(low = "#f7fbff", high = "#2171b5", name = "Count") +
    labs(title = "Random-forest classification performance",
         subtitle = sprintf("Multi-class parity (out-of-bag error = %.1f%%)", 100 * oob),
         x = "Predicted class", y = "Actual class") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(size = 13, face = "bold", colour = "black"),
          plot.subtitle = element_text(size = 10, colour = "gray35"),
          axis.text = element_text(size = 11, colour = "black"),
          panel.grid = element_blank())
  ggsave(file.path(OUT$suppfig234, "supp_fig5_rf_classification_performance.png"),
         p_conf, width = 6, height = 5, dpi = 600)
  message("[supp5] wrote supp_fig5_rf_classification_performance.png")
}

rf_panels()
message("[08] supplementary RF panels complete")
