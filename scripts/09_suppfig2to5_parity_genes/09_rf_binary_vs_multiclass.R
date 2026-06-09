# 09_rf_binary_vs_multiclass.R — bar plot of random-forest classification
# performance (out-of-bag accuracy): binary (parous vs nulliparous) vs
# multi-class (nulliparous / primiparous / biparous). Read from the SAVED
# model objects (not a refit), so it matches the published models.
# Writes outputs/09_suppfig2to5_parity_genes/supp_fig5_rf_binary_vs_multiclass.png

suppressPackageStartupMessages({ library(tidyverse); library(ggplot2) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

load_rf <- function(path) {
  e <- new.env(); load(path, envir = e)
  for (n in ls(e)) if (inherits(e[[n]], "randomForest")) return(e[[n]])
  stop("no randomForest object in ", path)
}

rf_acc <- function(rf, label) {
  conf <- rf$confusion[, rownames(rf$confusion), drop = FALSE]
  oob  <- if (!is.null(rf$err.rate)) rf$err.rate[nrow(rf$err.rate), "OOB"]
          else 1 - sum(diag(conf)) / sum(conf)
  data.frame(model = label, accuracy = 1 - oob)
}

main <- function() {
  rf_dir <- file.path(REPO_ROOT, "checkpoints", "rf_inputs")
  bin   <- load_rf(file.path(rf_dir, "binary_random_forest.RData"))
  multi <- load_rf(file.path(rf_dir, "multiclass_random_forest.RData"))

  df <- bind_rows(rf_acc(bin,   "Binary\n(2-class)"),
                  rf_acc(multi, "Multi-class\n(3-class)"))
  df$model <- factor(df$model, levels = c("Binary\n(2-class)", "Multi-class\n(3-class)"))
  message(sprintf("[rf-perf] %s",
                  paste(sprintf("%s = %.1f%%", gsub("\n", " ", df$model), 100 * df$accuracy),
                        collapse = " | ")))

  p <- ggplot(df, aes(model, accuracy, fill = model)) +
    geom_col(width = 0.55, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", 100 * accuracy)),
              vjust = -0.4, size = 5, fontface = "bold") +
    scale_fill_manual(values = c("Binary\n(2-class)" = "#2171b5",
                                 "Multi-class\n(3-class)" = "#6baed6"), guide = "none") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1), expand = expansion(mult = c(0, 0.08))) +
    labs(title = "Random-forest classification performance",
         subtitle = "Out-of-bag accuracy",
         x = NULL, y = "Classification accuracy") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(size = 14, face = "bold", colour = "black"),
          plot.subtitle = element_text(size = 10, colour = "gray35"),
          axis.text.x = element_text(size = 12, colour = "black"),
          panel.grid.major.x = element_blank())

  ggsave(file.path(OUT$suppfig234, "supp_fig5_rf_binary_vs_multiclass.png"),
         p, width = 6, height = 5.5, dpi = 600)
  message("[rf-perf] wrote supp_fig5_rf_binary_vs_multiclass.png")
}

if (sys.nframe() == 0L) main()
