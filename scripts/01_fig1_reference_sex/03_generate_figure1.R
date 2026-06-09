# 03_generate_figure1.R — renders the 6-panel Figure 1 PDF.
# Condensed from Set 1/Transfer/Final/generate_nature_biotech_figures_final_v4.R (lines 95-524).

suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(patchwork)
})
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

SEX_COLORS <- c(male = "#377EB8", female = "#E41A1C", mixed = "#4DAF4A")

nature_theme <- theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "sans"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9, color = "black"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 9, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    plot.margin = margin(10, 20, 10, 10, "pt")
  )

panel <- function(m, y, ylab, fmt = "%.2f") {
  ggplot(m, aes(x = sex, y = .data[[y]], fill = sex)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.5) +
    geom_text(aes(label = sprintf(fmt, .data[[y]])), vjust = -0.3, size = 3) +
    facet_wrap(~ dataset_label, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = SEX_COLORS, name = "Reference") +
    labs(x = NULL, y = ylab) +
    nature_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "none")
}

main <- function() {
  out_dir <- OUT$fig1
  metrics_df <- read.csv(file.path(out_dir, "normalized_metrics_table.csv"))
  metrics_df$sex <- factor(metrics_df$sex, levels = c("male","female","mixed"))

  a <- panel(metrics_df, "types_per_k_cells",         "Cell types detected\nper 1,000 ref cells", "%.2f")
  b <- panel(metrics_df, "entropy_per_k_cells",       "Entropy per\n1,000 ref cells",              "%.3f")
  c <- panel(metrics_df, "detected_per_sample_per_k", "Types per sample\nper 1,000 ref cells",     "%.3f")
  d <- panel(metrics_df, "diversity_index",           "Normalized\ndiversity index",               "%.3f")
  e <- panel(metrics_df, "mean_prop_per_k_cells",     "Mean proportion\nper 1,000 ref cells",      "%.4f")
  f <- panel(metrics_df, "sparsity_per_k_cells",      "% zero proportions\nper 1,000 ref cells",   "%.1f")

  fig1 <- (a | b | c) / (d | e | f) +
    plot_annotation(title = "Figure 1 | Sex-matched reference datasets optimize deconvolution",
                    tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 14))

  out_pdf <- file.path(out_dir, "Figure1_Sex_Reference_Performance.pdf")
  out_png <- file.path(out_dir, "Figure1_Sex_Reference_Performance.png")
  ggsave(out_pdf, fig1, width = 12, height = 8, device = cairo_pdf)
  ggsave(out_png, fig1, width = 12, height = 8, dpi = 300)
  message(sprintf("[fig1] Wrote %s", out_pdf))
}

if (sys.nframe() == 0L) main()
