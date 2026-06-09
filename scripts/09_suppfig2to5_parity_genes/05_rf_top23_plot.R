# 05_rf_top23_plot.R — Supp Fig 5: top-23 RF parity-predictive gene bar plot.
# Reads Random_Forest/top23_genes_with_names.csv produced by 02_random_forest.R.

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

main <- function() {
  out_dir <- OUT$suppfig234
  rf_file <- file.path(out_dir, "Random_Forest", "top23_genes_with_names.csv")
  df <- read.csv(rf_file, stringsAsFactors = FALSE)

  df$display_name <- ifelse(
    is.na(df$external_gene_name) | df$external_gene_name == "",
    df$gene, df$external_gene_name)

  df <- df %>% arrange(desc(mean_decrease_accuracy))

  p <- ggplot(df, aes(x = reorder(display_name, mean_decrease_accuracy),
                      y = mean_decrease_accuracy)) +
    geom_col(fill = "#118AB2", alpha = 0.85) +
    geom_text(aes(label = sprintf("%.2f", mean_decrease_accuracy)),
              hjust = -0.1, size = 3) +
    coord_flip() +
    labs(title  = "Supp Fig 5 | Top-23 parity-predictive genes (RF, multi-class)",
         subtitle = "Mean decrease in accuracy",
         x = NULL, y = "Mean Decrease in Accuracy") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 9))

  out_png <- file.path(out_dir, "supp_fig5_rf_top23.png")
  out_pdf <- file.path(out_dir, "supp_fig5_rf_top23.pdf")
  ggsave(out_pdf, p, width = 8, height = 7)
  ggsave(out_png, p, width = 8, height = 7, dpi = 600)
  message(sprintf("[supp2] %s", out_png))
}

if (sys.nframe() == 0L) main()
