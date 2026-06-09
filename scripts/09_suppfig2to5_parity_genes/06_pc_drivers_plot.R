# 06_pc_drivers_plot.R — Supp Fig 4: top driver genes for PC5/PC6/PC8 (parity-associated PCs).
# Reads PC_Drivers/*_Parity_Binary_improved_drivers.csv from 03_pc_driver_analysis.R.

suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

main <- function() {
  out_dir <- OUT$suppfig234
  pc_dir  <- file.path(out_dir, "PC_Drivers")

  read_pc <- function(pc_label) {
    f <- file.path(pc_dir, sprintf("%s_Parity_Binary_improved_drivers.csv", pc_label))
    df <- read.csv(f, stringsAsFactors = FALSE)
    df$PC <- pc_label
    df
  }

  pcs <- bind_rows(read_pc("PC5"), read_pc("PC6"), read_pc("PC8"))
  pcs$label <- ifelse(is.na(pcs$external_gene_name) | pcs$external_gene_name == "",
                      pcs$gene_name, pcs$external_gene_name)

  one_pc <- function(pc_label, color) {
    d <- pcs %>% filter(PC == pc_label) %>%
      arrange(desc(contribution_score)) %>% slice_head(n = 20)
    ggplot(d, aes(x = reorder(label, contribution_score),
                  y = contribution_score)) +
      geom_col(fill = color, alpha = 0.85) +
      coord_flip() +
      labs(title = sprintf("%s — top 20 parity drivers", pc_label),
           x = NULL, y = "Contribution score") +
      theme_minimal(base_size = 10) +
      theme(axis.text.y = element_text(size = 8))
  }

  p <- one_pc("PC5", "#EF476F") +
       one_pc("PC6", "#06D6A0") +
       one_pc("PC8", "#118AB2") +
       plot_annotation(title = "Supp Fig 4 | Parity-associated PC driver genes",
                       subtitle = sprintf("PC5 n=%d  |  PC6 n=%d  |  PC8 n=%d",
                                          sum(pcs$PC == "PC5"),
                                          sum(pcs$PC == "PC6"),
                                          sum(pcs$PC == "PC8")))

  out_png <- file.path(out_dir, "supp_fig4_pc_drivers.png")
  out_pdf <- file.path(out_dir, "supp_fig4_pc_drivers.pdf")
  ggsave(out_pdf, p, width = 16, height = 7)
  ggsave(out_png, p, width = 16, height = 7, dpi = 600)
  message(sprintf("[supp3] %s", out_png))
}

if (sys.nframe() == 0L) main()
