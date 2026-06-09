# 07_binary_deseq_volcano.R — Supp Fig 3: binary parity DESeq2 volcano
# (parous vs nulliparous) with age covariate.

suppressPackageStartupMessages({ library(tidyverse); library(ggrepel) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

main <- function() {
  out_dir <- OUT$suppfig234
  full_path <- file.path(out_dir, "DESeq2", "all_genes_parity_deseq2.csv")
  fallback  <- file.path(out_dir, "DESeq2", "significant_genes_parity.csv")
  df <- read.csv(if (file.exists(full_path)) full_path else fallback,
                 stringsAsFactors = FALSE)

  # Annotate each row by combined p-value + fold-change criteria
  df <- df %>% mutate(
    log10p = -log10(pmax(pvalue, 1e-300)),
    label  = ifelse(is.na(external_gene_name) | external_gene_name == "",
                    gene_id, external_gene_name),
    bucket = case_when(
      padj < 0.05                                 ~ "FDR < 0.05",
      padj < 0.20                                 ~ "FDR < 0.20",
      pvalue < 0.05 & abs(log2FoldChange) > 0.3   ~ "p < 0.05 & |FC| > 0.3",
      pvalue < 0.05                               ~ "p < 0.05",
      abs(log2FoldChange) > 0.3                   ~ "|FC| > 0.3",
      TRUE                                        ~ "n.s."
    ),
    bucket = factor(bucket,
                    levels = c("FDR < 0.05", "FDR < 0.20",
                               "p < 0.05 & |FC| > 0.3",
                               "p < 0.05", "|FC| > 0.3", "n.s."))
  )

  top_labels <- df %>% arrange(padj) %>% slice_head(n = 25)

  cols <- c("FDR < 0.05" = "#D7263D",
            "FDR < 0.20" = "#9B59B6",
            "p < 0.05 & |FC| > 0.3" = "#F18F01",
            "p < 0.05" = "#2E86AB",
            "|FC| > 0.3" = "#06D6A0",
            "n.s." = "grey75")

  p <- ggplot(df, aes(x = log2FoldChange, y = log10p, colour = bucket)) +
    geom_point(alpha = 0.55, size = 1.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#2E86AB") +
    geom_hline(yintercept = -log10(min(df$pvalue[df$padj < 0.05], na.rm = TRUE)),
               linetype = "solid", colour = "#D7263D") +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", colour = "gray50") +
    geom_text_repel(data = top_labels, aes(label = label),
                    size = 3, max.overlaps = 25,
                    box.padding = 0.4, point.padding = 0.3) +
    scale_colour_manual(values = cols, name = NULL) +
    labs(title = "Supp Fig 3 | DESeq2 binary parity (parous vs nulliparous)",
         subtitle = "Design: ~ age + parity_binary",
         x = "log2 fold change (Parous / Nulliparous)",
         y = "-log10(p-value)") +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom")

  out_png <- file.path(out_dir, "supp_fig3_binary_deseq_volcano.png")
  out_pdf <- file.path(out_dir, "supp_fig3_binary_deseq_volcano.pdf")
  ggsave(out_pdf, p, width = 9, height = 7)
  ggsave(out_png, p, width = 9, height = 7, dpi = 600)
  message(sprintf("[supp5] %s", out_png))
}

if (sys.nframe() == 0L) main()
