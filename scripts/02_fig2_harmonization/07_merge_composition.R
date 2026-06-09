# 07_merge_composition.R — Figure 2c: harmonized categories vs the original
# cell-type annotations merged into each, stacked by source reference dataset.
#
# Recreates the "harmonized cell type against the original names merged" panel:
# each bar is a harmonized category, length = number of original annotations it
# absorbed, split by which reference dataset they came from. Shows the varying
# cluster sizes (largest = DG/CA1/CA2/CA3_dorsal, then Sst_IN, ...) and the
# cross-dataset merging described in the Figure 2 legend.
#
# Writes outputs/02_fig2_harmonization/Figure2_PanelC_merge_composition.{png,pdf}

suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2)
})
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

MIN_CONSTITUENTS <- 2   # show only genuinely merged categories (>=2 originals)

DATASET_LABELS <- c(
  "mouse10x_2020"       = "Mouse 10x (2020)",
  "mouse_smartseq_2019" = "Mouse Smart-seq (2019)",
  "yao_hippo_10x"       = "Yao hippocampus 10x"
)
DATASET_COLOURS <- c(
  "Mouse 10x (2020)"       = "#1b9e77",
  "Mouse Smart-seq (2019)" = "#d95f02",
  "Yao hippocampus 10x"    = "#7570b3"
)

main <- function() {
  out_dir <- OUT$fig2
  map <- read.csv(file.path(out_dir, "refined_cell_type_mapping.csv"),
                  stringsAsFactors = FALSE)

  # constituents merged per (harmonized category, dataset)
  comp <- map %>%
    distinct(dataset, original_name, unified_name) %>%
    count(unified_name, dataset, name = "n")

  totals <- comp %>%
    group_by(unified_name) %>%
    summarise(total = sum(n), .groups = "drop") %>%
    filter(total >= MIN_CONSTITUENTS) %>%
    arrange(total)

  comp <- comp %>%
    filter(unified_name %in% totals$unified_name) %>%
    mutate(unified_name = factor(unified_name, levels = totals$unified_name),
           dataset_lab  = factor(DATASET_LABELS[dataset],
                                 levels = names(DATASET_COLOURS)))

  message(sprintf("[fig2c] %d merged harmonized categories (>=%d originals); %d shown",
                  nrow(totals), MIN_CONSTITUENTS, nrow(totals)))

  p <- ggplot(comp, aes(x = n, y = unified_name, fill = dataset_lab)) +
    geom_col(width = 0.78, colour = "white", linewidth = 0.25) +
    geom_text(data = totals,
              aes(x = total, y = factor(unified_name, levels = totals$unified_name),
                  label = total),
              inherit.aes = FALSE, hjust = -0.35, size = 3.1, colour = "black") +
    scale_fill_manual(values = DATASET_COLOURS, name = "Reference dataset") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(title = "Cell-type harmonization: original annotations merged per category",
         subtitle = sprintf("%d harmonized categories that merged ≥2 of the 201 original annotations across 3 references",
                            nrow(totals)),
         x = "Number of original cell-type annotations merged",
         y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title       = element_text(size = 13, face = "bold", colour = "black"),
          plot.subtitle    = element_text(size = 10, colour = "gray35"),
          axis.title.x     = element_text(size = 11, face = "bold", colour = "black"),
          axis.text.y      = element_text(size = 10, colour = "black"),
          axis.text.x      = element_text(size = 10, colour = "black"),
          legend.position  = "bottom",
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank())

  ggsave(file.path(out_dir, "Figure2_PanelC_merge_composition.png"),
         p, width = 8.5, height = 7.5, dpi = 600)
  ggsave(file.path(out_dir, "Figure2_PanelC_merge_composition.pdf"),
         p, width = 8.5, height = 7.5)
  message("[fig2c] wrote Figure2_PanelC_merge_composition.{png,pdf}")
}

if (sys.nframe() == 0L) main()
