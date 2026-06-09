# 06_regenerate_panel_a.R
# ========================
# Regenerates Figure 2 Panel A (cell-type reduction through harmonization)
# using the per-dataset-sum framing.
#
# What's different from the original Manuscript_Figures rendering:
#   * Uses cell_type_mapping_table.csv (canonical first-pass harmonization).
#   * Adds a fourth bar group, "Total (sum across datasets)", that makes
#     the 349 -> 27 reduction visible alongside the per-dataset bars
#     (169 / 125 / 55  ->  16 / 10 / 1). This avoids the apparent mismatch
#     between bar sums and the previously-cited 201 / 23 deduplicated totals.
#   * Saves both PDF and PNG into Repository/outputs/02_fig2_harmonization/.
#
# Run:
#   Rscript scripts/02_fig2_harmonization/06_regenerate_panel_a.R
#
# Outputs:
#   outputs/02_fig2_harmonization/Figure2_PanelA_per_dataset_sum.pdf
#   outputs/02_fig2_harmonization/Figure2_PanelA_per_dataset_sum.png
#   outputs/02_fig2_harmonization/figure2_panel_a_data.csv     (numbers behind the bars)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

main <- function() {
  out_dir <- OUT$fig2
  mapping_path <- file.path(out_dir, "cell_type_mapping_table.csv")
  if (!file.exists(mapping_path)) {
    # Fallback to the ReproductiveExperienceAndAgeDeconvolution snapshot
    fallback <- "f:/Parity/2/ReproductiveExperienceAndAgeDeconvolution/2_Cell_Harmonization/cell_type_mapping_table.csv"
    if (!file.exists(fallback)) {
      stop(sprintf("cell_type_mapping_table.csv not found at %s or %s",
                   mapping_path, fallback))
    }
    mapping_path <- fallback
  }
  mapping <- read.csv(mapping_path, stringsAsFactors = FALSE)

  # Per-dataset originals + harmonized counts
  per_dataset <- mapping %>%
    group_by(dataset) %>%
    summarise(
      Originals  = n_distinct(original_name),
      Harmonized = n_distinct(unified_name),
      .groups    = "drop"
    ) %>%
    mutate(
      Dataset_Clean = recode(
        dataset,
        "mouse10x_2020"        = "Mouse 10X 2020",
        "mouse_smartseq_2019"  = "Mouse Smart-seq 2019",
        "yao_hippo_10x"        = "Yao Hippocampus 10X"
      )
    )

  # Append the "Total (sum across datasets)" group:
  # Originals  = sum of per-dataset distinct originals  (= 349)
  # Harmonized = sum of per-dataset distinct harmonized (= 27)
  total_row <- data.frame(
    dataset       = "ALL",
    Originals     = sum(per_dataset$Originals),
    Harmonized    = sum(per_dataset$Harmonized),
    Dataset_Clean = "Total (sum across datasets)"
  )
  combined <- bind_rows(per_dataset, total_row) %>%
    mutate(Dataset_Clean = factor(
      Dataset_Clean,
      levels = c("Mouse 10X 2020", "Mouse Smart-seq 2019",
                 "Yao Hippocampus 10X", "Total (sum across datasets)")
    ))

  # Save the underlying numbers next to the figure
  write.csv(combined, file.path(out_dir, "figure2_panel_a_data.csv"),
            row.names = FALSE)

  # Reshape for plotting
  long <- combined %>%
    pivot_longer(cols = c(Originals, Harmonized),
                 names_to = "Stage", values_to = "Count") %>%
    mutate(Stage = factor(
      ifelse(Stage == "Originals", "Before harmonization", "After harmonization"),
      levels = c("Before harmonization", "After harmonization")
    ))

  panel_a <- ggplot(long,
                    aes(x = Dataset_Clean, y = Count, fill = Stage)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = Count),
              position = position_dodge(width = 0.8),
              vjust = -0.3, size = 4) +
    geom_vline(xintercept = 3.5, linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    annotate("text", x = 3.5, y = max(long$Count) * 1.05,
             label = "  per-dataset    |    sum",
             hjust = 0.5, size = 3.2, colour = "grey40") +
    scale_fill_manual(values = c("Before harmonization" = "#FC8D62",
                                  "After harmonization" = "#66C2A5")) +
    labs(
      title    = "Cell-type reduction through harmonization",
      subtitle = "Per-dataset bars sum to the rightmost group (349 originals -> 27 harmonized).",
      x        = "Reference dataset",
      y        = "Number of cell-type annotations",
      fill     = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(fill = NA, colour = "black", linewidth = 0.4),
      axis.text.x      = element_text(angle = 25, hjust = 1),
      legend.position  = "top",
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 11, colour = "grey30")
    )

  ggsave(file.path(out_dir, "Figure2_PanelA_per_dataset_sum.pdf"),
         panel_a, width = 8, height = 5)
  ggsave(file.path(out_dir, "Figure2_PanelA_per_dataset_sum.png"),
         panel_a, width = 8, height = 5, dpi = 300)

  message(sprintf(
    "[panel-a] %d originals (sum) -> %d harmonized (sum) across 3 references; per-dataset: %s",
    sum(per_dataset$Originals),
    sum(per_dataset$Harmonized),
    paste(sprintf("%s %d->%d",
                  per_dataset$Dataset_Clean,
                  per_dataset$Originals,
                  per_dataset$Harmonized),
          collapse = "; ")
  ))
}

if (sys.nframe() == 0L) main()
