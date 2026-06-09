#!/usr/bin/env Rscript

# ============================================================================
# Cell Type Specificity Analysis - Weighted Meta-Analysis Results
# Identify pathways specific to each cell type using weighted results
# ============================================================================

library(tidyverse)
library(ggplot2)

cat("================================================================\n")
cat("WEIGHTED CELL TYPE SPECIFICITY ANALYSIS\n")
cat("================================================================\n\n")

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------

# Paths via curated config + v7-adapter shim (02b_v7_adapter.R).
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))
input_dir  <- OUT$fig5ce
output_dir <- file.path(OUT$fig5ce, "cell_type_specificity")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
significance_threshold <- 0.05
nes_difference_threshold <- 1.0  # Minimum NES difference for preferential enrichment
specificity_threshold <- 0.8     # Minimum specificity score

cat("Input directory:", input_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("NES difference threshold:", nes_difference_threshold, "\n")
cat("Specificity threshold:", specificity_threshold, "\n\n")

# ------------------------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------------------------

cat("Loading weighted meta-analysis results...\n")

# Load weighted meta-analysis results
meta_file <- file.path(input_dir, "weighted_meta_analysis_significant_v7.csv")
meta_data <- read_csv(meta_file, show_col_types = FALSE)

cat("Data loaded:\n")
cat("- Total pathways:", nrow(meta_data), "\n")
cat("- Cell types found:", paste(unique(meta_data$meta_group), collapse = ", "), "\n")
cat("- Databases:", paste(unique(meta_data$database), collapse = ", "), "\n\n")

# ------------------------------------------------------------------------------
# IDENTIFY CELL TYPE-SPECIFIC PATHWAYS
# ------------------------------------------------------------------------------

cat("Identifying cell type-specific pathways...\n")

# Group pathways by ID to find which appear in multiple cell types
pathway_analysis <- meta_data %>%
  group_by(ID, Description, database) %>%
  summarise(
    n_cell_types = n(),
    cell_types = paste(meta_group, collapse = ";"),
    nes_values = list(weighted_mean_NES),
    fdr_values = list(weighted_meta_padj_corrected),
    meta_groups = list(meta_group),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    # Calculate NES statistics
    min_nes = min(nes_values),
    max_nes = max(nes_values),
    nes_range = max_nes - min_nes,
    mean_nes = mean(nes_values),
    
    # Determine specificity type
    specificity_type = case_when(
      n_cell_types == 1 ~ "Highly Specific",
      n_cell_types > 1 & nes_range > nes_difference_threshold ~ "Preferentially Enriched",
      n_cell_types > 1 & nes_range <= nes_difference_threshold ~ "Shared",
      TRUE ~ "Other"
    ),
    
    # Identify dominant cell type for preferentially enriched pathways
    dominant_cell_type = if (n_cell_types > 1 & nes_range > nes_difference_threshold) {
      meta_groups[[which.max(abs(nes_values))]]
    } else if (n_cell_types == 1) {
      meta_groups[[1]]
    } else {
      "Shared"
    },
    
    # Calculate specificity score
    specificity_score = if (n_cell_types == 1) {
      1.0
    } else {
      max(abs(nes_values)) / mean(abs(nes_values))
    }
  ) %>%
  ungroup()

# Summary of specificity types
specificity_summary <- pathway_analysis %>%
  count(specificity_type, name = "n_pathways") %>%
  mutate(percentage = round(n_pathways / sum(n_pathways) * 100, 1))

cat("Pathway specificity distribution:\n")
print(specificity_summary)
cat("\n")

# ------------------------------------------------------------------------------
# ANALYZE HIGHLY SPECIFIC PATHWAYS
# ------------------------------------------------------------------------------

cat("Analyzing highly specific pathways...\n")

# Get highly specific pathways
specific_pathways <- pathway_analysis %>%
  filter(specificity_type == "Highly Specific") %>%
  mutate(
    cell_type_clean = case_when(
      dominant_cell_type == "CA3_do_combined" ~ "CA3 Dorsal",
      dominant_cell_type == "376_Astro" ~ "Astrocytes",
      dominant_cell_type == "78_Sst_HPF" ~ "Sst+ HPF",
      TRUE ~ dominant_cell_type
    )
  ) %>%
  arrange(cell_type_clean, desc(abs(mean_nes)))

# Summary by cell type
cell_type_summary <- specific_pathways %>%
  group_by(cell_type_clean) %>%
  summarise(
    n_specific_pathways = n(),
    mean_nes = round(mean(mean_nes), 3),
    median_nes = round(median(mean_nes), 3),
    max_nes = round(max(abs(mean_nes)), 3),
    n_upregulated = sum(mean_nes > 0),
    n_downregulated = sum(mean_nes < 0),
    pct_upregulated = round(n_upregulated / n_specific_pathways * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_specific_pathways))

cat("Cell type-specific pathway summary:\n")
print(cell_type_summary)
cat("\n")

# ------------------------------------------------------------------------------
# COMPARE WITH PREVIOUS EQUAL-WEIGHTED ANALYSIS
# ------------------------------------------------------------------------------

cat("Comparing with previous equal-weighted specificity...\n")

# Load previous specificity results if available
previous_specificity_file <- file.path(unix_base_dir, "weighted_parity_analysis/results/enhanced_percentage_weighting/cell_type_specificity/specificity_summary.csv")

if (file.exists(previous_specificity_file)) {
  
  previous_summary <- read_csv(previous_specificity_file, show_col_types = FALSE) %>%
    mutate(analysis_type = "Equal-weighted") %>%
    select(meta_group, n_specific_pathways, analysis_type) %>%
    rename(cell_type = meta_group)
  
  # Current summary for comparison
  current_summary <- cell_type_summary %>%
    mutate(analysis_type = "Weighted (5%=50%, 10%=40%, 25%=10%)") %>%
    select(cell_type_clean, n_specific_pathways, analysis_type) %>%
    rename(cell_type = cell_type_clean)
  
  # Combine for comparison
  comparison_data <- bind_rows(previous_summary, current_summary) %>%
    mutate(
      cell_type = case_when(
        cell_type == "CA3_do_combined" ~ "CA3 Dorsal",
        cell_type == "376_Astro" ~ "Astrocytes",
        cell_type == "78_Sst_HPF" ~ "Sst+ HPF",
        TRUE ~ cell_type
      )
    )
  
  # Create comparison table
  comparison_table <- comparison_data %>%
    pivot_wider(names_from = analysis_type, 
                values_from = n_specific_pathways,
                names_prefix = "n_pathways_") %>%
    mutate(
      change = `n_pathways_Weighted (5%=50%, 10%=40%, 25%=10%)` - `n_pathways_Equal-weighted`,
      pct_change = round(change / `n_pathways_Equal-weighted` * 100, 1)
    )
  
  cat("Comparison of specific pathway counts:\n")
  print(comparison_table)
  cat("\n")
  
  write_csv(comparison_table, file.path(output_dir, "weighted_vs_equal_specificity_comparison.csv"))
  
} else {
  cat("Previous specificity results not found for comparison.\n\n")
}

# ------------------------------------------------------------------------------
# SAVE RESULTS
# ------------------------------------------------------------------------------

cat("Saving specificity analysis results...\n")

# Save all pathway analysis results
write_csv(pathway_analysis, file.path(output_dir, "weighted_pathway_specificity_analysis.csv"))

# Save specific pathways
write_csv(specific_pathways, file.path(output_dir, "weighted_highly_specific_pathways.csv"))

# Save cell type summary
write_csv(cell_type_summary, file.path(output_dir, "weighted_cell_type_summary.csv"))

# Save specificity summary
write_csv(specificity_summary, file.path(output_dir, "weighted_specificity_summary.csv"))

# Create detailed tables for each cell type
for(cell_type in unique(specific_pathways$cell_type_clean)) {
  cell_specific <- specific_pathways %>%
    filter(cell_type_clean == cell_type) %>%
    arrange(desc(abs(mean_nes)))
  
  filename <- paste0(gsub("[^A-Za-z0-9_]", "_", cell_type), "_weighted_specific_pathways.csv")
  write_csv(cell_specific, file.path(output_dir, filename))
  
  cat("Saved", nrow(cell_specific), "specific pathways for", cell_type, "\n")
}

# ------------------------------------------------------------------------------
# CREATE VISUALIZATION
# ------------------------------------------------------------------------------

cat("\nCreating specificity visualization...\n")

# Create specificity overview plot
specificity_plot <- ggplot(cell_type_summary, aes(x = reorder(cell_type_clean, n_specific_pathways), 
                                                  y = n_specific_pathways, 
                                                  fill = cell_type_clean)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_text(aes(label = paste0(n_specific_pathways, "\n(", pct_upregulated, "% up)")), 
            hjust = -0.1, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("CA3 Dorsal" = "#E31A1C", 
                              "Astrocytes" = "#FF7F00", 
                              "Sst+ HPF" = "#33A02C")) +
  coord_flip() +
  labs(
    title = "Cell Type-Specific Pathways: Weighted Meta-Analysis",
    subtitle = "Pathways enriched exclusively in each cell type (5%=50%, 10%=40%, 25%=10% weighting)",
    x = "Cell Type",
    y = "Number of Highly Specific Pathways",
    caption = "Excludes 50% gene lists | Weighted meta-analysis FDR < 0.05"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  ) +
  ylim(c(0, max(cell_type_summary$n_specific_pathways) * 1.2))

# Save plot
ggsave(file.path(output_dir, "weighted_cell_type_specificity_overview.png"), specificity_plot, 
       width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "weighted_cell_type_specificity_overview.pdf"), specificity_plot, 
       width = 10, height = 6, bg = "white")

cat("Saved: weighted_cell_type_specificity_overview.png/pdf\n")

# ------------------------------------------------------------------------------
# FINAL SUMMARY
# ------------------------------------------------------------------------------

cat("\n================================================================\n")
cat("WEIGHTED CELL TYPE SPECIFICITY ANALYSIS COMPLETE\n")
cat("================================================================\n\n")

cat("Key findings:\n")
cat("1. Highly specific pathways by cell type:\n")
for(i in 1:nrow(cell_type_summary)) {
  row <- cell_type_summary[i,]
  cat("   -", row$cell_type_clean, ":", row$n_specific_pathways, "pathways\n")
}

cat("\n2. Specificity distribution:\n")
for(i in 1:nrow(specificity_summary)) {
  row <- specificity_summary[i,]
  cat("   -", row$specificity_type, ":", row$n_pathways, "pathways (", row$percentage, "%)\n", sep = "")
}

cat("\n3. Weighting effect:\n")
total_specific <- sum(cell_type_summary$n_specific_pathways)
total_pathways <- nrow(meta_data)
specificity_ratio <- round(total_specific / total_pathways * 100, 1)
cat("   - Total specific pathways:", total_specific, "out of", total_pathways, "\n")
cat("   - Specificity ratio:", specificity_ratio, "%\n")

cat("\nFiles created in:", output_dir, "\n")
cat("- weighted_pathway_specificity_analysis.csv\n")
cat("- weighted_highly_specific_pathways.csv\n")
cat("- weighted_cell_type_summary.csv\n")
cat("- weighted_specificity_summary.csv\n")
if (file.exists(previous_specificity_file)) {
  cat("- weighted_vs_equal_specificity_comparison.csv\n")
}
cat("- weighted_cell_type_specificity_overview.png/pdf\n")
cat("- Individual cell type files\n\n")

cat("🎯 Weighted cell type specificity analysis complete!\n")
cat("📊 Results show effect of emphasizing specific gene lists (5% and 10%).\n")
cat("🔍 Compare with previous equal-weighted results to see the improvement.\n")