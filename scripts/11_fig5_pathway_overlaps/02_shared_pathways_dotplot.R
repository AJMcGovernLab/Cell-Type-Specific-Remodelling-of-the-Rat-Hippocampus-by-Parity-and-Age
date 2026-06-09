#!/usr/bin/env Rscript

# ============================================================================
# Dotplot for Top 20 Pathways Shared Across All Three Cell Types
# Creates a dotplot showing pathways present in all cell types
# ============================================================================

library(tidyverse)
library(ggplot2)
library(viridis)

cat("================================================================\n")
cat("CREATING DOTPLOT FOR PATHWAYS SHARED ACROSS ALL CELL TYPES\n")
cat("================================================================\n\n")

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------

# Paths via curated config + v7-adapter shim (02b_v7_adapter.R).
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))
input_dir  <- OUT$fig5ce
venn_dir   <- OUT$fig5fi
output_dir <- OUT$fig5fi

# Parameters
top_n <- 20  # Number of top pathways to display

cat("Input directory:", input_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Displaying top", top_n, "pathways\n\n")

# ------------------------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------------------------

cat("Loading weighted meta-analysis results...\n")

# Load weighted meta-analysis results
meta_file <- file.path(input_dir, "weighted_meta_analysis_significant_v7.csv")
meta_data <- read_csv(meta_file, show_col_types = FALSE)

# Clean cell type names
meta_data <- meta_data %>%
  mutate(
    cell_type_clean = case_when(
      meta_group == "CA3_do_combined" ~ "CA3 Dorsal",
      meta_group == "376_Astro" ~ "Astrocytes",
      meta_group == "78_Sst_HPF" ~ "Sst+ HPF",
      TRUE ~ meta_group
    )
  )

# Load pathway lists to identify shared pathways
ca3_only_file <- file.path(venn_dir, "CA3_Dorsal_only_pathways.csv")
astro_only_file <- file.path(venn_dir, "Astrocytes_only_pathways.csv")
sst_only_file <- file.path(venn_dir, "Sst_HPF_only_pathways.csv")
ca3_astro_file <- file.path(venn_dir, "CA3_Dorsal_and_Astrocytes_pathways.csv")
ca3_sst_file <- file.path(venn_dir, "CA3_Dorsal_and_Sst_HPF_pathways.csv")
astro_sst_file <- file.path(venn_dir, "Astrocytes_and_Sst_HPF_pathways.csv")

# Read all pathway lists (column is named `pathway_id` in the source CSVs).
read_ids <- function(path) {
  if (!file.exists(path)) return(character(0))
  df <- read_csv(path, show_col_types = FALSE)
  col <- intersect(c("pathway_id", "ID"), colnames(df))[1]
  if (is.na(col)) character(0) else df[[col]]
}
ca3_only   <- read_ids(ca3_only_file)
astro_only <- read_ids(astro_only_file)
sst_only   <- read_ids(sst_only_file)
ca3_astro  <- read_ids(ca3_astro_file)
ca3_sst    <- read_ids(ca3_sst_file)
astro_sst  <- read_ids(astro_sst_file)

# Combine all pathways that are NOT shared by all three
not_all_three <- c(ca3_only, astro_only, sst_only, ca3_astro, ca3_sst, astro_sst)

cat("Identified", length(not_all_three), "pathways that are NOT shared by all three cell types\n")

# ------------------------------------------------------------------------------
# IDENTIFY PATHWAYS SHARED BY ALL THREE CELL TYPES
# ------------------------------------------------------------------------------

cat("\nIdentifying pathways shared by all three cell types...\n")

# Get pathways that are in the data but NOT in the exclusion list
# Column in source CSV is `databases` (plural)
shared_by_all <- meta_data %>%
  filter(!ID %in% not_all_three) %>%
  select(ID, Description, database, cell_type_clean,
         weighted_mean_NES, weighted_meta_padj_corrected) %>%
  distinct()

# Check which pathways appear in all three cell types
pathway_counts <- shared_by_all %>%
  group_by(ID, Description, database) %>%
  summarise(
    n_cell_types = n_distinct(cell_type_clean),
    .groups = 'drop'
  )

all_three_pathways <- pathway_counts %>%
  filter(n_cell_types == 3) %>%
  pull(ID)

cat("Found", length(all_three_pathways), "pathways shared by all three cell types\n")

# ------------------------------------------------------------------------------
# PREPARE DATA FOR DOTPLOT
# ------------------------------------------------------------------------------

cat("\nPreparing data for dotplot...\n")

# Get data for all three cell types for these pathways
dotplot_data <- meta_data %>%
  filter(ID %in% all_three_pathways) %>%
  select(ID, Description, database, cell_type_clean,
         NES = weighted_mean_NES,
         FDR = weighted_meta_padj_corrected) %>%
  mutate(
    # Transform FDR for visualization
    neg_log10_fdr = -log10(FDR + 1e-300),
    
    # Determine regulation direction
    regulation = case_when(
      NES > 0 ~ "Up",
      NES < 0 ~ "Down",
      TRUE ~ "Neutral"
    )
  )

# Calculate average significance across cell types for ranking
pathway_ranking <- dotplot_data %>%
  group_by(ID, Description) %>%
  summarise(
    avg_neg_log10_fdr = mean(neg_log10_fdr),
    max_abs_nes = max(abs(NES)),
    .groups = 'drop'
  ) %>%
  arrange(desc(avg_neg_log10_fdr), desc(max_abs_nes))

# Select top N pathways
top_pathways <- pathway_ranking %>%
  slice_head(n = top_n) %>%
  pull(ID)

# Filter data to top pathways
plot_data <- dotplot_data %>%
  filter(ID %in% top_pathways) %>%
  mutate(
    # Create pathway label (truncate if needed)
    pathway_label = str_trunc(Description, 50, "right"),
    # Order pathways by average significance
    pathway_label = factor(pathway_label,
                          levels = pathway_ranking %>%
                            filter(ID %in% top_pathways) %>%
                            mutate(label = str_trunc(Description, 50, "right")) %>%
                            pull(label) %>%
                            rev())
  )

# ------------------------------------------------------------------------------
# CREATE DOTPLOT
# ------------------------------------------------------------------------------

cat("\nCreating dotplot...\n")

# Create the dotplot
p <- ggplot(plot_data, aes(x = cell_type_clean, y = pathway_label)) +
  
  # Add dots - size by significance, color by NES
  geom_point(aes(size = neg_log10_fdr, color = NES), alpha = 0.9) +
  
  # Color scale for NES
  scale_color_gradient2(
    low = "#1F78B4",      # Blue for negative
    mid = "white",
    high = "#E31A1C",     # Red for positive
    midpoint = 0,
    limits = c(-max(abs(plot_data$NES)), max(abs(plot_data$NES))),
    name = "NES"
  ) +
  
  # Size scale for significance
  scale_size_continuous(
    range = c(3, 10),
    name = "-log10(FDR)",
    breaks = c(5, 10, 15, 20)
  ) +
  
  # Theme and labels
  labs(
    title = "Top 20 Pathways Shared Across All Three Cell Types",
    subtitle = "Pathways significantly enriched in CA3 Dorsal, Astrocytes, and Sst+ HPF",
    x = NULL,
    y = NULL,
    caption = "Size: -log10(FDR) | Color: Normalized Enrichment Score (NES)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
    plot.caption = element_text(size = 9, hjust = 0.5, color = "gray50"),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "gray95", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

# Save the plot
ggsave(
  filename = file.path(output_dir, "top20_shared_pathways_dotplot.pdf"),
  plot = p,
  width = 10,
  height = 10,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "top20_shared_pathways_dotplot.png"),
  plot = p,
  width = 10,
  height = 10,
  dpi = 300
)

cat("Saved dotplot\n")

# ------------------------------------------------------------------------------
# CREATE ALTERNATIVE HEATMAP-STYLE VISUALIZATION
# ------------------------------------------------------------------------------

cat("\nCreating alternative heatmap-style visualization...\n")

# Prepare data for heatmap
heatmap_data <- plot_data %>%
  select(pathway_label, cell_type_clean, NES, neg_log10_fdr)

# Create heatmap-style plot
p_heatmap <- ggplot(heatmap_data, aes(x = cell_type_clean, y = pathway_label)) +
  
  # Add tiles colored by NES
  geom_tile(aes(fill = NES), color = "white", linewidth = 0.5) +
  
  # Add text labels showing FDR
  geom_text(aes(label = sprintf("%.1f", neg_log10_fdr)), 
            size = 3, color = "black") +
  
  # Color scale
  scale_fill_gradient2(
    low = "#1F78B4",
    mid = "white", 
    high = "#E31A1C",
    midpoint = 0,
    limits = c(-max(abs(heatmap_data$NES)), max(abs(heatmap_data$NES))),
    name = "NES"
  ) +
  
  # Theme and labels
  labs(
    title = "Top 20 Pathways Shared Across All Three Cell Types",
    subtitle = "Heatmap showing NES values with -log10(FDR) labels",
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

# Save heatmap
ggsave(
  filename = file.path(output_dir, "top20_shared_pathways_heatmap.pdf"),
  plot = p_heatmap,
  width = 10,
  height = 10,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "top20_shared_pathways_heatmap.png"),
  plot = p_heatmap,
  width = 10,
  height = 10,
  dpi = 300
)

cat("Saved heatmap visualization\n")

# ------------------------------------------------------------------------------
# SUMMARY STATISTICS
# ------------------------------------------------------------------------------

cat("\n================================================================\n")
cat("SUMMARY OF SHARED PATHWAYS\n")
cat("================================================================\n")

# Display top 20 pathways
cat("\nTop 20 pathways shared across all three cell types:\n")
cat("----------------------------------------------------\n")

display_pathways <- pathway_ranking %>%
  filter(ID %in% top_pathways) %>%
  slice_head(n = 20) %>%
  mutate(
    rank = row_number(),
    pathway = str_trunc(Description, 60, "right")
  )

for(i in 1:nrow(display_pathways)) {
  cat(sprintf("%2d. %s (avg -log10 FDR: %.2f)\n", 
              display_pathways$rank[i],
              display_pathways$pathway[i],
              display_pathways$avg_neg_log10_fdr[i]))
}

# Overall statistics
cat("\n----------------------------------------------------\n")
cat("Total pathways shared by all three cell types:", length(all_three_pathways), "\n")
cat("Pathways displayed in plot:", top_n, "\n")

cat("\n================================================================\n")
cat("DOTPLOTS COMPLETED\n")
cat("================================================================\n")
cat("Output files:\n")
cat("  - top20_shared_pathways_dotplot.pdf/png\n")
cat("  - top20_shared_pathways_heatmap.pdf/png\n")
cat("\nAll plots have been generated successfully!\n")