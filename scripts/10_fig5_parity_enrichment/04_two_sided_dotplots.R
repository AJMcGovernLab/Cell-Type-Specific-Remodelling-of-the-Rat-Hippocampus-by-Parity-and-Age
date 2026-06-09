#!/usr/bin/env Rscript

# ============================================================================
# Two-Sided Dot Plot for Pairwise Shared Pathways
# Creates mirrored dot plots showing pathway enrichment in two cell types
# ============================================================================

library(tidyverse)
library(ggplot2)
library(scales)

cat("================================================================\n")
cat("CREATING TWO-SIDED DOT PLOTS FOR PAIRWISE SHARED PATHWAYS\n")
cat("================================================================\n\n")

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------

# Paths via curated config + v7-adapter shim (02b_v7_adapter.R).
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))
input_dir  <- OUT$fig5ce
venn_dir   <- OUT$fig5fi
output_dir <- file.path(OUT$fig5ce, "two_sided_dotplots")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Input directory:", input_dir, "\n")
cat("Venn analysis directory:", venn_dir, "\n")
cat("Output directory:", output_dir, "\n\n")

# ------------------------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------------------------

cat("Loading weighted meta-analysis results and shared pathway data...\n")

# Load weighted meta-analysis results
meta_file <- file.path(input_dir, "weighted_meta_analysis_significant_v7.csv")
meta_data <- read_csv(meta_file, show_col_types = FALSE)

# Clean cell type names in meta_data
meta_data <- meta_data %>%
  mutate(
    cell_type_clean = case_when(
      meta_group == "CA3_do_combined" ~ "CA3 Dorsal",
      meta_group == "376_Astro" ~ "Astrocytes",
      meta_group == "78_Sst_HPF" ~ "Sst+ HPF",
      TRUE ~ meta_group
    )
  )

# Load shared pathway lists
ca3_astro_file <- file.path(venn_dir, "CA3_Dorsal_and_Astrocytes_pathways.csv")
ca3_sst_file <- file.path(venn_dir, "CA3_Dorsal_and_Sst_HPF_pathways.csv")
astro_sst_file <- file.path(venn_dir, "Astrocytes_and_Sst_HPF_pathways.csv")

ca3_astro_pathways <- read_csv(ca3_astro_file, show_col_types = FALSE)
ca3_sst_pathways <- read_csv(ca3_sst_file, show_col_types = FALSE)
astro_sst_pathways <- read_csv(astro_sst_file, show_col_types = FALSE)

cat("Shared pathway counts:\n")
cat("- CA3 Dorsal & Astrocytes:", nrow(ca3_astro_pathways), "pathways\n")
cat("- CA3 Dorsal & Sst+ HPF:", nrow(ca3_sst_pathways), "pathways\n")
cat("- Astrocytes & Sst+ HPF:", nrow(astro_sst_pathways), "pathways\n\n")

# ------------------------------------------------------------------------------
# FUNCTION TO CREATE TWO-SIDED DOT PLOT
# ------------------------------------------------------------------------------

create_two_sided_dotplot <- function(pathway_list, cell_type1, cell_type2, plot_title) {
  
  cat("Creating two-sided dot plot for", cell_type1, "and", cell_type2, "\n")
  
  # Get data for both cell types
  cell1_data <- meta_data %>%
    filter(ID %in% pathway_list$pathway_id, cell_type_clean == cell_type1) %>%
    select(ID, Description, database, 
           NES = weighted_mean_NES, 
           FDR = weighted_meta_padj_corrected) %>%
    mutate(cell_type = cell_type1)
  
  cell2_data <- meta_data %>%
    filter(ID %in% pathway_list$pathway_id, cell_type_clean == cell_type2) %>%
    select(ID, Description, database, 
           NES = weighted_mean_NES, 
           FDR = weighted_meta_padj_corrected) %>%
    mutate(cell_type = cell_type2)
  
  # Combine data
  plot_data <- bind_rows(cell1_data, cell2_data)
  
  if(nrow(plot_data) == 0) {
    cat("No data found for this comparison, skipping...\n")
    return(NULL)
  }
  
  # Calculate average FDR for ordering pathways
  pathway_order <- plot_data %>%
    group_by(ID, Description) %>%
    summarise(
      avg_neg_log10_fdr = mean(-log10(FDR + 1e-300)),
      .groups = 'drop'
    ) %>%
    arrange(desc(avg_neg_log10_fdr))
  
  # Prepare plot data with transformations
  plot_data <- plot_data %>%
    mutate(
      # Transform FDR to -log10
      neg_log10_fdr = -log10(FDR + 1e-300),
      
      # Mirror the FDR values for left/right placement
      x_position = ifelse(cell_type == cell_type1, 
                         -neg_log10_fdr,  # Left side (negative)
                         neg_log10_fdr),   # Right side (positive)
      
      # Determine regulation direction
      regulation = case_when(
        NES > 0 ~ "Upregulated",
        NES < 0 ~ "Downregulated",
        TRUE ~ "Neutral"
      ),
      
      # Absolute NES for size
      abs_NES = abs(NES),
      
      # Create pathway label (truncate if too long)
      pathway_label = str_trunc(Description, 50, "right")
    ) %>%
    # Order pathways by average significance
    mutate(pathway_label = factor(pathway_label, 
                                  levels = pathway_order$Description %>% 
                                           str_trunc(50, "right") %>% 
                                           rev()))
  
  # Determine x-axis limits
  max_fdr <- max(abs(plot_data$x_position), na.rm = TRUE)
  x_limits <- c(-max_fdr * 1.1, max_fdr * 1.1)
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = x_position, y = pathway_label)) +
    
    # Add vertical line at 0
    geom_vline(xintercept = 0, color = "gray50", linewidth = 0.5) +
    
    # Add dots
    geom_point(aes(size = abs_NES, color = regulation), alpha = 0.8) +
    
    # Color scale
    scale_color_manual(
      values = c("Upregulated" = "#E31A1C", 
                "Downregulated" = "#1F78B4",
                "Neutral" = "gray50"),
      name = "Regulation"
    ) +
    
    # Size scale
    scale_size_continuous(
      range = c(2, 8),
      breaks = c(1, 2, 3, 4, 5),
      name = "Absolute NES"
    ) +
    
    # X-axis scale with custom labels
    scale_x_continuous(
      limits = x_limits,
      breaks = seq(-10, 10, by = 2),
      labels = function(x) sprintf("%.0f", abs(x)),
      expand = c(0.02, 0)
    ) +
    
    # Labels and theme
    labs(
      title = plot_title,
      subtitle = paste0("Shared pathways between ", cell_type1, " and ", cell_type2),
      x = "-log10(FDR)",
      y = NULL,
      caption = paste0(cell_type1, " (left) | ", cell_type2, " (right)")
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0.5, color = "gray50"),
      axis.title.x = element_text(size = 11),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 8),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
    )
  
  return(p)
}

# ------------------------------------------------------------------------------
# CREATE PLOTS FOR EACH PAIRWISE COMPARISON
# ------------------------------------------------------------------------------

cat("\nGenerating two-sided dot plots...\n")

# 1. CA3 Dorsal & Astrocytes
p1 <- create_two_sided_dotplot(
  ca3_astro_pathways,
  "CA3 Dorsal", 
  "Astrocytes",
  "Pairwise Shared Pathways: CA3 Dorsal & Astrocytes"
)

if(!is.null(p1)) {
  # Adjust height based on number of pathways
  plot_height <- max(6, min(12, nrow(ca3_astro_pathways) * 0.3))
  
  ggsave(
    filename = file.path(output_dir, "CA3_Astro_two_sided_dotplot.pdf"),
    plot = p1,
    width = 10,
    height = plot_height,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(output_dir, "CA3_Astro_two_sided_dotplot.png"),
    plot = p1,
    width = 10,
    height = plot_height,
    dpi = 300
  )
  
  cat("Saved CA3 Dorsal & Astrocytes plot\n")
}

# 2. CA3 Dorsal & Sst+ HPF
p2 <- create_two_sided_dotplot(
  ca3_sst_pathways,
  "CA3 Dorsal",
  "Sst+ HPF",
  "Pairwise Shared Pathways: CA3 Dorsal & Sst+ HPF"
)

if(!is.null(p2)) {
  plot_height <- max(6, min(12, nrow(ca3_sst_pathways) * 0.3))
  
  ggsave(
    filename = file.path(output_dir, "CA3_Sst_two_sided_dotplot.pdf"),
    plot = p2,
    width = 10,
    height = plot_height,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(output_dir, "CA3_Sst_two_sided_dotplot.png"),
    plot = p2,
    width = 10,
    height = plot_height,
    dpi = 300
  )
  
  cat("Saved CA3 Dorsal & Sst+ HPF plot\n")
}

# 3. Astrocytes & Sst+ HPF - Limited to top 10 up and top 10 down
# Since this has many pathways (110), we'll create a filtered version
cat("\nFiltering Astrocytes & Sst+ HPF pathways to top 10 up and top 10 down...\n")

# Get data for both cell types for filtering
astro_data_filter <- meta_data %>%
  filter(ID %in% astro_sst_pathways$pathway_id, cell_type_clean == "Astrocytes")

sst_data_filter <- meta_data %>%
  filter(ID %in% astro_sst_pathways$pathway_id, cell_type_clean == "Sst+ HPF")

# Combine and calculate average metrics for ranking
combined_metrics <- bind_rows(astro_data_filter, sst_data_filter) %>%
  group_by(ID, Description, database) %>%
  summarise(
    avg_NES = mean(weighted_mean_NES),
    avg_neg_log10_fdr = mean(-log10(weighted_meta_padj_corrected + 1e-300)),
    is_significant = all(weighted_meta_padj_corrected < 0.05),
    .groups = 'drop'
  ) %>%
  filter(is_significant)  # Only keep pathways significant in both cell types

# Get top 10 upregulated and top 10 downregulated
top_upregulated <- combined_metrics %>%
  filter(avg_NES > 0) %>%
  arrange(desc(avg_neg_log10_fdr), desc(avg_NES)) %>%
  slice_head(n = 10)

top_downregulated <- combined_metrics %>%
  filter(avg_NES < 0) %>%
  arrange(desc(avg_neg_log10_fdr), desc(abs(avg_NES))) %>%
  slice_head(n = 10)

# Combine the filtered pathways (rename ID -> pathway_id to match other inputs)
filtered_pathways <- bind_rows(top_upregulated, top_downregulated) %>%
  select(pathway_id = ID, Description, database)

cat("Selected", nrow(top_upregulated), "upregulated and", nrow(top_downregulated), "downregulated pathways\n")

# Create plot with filtered pathways
p3 <- create_two_sided_dotplot(
  filtered_pathways,
  "Astrocytes",
  "Sst+ HPF",
  "Top Shared Pathways: Astrocytes & Sst+ HPF (Top 10 Up/Down)"
)

if(!is.null(p3)) {
  plot_height <- 8  # Fixed height for ~20 pathways
  
  ggsave(
    filename = file.path(output_dir, "Astro_Sst_two_sided_dotplot.pdf"),
    plot = p3,
    width = 10,
    height = plot_height,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(output_dir, "Astro_Sst_two_sided_dotplot.png"),
    plot = p3,
    width = 10,
    height = plot_height,
    dpi = 300
  )
  
  cat("Saved Astrocytes & Sst+ HPF plot (filtered)\n")
}

# ------------------------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------------------------

cat("\n================================================================\n")
cat("TWO-SIDED DOT PLOTS COMPLETED\n")
cat("================================================================\n")
cat("Output directory:", output_dir, "\n")
cat("Files created:\n")
system(paste("ls -la", output_dir))
cat("\nAll plots have been generated successfully!\n")