#!/usr/bin/env Rscript

# ============================================================================
# Percentage-Based Venn Diagrams for Individual Enrichment Results
# Create separate Venn diagrams for 5%, 10%, and 25% gene lists
# ============================================================================

library(tidyverse)
library(VennDiagram)
library(RColorBrewer)
library(gridExtra)
library(grid)

cat("================================================================\n")
cat("CREATING PERCENTAGE-BASED VENN DIAGRAMS\n")
cat("================================================================\n\n")

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------

# Paths via curated config. NOTE: this script reads the v6 percentage-
# weighted enrichment file. Under the curated GSEA pipeline that input
# does not exist; this script is preserved for legacy reference but
# will short-circuit if its input file is missing.
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))
input_file <- file.path(CHECKPOINT_DIR, "enrichment_parity",
                        "individual_enhanced_results_v6.csv")
output_dir <- file.path(OUT$fig5fi, "percentage_venns")
if (!file.exists(input_file)) {
  message(sprintf("[skip] %s not staged; percentage Venns not regenerated.",
                  input_file))
  if (sys.nframe() == 0L) quit(save = "no", status = 0)
}

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
significance_threshold <- 0.05
percentages_to_analyze <- c("5pct", "10pct", "25pct")

cat("Input file:", input_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Percentages to analyze:", paste(percentages_to_analyze, collapse = ", "), "\n")
cat("Significance threshold:", significance_threshold, "\n\n")

# ------------------------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------------------------

cat("Loading individual enrichment results...\n")

# Load individual enrichment results
individual_results <- read_csv(input_file, show_col_types = FALSE)

cat("Data loaded:\n")
cat("- Total individual results:", nrow(individual_results), "\n")
cat("- Cell types:", paste(unique(individual_results$cell_type), collapse = ", "), "\n")
cat("- Percentages:", paste(unique(individual_results$percentage), collapse = ", "), "\n")
cat("- Databases:", paste(unique(individual_results$database), collapse = ", "), "\n\n")

# Filter for significant results only
significant_results <- individual_results %>%
  filter(p.adjust < significance_threshold)

cat("After filtering for significance (FDR < ", significance_threshold, "):\n")
cat("- Significant results:", nrow(significant_results), "\n\n")

# ------------------------------------------------------------------------------
# DEFINE CELL TYPE MAPPING
# ------------------------------------------------------------------------------

# Map individual cell types to meta-groups
cell_type_mapping <- list(
  "CA3_do_combined" = c("356_CA3_do", "358_CA3_do"),
  "376_Astro" = c("376_Astro"),
  "78_Sst_HPF" = c("78_Sst_HPF")
)

# Add meta-group column
significant_results <- significant_results %>%
  mutate(
    meta_group = case_when(
      cell_type %in% cell_type_mapping$CA3_do_combined ~ "CA3_do_combined",
      cell_type %in% cell_type_mapping$`376_Astro` ~ "376_Astro",
      cell_type %in% cell_type_mapping$`78_Sst_HPF` ~ "78_Sst_HPF",
      TRUE ~ cell_type
    ),
    cell_type_clean = case_when(
      meta_group == "CA3_do_combined" ~ "CA3 Dorsal",
      meta_group == "376_Astro" ~ "Astrocytes",
      meta_group == "78_Sst_HPF" ~ "Sst+ HPF",
      TRUE ~ meta_group
    )
  )

# ------------------------------------------------------------------------------
# FUNCTION TO CREATE VENN DIAGRAM FOR EACH PERCENTAGE
# ------------------------------------------------------------------------------

create_percentage_venn <- function(data, percentage_name) {
  
  cat("Creating Venn diagram for", percentage_name, "gene lists...\n")
  
  # Filter data for this percentage
  pct_data <- data %>%
    filter(percentage == percentage_name)
  
  if(nrow(pct_data) == 0) {
    cat("No data found for", percentage_name, ", skipping...\n")
    return(NULL)
  }
  
  # Get unique pathways for each meta-group
  ca3_pathways <- pct_data %>%
    filter(meta_group == "CA3_do_combined") %>%
    pull(ID) %>%
    unique()
  
  astro_pathways <- pct_data %>%
    filter(meta_group == "376_Astro") %>%
    pull(ID) %>%
    unique()
  
  sst_pathways <- pct_data %>%
    filter(meta_group == "78_Sst_HPF") %>%
    pull(ID) %>%
    unique()
  
  cat("  Pathway counts for", percentage_name, ":\n")
  cat("  - CA3 Dorsal:", length(ca3_pathways), "pathways\n")
  cat("  - Astrocytes:", length(astro_pathways), "pathways\n")
  cat("  - Sst+ HPF:", length(sst_pathways), "pathways\n")
  
  # Calculate overlaps
  ca3_only <- setdiff(setdiff(ca3_pathways, astro_pathways), sst_pathways)
  astro_only <- setdiff(setdiff(astro_pathways, ca3_pathways), sst_pathways)
  sst_only <- setdiff(setdiff(sst_pathways, ca3_pathways), astro_pathways)
  
  ca3_astro <- setdiff(intersect(ca3_pathways, astro_pathways), sst_pathways)
  ca3_sst <- setdiff(intersect(ca3_pathways, sst_pathways), astro_pathways)
  astro_sst <- setdiff(intersect(astro_pathways, sst_pathways), ca3_pathways)
  
  all_three <- intersect(intersect(ca3_pathways, astro_pathways), sst_pathways)
  
  # Create overlap summary for this percentage
  overlap_summary <- tibble(
    percentage = percentage_name,
    region = c(
      "CA3 Dorsal only",
      "Astrocytes only", 
      "Sst+ HPF only",
      "CA3 Dorsal & Astrocytes",
      "CA3 Dorsal & Sst+ HPF",
      "Astrocytes & Sst+ HPF",
      "All three cell types"
    ),
    count = c(
      length(ca3_only),
      length(astro_only),
      length(sst_only),
      length(ca3_astro),
      length(ca3_sst),
      length(astro_sst),
      length(all_three)
    )
  )
  
  total_unique <- length(unique(c(ca3_pathways, astro_pathways, sst_pathways)))
  overlap_summary$percentage_of_total <- round(overlap_summary$count / total_unique * 100, 1)
  
  cat("  Overlap summary for", percentage_name, ":\n")
  print(overlap_summary)
  cat("\n")
  
  # Set colors for each cell type
  venn_colors <- c("#E31A1C", "#FF7F00", "#33A02C")  # CA3, Astro, Sst+ respectively
  
  # Create percentage-specific title
  pct_clean <- gsub("pct", "%", percentage_name)
  
  # Create the Venn diagram
  venn_plot <- venn.diagram(
    x = list(
      "CA3 Dorsal" = ca3_pathways,
      "Astrocytes" = astro_pathways,
      "Sst+ HPF" = sst_pathways
    ),
    category.names = c(
      paste0("CA3 Dorsal\n(", length(ca3_pathways), " pathways)"), 
      paste0("Astrocytes\n(", length(astro_pathways), " pathways)"), 
      paste0("Sst+ HPF\n(", length(sst_pathways), " pathways)")
    ),
    filename = NULL,
    
    # Colors
    fill = venn_colors,
    alpha = 0.6,
    
    # Circle properties
    lwd = 2.5,
    col = "black",
    
    # Text properties
    cex = 1.4,
    fontfamily = "sans",
    fontface = "bold",
    
    # Category label properties
    cat.cex = 1.1,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    
    # Main title
    main = paste0("Pathway Enrichment: ", pct_clean, " Gene Lists"),
    main.cex = 1.6,
    main.fontface = "bold",
    main.pos = c(0.5, 0.95),
    
    # Subtitle
    sub = paste0("Individual cell type enrichment | Total unique pathways: ", total_unique, " | FDR < 0.05"),
    sub.cex = 1.0,
    sub.pos = c(0.5, 0.05),
    
    # Margins
    margin = 0.2,
    
    # Disable log file
    disable.logging = TRUE
  )
  
  # Save the Venn diagram
  pct_filename <- tolower(gsub("pct", "_percent", percentage_name))
  
  png(file.path(output_dir, paste0("pathway_venn_", pct_filename, ".png")), 
      width = 12, height = 10, units = "in", res = 300, bg = "white")
  grid.draw(venn_plot)
  dev.off()
  
  pdf(file.path(output_dir, paste0("pathway_venn_", pct_filename, ".pdf")), 
      width = 12, height = 10, bg = "white")
  grid.draw(venn_plot)
  dev.off()
  
  cat("  Saved: pathway_venn_", pct_filename, ".png/pdf\n\n")
  
  return(overlap_summary)
}

# ------------------------------------------------------------------------------
# CREATE VENN DIAGRAMS FOR EACH PERCENTAGE
# ------------------------------------------------------------------------------

cat("Creating Venn diagrams for each percentage...\n\n")

# Store overlap summaries for comparison
all_overlaps <- list()

for(pct in percentages_to_analyze) {
  overlap_data <- create_percentage_venn(significant_results, pct)
  if(!is.null(overlap_data)) {
    all_overlaps[[pct]] <- overlap_data
  }
}

# ------------------------------------------------------------------------------
# CREATE COMPARISON ANALYSIS
# ------------------------------------------------------------------------------

cat("Creating comparison analysis across percentages...\n")

if(length(all_overlaps) > 0) {
  
  # Combine all overlap data
  combined_overlaps <- bind_rows(all_overlaps)
  
  # Create comparison table
  comparison_table <- combined_overlaps %>%
    select(percentage, region, count, percentage_of_total) %>%
    pivot_wider(names_from = percentage, 
                values_from = c(count, percentage_of_total),
                names_sep = "_") %>%
    arrange(region)
  
  cat("Comparison across percentages:\n")
  print(comparison_table)
  cat("\n")
  
  # Save comparison table
  write_csv(comparison_table, file.path(output_dir, "percentage_comparison_summary.csv"))
  write_csv(combined_overlaps, file.path(output_dir, "all_percentage_overlaps.csv"))
  
  # Create trend analysis
  trend_analysis <- combined_overlaps %>%
    group_by(region) %>%
    summarise(
      pct_5_count = count[percentage == "5pct"],
      pct_10_count = count[percentage == "10pct"],
      pct_25_count = count[percentage == "25pct"],
      trend_5_to_25 = pct_25_count - pct_5_count,
      pct_change = round((pct_25_count - pct_5_count) / pmax(pct_5_count, 1) * 100, 1),
      .groups = "drop"
    ) %>%
    arrange(desc(abs(trend_5_to_25)))
  
  cat("Trend analysis (5% to 25% change):\n")
  print(trend_analysis)
  cat("\n")
  
  write_csv(trend_analysis, file.path(output_dir, "percentage_trend_analysis.csv"))
  
} else {
  cat("No overlap data available for comparison.\n")
}

# ------------------------------------------------------------------------------
# CREATE SUMMARY STATISTICS
# ------------------------------------------------------------------------------

cat("Creating summary statistics by percentage...\n")

# Summary statistics for each percentage
percentage_stats <- significant_results %>%
  filter(percentage %in% percentages_to_analyze) %>%
  group_by(percentage, meta_group) %>%
  summarise(
    n_pathways = n_distinct(ID),
    mean_nes = round(mean(NES, na.rm = TRUE), 3),
    median_fdr = round(median(p.adjust, na.rm = TRUE), 6),
    n_databases = n_distinct(database),
    .groups = "drop"
  ) %>%
  arrange(percentage, desc(n_pathways))

cat("Pathway statistics by percentage and cell type:\n")
print(percentage_stats)
cat("\n")

write_csv(percentage_stats, file.path(output_dir, "percentage_pathway_statistics.csv"))

# Overall summary
overall_summary <- significant_results %>%
  filter(percentage %in% percentages_to_analyze) %>%
  group_by(percentage) %>%
  summarise(
    total_pathways = n(),
    unique_pathways = n_distinct(ID),
    n_cell_types = n_distinct(meta_group),
    mean_nes = round(mean(NES, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(percentage)

cat("Overall summary by percentage:\n")
print(overall_summary)
cat("\n")

write_csv(overall_summary, file.path(output_dir, "percentage_overall_summary.csv"))

# ------------------------------------------------------------------------------
# FINAL SUMMARY
# ------------------------------------------------------------------------------

cat("================================================================\n")
cat("PERCENTAGE-BASED VENN DIAGRAM ANALYSIS COMPLETE\n")
cat("================================================================\n\n")

cat("Venn diagrams created for percentages:\n")
for(pct in percentages_to_analyze) {
  pct_clean <- gsub("pct", "%", pct)
  pct_filename <- tolower(gsub("pct", "_percent", pct))
  cat("- ", pct_clean, " gene lists: pathway_venn_", pct_filename, ".png/pdf\n", sep = "")
}

cat("\nKey insights:\n")
if(exists("trend_analysis") && nrow(trend_analysis) > 0) {
  # Find the most interesting trends
  biggest_increase <- trend_analysis %>% slice_max(trend_5_to_25, n = 1)
  biggest_decrease <- trend_analysis %>% slice_min(trend_5_to_25, n = 1)
  
  cat("- Biggest increase from 5% to 25%:", biggest_increase$region, "(+", biggest_increase$trend_5_to_25, " pathways)\n")
  cat("- Biggest decrease from 5% to 25%:", biggest_decrease$region, "(", biggest_decrease$trend_5_to_25, " pathways)\n")
}

cat("\nFiles created in:", output_dir, "\n")
cat("- pathway_venn_[percentage].png/pdf (individual Venn diagrams)\n")
cat("- percentage_comparison_summary.csv (cross-percentage comparison)\n")
cat("- percentage_trend_analysis.csv (trend analysis)\n")
cat("- percentage_pathway_statistics.csv (detailed statistics)\n")
cat("- all_percentage_overlaps.csv (complete overlap data)\n\n")

cat("🎯 Percentage-based Venn analysis complete!\n")
cat("📊 Shows how pathway overlap changes with gene list specificity.\n")
cat("🔍 Higher percentages likely show more shared/core effects.\n")
cat("🎲 Lower percentages likely show more cell type-specific effects.\n")