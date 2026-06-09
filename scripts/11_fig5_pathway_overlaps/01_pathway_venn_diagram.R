#!/usr/bin/env Rscript

# ============================================================================
# Pathway Venn Diagram for Weighted Meta-Analysis Results
# Create Venn diagram showing pathway overlap between cell types
# ============================================================================

library(tidyverse)
library(VennDiagram)
library(RColorBrewer)
library(gridExtra)
library(grid)

cat("================================================================\n")
cat("CREATING PATHWAY VENN DIAGRAM\n")
cat("================================================================\n\n")

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------

# Paths via curated config + v7-adapter shim (02b_v7_adapter.R).
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))
input_dir  <- OUT$fig5ce
output_dir <- OUT$fig5fi

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Input directory:", input_dir, "\n")
cat("Output directory:", output_dir, "\n\n")

# ------------------------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------------------------

cat("Loading weighted meta-analysis results...\n")

# Load weighted meta-analysis results
meta_file <- file.path(input_dir, "weighted_meta_analysis_significant_v7.csv")
meta_data <- read_csv(meta_file, show_col_types = FALSE)

cat("Data loaded:\n")
cat("- Total significant pathways:", nrow(meta_data), "\n")
cat("- Cell types found:", paste(unique(meta_data$meta_group), collapse = ", "), "\n")
cat("- Databases:", paste(unique(meta_data$database), collapse = ", "), "\n\n")

# ------------------------------------------------------------------------------
# EXTRACT PATHWAY LISTS BY CELL TYPE
# ------------------------------------------------------------------------------

cat("Extracting pathway lists for each cell type...\n")

# Create pathway lists for each cell type
ca3_pathways <- meta_data %>%
  filter(meta_group == "CA3_do_combined") %>%
  pull(ID) %>%
  unique()

astro_pathways <- meta_data %>%
  filter(meta_group == "376_Astro") %>%
  pull(ID) %>%
  unique()

sst_pathways <- meta_data %>%
  filter(meta_group == "78_Sst_HPF") %>%
  pull(ID) %>%
  unique()

cat("Pathway counts by cell type:\n")
cat("- CA3 Dorsal:", length(ca3_pathways), "pathways\n")
cat("- Astrocytes:", length(astro_pathways), "pathways\n")
cat("- Sst+ HPF:", length(sst_pathways), "pathways\n\n")

# ------------------------------------------------------------------------------
# CALCULATE OVERLAPS
# ------------------------------------------------------------------------------

cat("Calculating pathway overlaps...\n")

# Calculate all possible intersections
ca3_only <- setdiff(setdiff(ca3_pathways, astro_pathways), sst_pathways)
astro_only <- setdiff(setdiff(astro_pathways, ca3_pathways), sst_pathways)
sst_only <- setdiff(setdiff(sst_pathways, ca3_pathways), astro_pathways)

ca3_astro <- setdiff(intersect(ca3_pathways, astro_pathways), sst_pathways)
ca3_sst <- setdiff(intersect(ca3_pathways, sst_pathways), astro_pathways)
astro_sst <- setdiff(intersect(astro_pathways, sst_pathways), ca3_pathways)

all_three <- intersect(intersect(ca3_pathways, astro_pathways), sst_pathways)

# Create overlap summary
overlap_summary <- tibble(
  Region = c(
    "CA3 Dorsal only",
    "Astrocytes only", 
    "Sst+ HPF only",
    "CA3 Dorsal & Astrocytes",
    "CA3 Dorsal & Sst+ HPF",
    "Astrocytes & Sst+ HPF",
    "All three cell types"
  ),
  Count = c(
    length(ca3_only),
    length(astro_only),
    length(sst_only),
    length(ca3_astro),
    length(ca3_sst),
    length(astro_sst),
    length(all_three)
  ),
  Percentage = round(Count / nrow(meta_data) * 100, 1)
)

cat("Pathway overlap summary:\n")
print(overlap_summary)
cat("\n")

# ------------------------------------------------------------------------------
# CREATE VENN DIAGRAM
# ------------------------------------------------------------------------------

cat("Creating Venn diagram...\n")

# Set colors for each cell type
venn_colors <- c("#E31A1C", "#FF7F00", "#33A02C")  # CA3, Astro, Sst+ respectively

# Create the Venn diagram
venn_plot <- venn.diagram(
  x = list(
    "CA3 Dorsal" = ca3_pathways,
    "Astrocytes" = astro_pathways,
    "Sst+ HPF" = sst_pathways
  ),
  category.names = c("CA3 Dorsal\n(CA3_do_combined)", "Astrocytes\n(376_Astro)", "Sst+ HPF\n(78_Sst_HPF)"),
  filename = NULL,  # Don't save automatically
  
  # Colors
  fill = venn_colors,
  alpha = 0.7,
  
  # Circle properties
  lwd = 2,
  col = "black",
  
  # Text properties
  cex = 1.2,
  fontfamily = "sans",
  fontface = "bold",
  
  # Category label properties
  cat.cex = 1.1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  
  # Margins
  margin = 0.2,
  
  # Disable log file
  disable.logging = TRUE
)

# Save the Venn diagram
png(file.path(output_dir, "pathway_venn_diagram.png"), width = 10, height = 8, 
    units = "in", res = 300, bg = "white")
grid.draw(venn_plot)
dev.off()

pdf(file.path(output_dir, "pathway_venn_diagram.pdf"), width = 10, height = 8, 
    bg = "white")
grid.draw(venn_plot)
dev.off()

cat("Saved: pathway_venn_diagram.png\n")
cat("Saved: pathway_venn_diagram.pdf\n\n")

# ------------------------------------------------------------------------------
# CREATE ENHANCED VENN WITH ANNOTATIONS
# ------------------------------------------------------------------------------

cat("Creating enhanced Venn diagram with detailed annotations...\n")

# Create enhanced version with more detailed labeling
enhanced_venn <- venn.diagram(
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
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  
  # Main title
  main = "Pathway Enrichment Overlap Between Cell Types",
  main.cex = 1.6,
  main.fontface = "bold",
  main.pos = c(0.5, 0.95),
  
  # Subtitle
  sub = paste0("Weighted Meta-Analysis Results (5%=50%, 10%=40%, 25%=10%) | Total: ", 
               nrow(meta_data), " significant pathways | FDR < 0.05"),
  sub.cex = 1.0,
  sub.pos = c(0.5, 0.05),
  
  # Margins
  margin = 0.2,
  
  # Disable log file
  disable.logging = TRUE
)

# Save enhanced version
png(file.path(output_dir, "pathway_venn_diagram_enhanced.png"), width = 12, height = 10, 
    units = "in", res = 300, bg = "white")
grid.draw(enhanced_venn)
dev.off()

pdf(file.path(output_dir, "pathway_venn_diagram_enhanced.pdf"), width = 12, height = 10, 
    bg = "white")
grid.draw(enhanced_venn)
dev.off()

cat("Saved: pathway_venn_diagram_enhanced.png\n")
cat("Saved: pathway_venn_diagram_enhanced.pdf\n\n")

# ------------------------------------------------------------------------------
# CREATE VENN BY DATABASE
# ------------------------------------------------------------------------------

cat("Creating Venn diagrams by database...\n")

databases <- unique(meta_data$database)

for(db in databases) {
  
  cat("Creating Venn for", db, "database...\n")
  
  # Filter data for this database
  db_data <- meta_data %>% filter(database == db)
  
  if(nrow(db_data) == 0) {
    cat("No data for", db, "database, skipping...\n")
    next
  }
  
  # Extract pathway lists for this database
  ca3_db <- db_data %>% filter(meta_group == "CA3_do_combined") %>% pull(ID) %>% unique()
  astro_db <- db_data %>% filter(meta_group == "376_Astro") %>% pull(ID) %>% unique()
  sst_db <- db_data %>% filter(meta_group == "78_Sst_HPF") %>% pull(ID) %>% unique()
  
  if(length(ca3_db) == 0 && length(astro_db) == 0 && length(sst_db) == 0) {
    cat("No pathways found for", db, "database, skipping...\n")
    next
  }
  
  # Create Venn for this database
  db_venn <- venn.diagram(
    x = list(
      "CA3 Dorsal" = ca3_db,
      "Astrocytes" = astro_db,
      "Sst+ HPF" = sst_db
    ),
    category.names = c(
      paste0("CA3 Dorsal\n(", length(ca3_db), ")"), 
      paste0("Astrocytes\n(", length(astro_db), ")"), 
      paste0("Sst+ HPF\n(", length(sst_db), ")")
    ),
    filename = NULL,
    
    fill = venn_colors,
    alpha = 0.6,
    lwd = 2,
    col = "black",
    cex = 1.2,
    fontfamily = "sans",
    fontface = "bold",
    cat.cex = 1.0,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    
    main = paste0("Pathway Overlap: ", db, " Database"),
    main.cex = 1.4,
    main.fontface = "bold",
    
    sub = paste0("Weighted Meta-Analysis | Total: ", nrow(db_data), " pathways"),
    sub.cex = 0.9,
    
    margin = 0.15,
    disable.logging = TRUE
  )
  
  # Save database-specific Venn
  db_filename <- gsub("[^A-Za-z0-9_]", "_", tolower(db))
  
  png(file.path(output_dir, paste0("pathway_venn_", db_filename, ".png")), 
      width = 10, height = 8, units = "in", res = 300, bg = "white")
  grid.draw(db_venn)
  dev.off()
  
  cat("Saved: pathway_venn_", db_filename, ".png\n")
}

# ------------------------------------------------------------------------------
# SAVE OVERLAP ANALYSIS
# ------------------------------------------------------------------------------

cat("\nSaving overlap analysis results...\n")

# Save overlap summary
write_csv(overlap_summary, file.path(output_dir, "pathway_overlap_summary.csv"))

# Create detailed overlap lists
overlap_details <- list(
  "CA3_Dorsal_only" = ca3_only,
  "Astrocytes_only" = astro_only,
  "Sst_HPF_only" = sst_only,
  "CA3_Dorsal_and_Astrocytes" = ca3_astro,
  "CA3_Dorsal_and_Sst_HPF" = ca3_sst,
  "Astrocytes_and_Sst_HPF" = astro_sst,
  "All_three_cell_types" = all_three
)

# Save detailed pathway lists for each overlap region
for(region_name in names(overlap_details)) {
  pathways <- overlap_details[[region_name]]
  
  if(length(pathways) > 0) {
    # Get pathway descriptions
    pathway_details <- meta_data %>%
      filter(ID %in% pathways) %>%
      select(ID, Description, database) %>%
      distinct() %>%
      arrange(database, Description)
    
    write_csv(pathway_details, file.path(output_dir, paste0(region_name, "_pathways.csv")))
  }
}

cat("Saved detailed pathway lists for each overlap region\n\n")

# ------------------------------------------------------------------------------
# FINAL SUMMARY
# ------------------------------------------------------------------------------

cat("================================================================\n")
cat("PATHWAY VENN DIAGRAM ANALYSIS COMPLETE\n")
cat("================================================================\n\n")

cat("Overlap Analysis Summary:\n")
total_unique_pathways <- length(unique(c(ca3_pathways, astro_pathways, sst_pathways)))
cat("- Total unique pathways across all cell types:", total_unique_pathways, "\n")
cat("- Shared by all three cell types:", length(all_three), "(", 
    round(length(all_three)/total_unique_pathways*100, 1), "%)\n")
cat("- Shared by exactly two cell types:", 
    length(ca3_astro) + length(ca3_sst) + length(astro_sst), "\n")
cat("- Unique to single cell types:", 
    length(ca3_only) + length(astro_only) + length(sst_only), "\n\n")

cat("Cell type-specific pathway counts:\n")
cat("- CA3 Dorsal specific:", length(ca3_only), "\n")
cat("- Astrocytes specific:", length(astro_only), "\n") 
cat("- Sst+ HPF specific:", length(sst_only), "\n\n")

cat("Files created in:", output_dir, "\n")
cat("- pathway_venn_diagram.png/pdf (simple version)\n")
cat("- pathway_venn_diagram_enhanced.png/pdf (detailed version)\n")
cat("- pathway_venn_[database].png (by database)\n")
cat("- pathway_overlap_summary.csv (overlap statistics)\n")
cat("- [region]_pathways.csv (detailed pathway lists for each overlap)\n\n")

cat("🎯 Pathway Venn diagram analysis complete!\n")
cat("📊 Visualizations show the extent of pathway sharing between cell types.\n")
cat("🔍 Use these results to understand functional overlap and specificity.\n")