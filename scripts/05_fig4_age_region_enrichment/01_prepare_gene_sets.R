#!/usr/bin/env Rscript

# SCRIPT 2: PREPARE GENE SETS FOR FUNCTIONAL ENRICHMENT ANALYSIS
# This script converts gene IDs and prepares gene sets for enrichment analysis

library(dplyr)
library(readr)

cat("================================================================\n")
cat("PREPARING GENE SETS FOR FUNCTIONAL ENRICHMENT ANALYSIS\n")
cat("================================================================\n\n")

# Auto-detect path format (Windows WSL vs Linux)
if (file.exists("\\\\wsl.localhost\\Ubuntu\\home\\ajukearth\\Parity")) {
  base_dir <- "\\\\wsl.localhost\\Ubuntu\\home\\ajukearth\\Parity"
  cat("Detected Windows WSL environment\n")
} else if (file.exists("/home/ajukearth/Parity")) {
  base_dir <- "/home/ajukearth/Parity"
  cat("Detected Linux/WSL environment\n")
} else {
  stop("Cannot find Parity directory. Please check working directory.")
}

output_dir <- file.path(base_dir, "Manuscript_Figures/Improved_Integration_Analysis/05_Comprehensive_Functional_Enrichment")
gene_mapping_dir <- file.path(base_dir, "Manuscript_Figures/Full integration/01_Data_Preparation/data")

cat(sprintf("Working from: %s\n", getwd()))
cat(sprintf("Output dir: %s\n", output_dir))

cat("\n===========================================\n")
cat("STEP 1: LOAD OVERLAP ANALYSIS RESULTS\n")
cat("===========================================\n")

# Load comprehensive overlap results
overlap_details_file <- file.path(output_dir, "results/summary_tables/all_effect_overlaps_details.csv")
overlap_summary_file <- file.path(output_dir, "results/summary_tables/all_effect_overlaps_summary.csv")

overlap_details <- read.csv(overlap_details_file, stringsAsFactors = FALSE)
overlap_summary <- read.csv(overlap_summary_file, stringsAsFactors = FALSE)

cat(sprintf("✓ Loaded overlap details: %d genes across all effects\n", nrow(overlap_details)))
cat(sprintf("✓ Loaded overlap summary: %d cell type-effect combinations\n", nrow(overlap_summary)))

# Summary by effect type
effect_summary <- overlap_details %>%
  group_by(effect_type) %>%
  summarise(
    n_genes = n(),
    n_cell_types = n_distinct(cell_type),
    .groups = "drop"
  )

print(effect_summary)

cat("\n===========================================\n")
cat("STEP 2: LOAD GENE ID MAPPING\n")
cat("===========================================\n")

# Load gene ID mapping file (ENSRNOG to gene symbols)
gene_mapping_file <- file.path(gene_mapping_dir, "gene_id_mapping.csv")

if (file.exists(gene_mapping_file)) {
  gene_mapping <- read.csv(gene_mapping_file, stringsAsFactors = FALSE)
  cat(sprintf("✓ Loaded gene mapping: %d mappings\n", nrow(gene_mapping)))
  
  # Check mapping columns
  cat("Gene mapping columns:", paste(colnames(gene_mapping), collapse = ", "), "\n")
  
  # Ensure we have the right column names
  if ("ensembl_gene_id" %in% colnames(gene_mapping)) {
    mapping_id_col <- "ensembl_gene_id"
  } else if ("gene_id" %in% colnames(gene_mapping)) {
    mapping_id_col <- "gene_id"
  } else {
    cat("⚠ Could not find gene ID column in mapping file\n")
    mapping_id_col <- colnames(gene_mapping)[1]
  }
  
  if ("external_gene_name" %in% colnames(gene_mapping)) {
    mapping_symbol_col <- "external_gene_name"
  } else if ("gene_name" %in% colnames(gene_mapping)) {
    mapping_symbol_col <- "gene_name"
  } else {
    cat("⚠ Could not find gene symbol column in mapping file\n")
    mapping_symbol_col <- colnames(gene_mapping)[2]
  }
  
  cat(sprintf("Using ID column: %s, Symbol column: %s\n", mapping_id_col, mapping_symbol_col))
  
} else {
  cat("⚠ Gene mapping file not found, will use gene IDs as-is\n")
  gene_mapping <- NULL
}

cat("\n===========================================\n")
cat("STEP 3: GENE ID CONVERSION FUNCTION\n")
cat("===========================================\n")

# Function to convert ENSRNOG IDs to gene symbols
convert_gene_ids <- function(gene_ids, mapping_data = gene_mapping) {
  if (is.null(mapping_data)) {
    return(data.frame(
      ensrnog_id = gene_ids,
      gene_symbol = gene_ids,
      stringsAsFactors = FALSE
    ))
  }
  
  # Create conversion data frame
  conversion_df <- data.frame(
    ensrnog_id = gene_ids,
    stringsAsFactors = FALSE
  )
  
  # Map to gene symbols
  conversion_df <- merge(conversion_df, 
                        mapping_data[, c(mapping_id_col, mapping_symbol_col)], 
                        by.x = "ensrnog_id", 
                        by.y = mapping_id_col, 
                        all.x = TRUE)
  
  # Use gene symbol if available, otherwise use ENSRNOG ID
  conversion_df$gene_symbol <- ifelse(
    is.na(conversion_df[[mapping_symbol_col]]) | conversion_df[[mapping_symbol_col]] == "",
    conversion_df$ensrnog_id,
    conversion_df[[mapping_symbol_col]]
  )
  
  # Clean up column names
  conversion_df <- conversion_df[, c("ensrnog_id", "gene_symbol")]
  
  return(conversion_df)
}

cat("✓ Gene ID conversion function created\n")

cat("\n===========================================\n")
cat("STEP 4: PREPARE GENE SETS BY EFFECT AND CELL TYPE\n")
cat("===========================================\n")

# Process each effect type and cell type combination
gene_sets_summary <- data.frame()

for (effect in c("age", "region", "age_region")) {
  cat(sprintf("\nProcessing %s effects...\n", toupper(effect)))
  
  effect_data <- overlap_details[overlap_details$effect_type == effect, ]
  
  if (nrow(effect_data) > 0) {
    for (cell_type in unique(effect_data$cell_type)) {
      cell_genes <- effect_data[effect_data$cell_type == cell_type, ]
      
      if (nrow(cell_genes) > 0) {
        # Convert gene IDs
        gene_conversion <- convert_gene_ids(cell_genes$gene_id)
        
        # Merge with effect data
        cell_genes_mapped <- merge(cell_genes, gene_conversion, 
                                  by.x = "gene_id", by.y = "ensrnog_id", all.x = TRUE)
        
        # Analyze directions (up/down regulation)
        if ("age_log2fc" %in% colnames(cell_genes_mapped)) {
          upregulated_genes <- cell_genes_mapped[!is.na(cell_genes_mapped$age_log2fc) & 
                                               cell_genes_mapped$age_log2fc > 0, ]
          downregulated_genes <- cell_genes_mapped[!is.na(cell_genes_mapped$age_log2fc) & 
                                                 cell_genes_mapped$age_log2fc < 0, ]
        } else {
          upregulated_genes <- data.frame()
          downregulated_genes <- data.frame()
        }
        
        # Clean file names for saving
        clean_cell_name <- gsub("[^A-Za-z0-9_]", "_", cell_type)
        
        # Save complete gene set
        complete_file <- file.path(output_dir, "data/gene_mappings", 
                                  paste0(effect, "_", clean_cell_name, "_complete_genes.csv"))
        write.csv(cell_genes_mapped, complete_file, row.names = FALSE)
        
        # Save gene symbols for enrichment
        symbols_file <- file.path(output_dir, "data/gene_mappings", 
                                 paste0(effect, "_", clean_cell_name, "_gene_symbols.txt"))
        write.table(cell_genes_mapped$gene_symbol, symbols_file, 
                   row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        # Save ENSRNOG IDs for enrichment
        ensrnog_file <- file.path(output_dir, "data/gene_mappings", 
                                 paste0(effect, "_", clean_cell_name, "_ensrnog_ids.txt"))
        write.table(cell_genes_mapped$gene_id, ensrnog_file, 
                   row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        # Save direction-specific gene sets if available
        if (nrow(upregulated_genes) > 0) {
          up_symbols_file <- file.path(output_dir, "data/gene_mappings", 
                                      paste0(effect, "_", clean_cell_name, "_upregulated_symbols.txt"))
          write.table(upregulated_genes$gene_symbol, up_symbols_file, 
                     row.names = FALSE, col.names = FALSE, quote = FALSE)
        }
        
        if (nrow(downregulated_genes) > 0) {
          down_symbols_file <- file.path(output_dir, "data/gene_mappings", 
                                        paste0(effect, "_", clean_cell_name, "_downregulated_symbols.txt"))
          write.table(downregulated_genes$gene_symbol, down_symbols_file, 
                     row.names = FALSE, col.names = FALSE, quote = FALSE)
        }
        
        # Record summary
        summary_row <- data.frame(
          effect_type = effect,
          cell_type = cell_type,
          total_genes = nrow(cell_genes_mapped),
          upregulated_genes = nrow(upregulated_genes),
          downregulated_genes = nrow(downregulated_genes),
          mapped_symbols = sum(!is.na(cell_genes_mapped$gene_symbol) & 
                              cell_genes_mapped$gene_symbol != cell_genes_mapped$gene_id),
          mapping_rate = round(100 * sum(!is.na(cell_genes_mapped$gene_symbol) & 
                                        cell_genes_mapped$gene_symbol != cell_genes_mapped$gene_id) / 
                              nrow(cell_genes_mapped), 1),
          stringsAsFactors = FALSE
        )
        
        gene_sets_summary <- rbind(gene_sets_summary, summary_row)
        
        cat(sprintf("  %s: %d total genes (%d up, %d down), %.1f%% mapped to symbols\n",
                   cell_type, nrow(cell_genes_mapped), nrow(upregulated_genes), 
                   nrow(downregulated_genes), summary_row$mapping_rate))
      }
    }
  }
}

cat("\n===========================================\n")
cat("STEP 5: CREATE COMBINED GENE SETS\n")
cat("===========================================\n")

# Create combined gene sets across cell types for each effect
for (effect in c("age", "region", "age_region")) {
  effect_data <- overlap_details[overlap_details$effect_type == effect, ]
  
  if (nrow(effect_data) > 0) {
    # Get all unique genes for this effect
    all_effect_genes <- unique(effect_data$gene_id)
    
    # Convert to symbols
    effect_conversion <- convert_gene_ids(all_effect_genes)
    
    # Save combined gene sets
    combined_symbols_file <- file.path(output_dir, "data/gene_mappings", 
                                      paste0(effect, "_all_cell_types_gene_symbols.txt"))
    write.table(effect_conversion$gene_symbol, combined_symbols_file, 
               row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    combined_ensrnog_file <- file.path(output_dir, "data/gene_mappings", 
                                      paste0(effect, "_all_cell_types_ensrnog_ids.txt"))
    write.table(effect_conversion$ensrnog_id, combined_ensrnog_file, 
               row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    cat(sprintf("✓ Created combined %s gene set: %d genes\n", effect, length(all_effect_genes)))
  }
}

cat("\n===========================================\n")
cat("STEP 6: PREPARE BACKGROUND GENE SETS\n")
cat("===========================================\n")

# Create background gene sets from all tested genes
age_genes_file <- file.path(base_dir, "Manuscript_Figures/Improved_Integration_Analysis/01_Data_Preparation/age_effect_genes.csv")
region_genes_file <- file.path(base_dir, "Manuscript_Figures/Improved_Integration_Analysis/01_Data_Preparation/region_effect_genes.csv")
interaction_genes_file <- file.path(base_dir, "Manuscript_Figures/Improved_Integration_Analysis/01_Data_Preparation/age_region_interaction_genes.csv")

all_tested_genes <- c()

if (file.exists(age_genes_file)) {
  age_genes <- read.csv(age_genes_file, stringsAsFactors = FALSE)
  all_tested_genes <- c(all_tested_genes, age_genes$gene_id)
  cat(sprintf("✓ Added %d age effect genes to background\n", nrow(age_genes)))
}

if (file.exists(region_genes_file)) {
  region_genes <- read.csv(region_genes_file, stringsAsFactors = FALSE)
  all_tested_genes <- c(all_tested_genes, region_genes$gene_id)
  cat(sprintf("✓ Added %d region effect genes to background\n", nrow(region_genes)))
}

if (file.exists(interaction_genes_file)) {
  interaction_genes <- read.csv(interaction_genes_file, stringsAsFactors = FALSE)
  all_tested_genes <- c(all_tested_genes, interaction_genes$gene_id)
  cat(sprintf("✓ Added %d interaction effect genes to background\n", nrow(interaction_genes)))
}

# Create unique background set
background_genes <- unique(all_tested_genes)
background_conversion <- convert_gene_ids(background_genes)

# Save background gene sets
background_symbols_file <- file.path(output_dir, "data/gene_mappings", "background_gene_symbols.txt")
write.table(background_conversion$gene_symbol, background_symbols_file, 
           row.names = FALSE, col.names = FALSE, quote = FALSE)

background_ensrnog_file <- file.path(output_dir, "data/gene_mappings", "background_ensrnog_ids.txt")
write.table(background_conversion$ensrnog_id, background_ensrnog_file, 
           row.names = FALSE, col.names = FALSE, quote = FALSE)

cat(sprintf("✓ Created background gene set: %d unique genes\n", length(background_genes)))

cat("\n===========================================\n")
cat("STEP 7: SAVE SUMMARY RESULTS\n")
cat("===========================================\n")

# Save gene sets summary
summary_file <- file.path(output_dir, "results/summary_tables/gene_sets_preparation_summary.csv")
write.csv(gene_sets_summary, summary_file, row.names = FALSE)

# Create overall summary
overall_summary <- data.frame(
  category = c("Total overlap genes", "Age effect genes", "Region effect genes", 
              "Age×Region interaction genes", "Background genes", "Gene sets created"),
  count = c(nrow(overlap_details), 
           sum(overlap_details$effect_type == "age"),
           sum(overlap_details$effect_type == "region"),
           sum(overlap_details$effect_type == "age_region"),
           length(background_genes),
           nrow(gene_sets_summary)),
  stringsAsFactors = FALSE
)

overall_summary_file <- file.path(output_dir, "results/summary_tables/gene_preparation_overall_summary.csv")
write.csv(overall_summary, overall_summary_file, row.names = FALSE)

cat("\n===========================================\n")
cat("FINAL SUMMARY\n")
cat("===========================================\n")

print(overall_summary)

cat("\nGene Sets Created by Effect Type:\n")
cat("---------------------------------\n")
effect_sets_summary <- gene_sets_summary %>%
  group_by(effect_type) %>%
  summarise(
    n_cell_types = n(),
    total_genes = sum(total_genes),
    avg_mapping_rate = round(mean(mapping_rate), 1),
    .groups = "drop"
  )

print(effect_sets_summary)

cat("\nTop Cell Types by Gene Count:\n")
cat("-----------------------------\n")
top_gene_sets <- gene_sets_summary %>%
  arrange(desc(total_genes)) %>%
  select(effect_type, cell_type, total_genes, upregulated_genes, downregulated_genes, mapping_rate) %>%
  head(10)

print(top_gene_sets)

cat("\n✓ Gene set preparation completed!\n")
cat(sprintf("Gene sets saved in: %s/data/gene_mappings/\n", output_dir))
cat(sprintf("Summary saved in: %s/results/summary_tables/\n", output_dir))

# Count files created
gene_files <- list.files(file.path(output_dir, "data/gene_mappings"), pattern = "\\.(txt|csv)$")
cat(sprintf("Files created: %d gene set files\n", length(gene_files)))

cat("\nReady for functional enrichment analysis!\n")