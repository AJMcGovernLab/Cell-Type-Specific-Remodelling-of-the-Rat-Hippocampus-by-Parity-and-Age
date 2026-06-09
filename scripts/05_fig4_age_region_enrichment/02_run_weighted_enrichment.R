#!/usr/bin/env Rscript

# CLASSICAL EFFECTS ENRICHMENT ANALYSIS
# Test enrichment of age, region, and age×region genes in harmonized cell types

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(tidyr)

cat("==========================================================\n")
cat("CLASSICAL EFFECTS → HARMONIZED CELL TYPES ENRICHMENT\n")
cat("==========================================================\n\n")

# Auto-detect path format (Windows WSL vs Linux)
if (file.exists("\\\\wsl.localhost\\Ubuntu\\home\\ajukearth\\Parity")) {
  # Running from Windows with WSL access
  data_dir <- "\\\\wsl.localhost\\Ubuntu\\home\\ajukearth\\Parity\\Manuscript_Figures\\Improved_Integration_Analysis\\01_Data_Preparation"
  output_dir <- "\\\\wsl.localhost\\Ubuntu\\home\\ajukearth\\Parity\\Manuscript_Figures\\Improved_Integration_Analysis\\02_Classical_Effects_Analysis"
  cat("Detected Windows WSL environment\n")
} else if (file.exists("/home/ajukearth/Parity")) {
  # Running from Linux/WSL directly
  data_dir <- "/home/ajukearth/Parity/Manuscript_Figures/Improved_Integration_Analysis/01_Data_Preparation"
  output_dir <- "/home/ajukearth/Parity/Manuscript_Figures/Improved_Integration_Analysis/02_Classical_Effects_Analysis"
  cat("Detected Linux/WSL environment\n")
} else {
  stop("Cannot find Parity directory. Please check working directory.")
}

cat(sprintf("Data directory: %s\n", data_dir))
cat(sprintf("Working from: %s\n", getwd()))

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load prepared data
cat("Loading prepared data...\n")
cat("----------------------------------------\n")

# Load effect-specific gene lists
effect_genes <- readRDS(file.path(data_dir, "effect_specific_gene_lists.rds"))
cat(sprintf("  Age genes: %d\n", length(effect_genes$age)))
cat(sprintf("  Region genes: %d\n", length(effect_genes$region)))
cat(sprintf("  Age×Region genes: %d\n", length(effect_genes$age_region)))

# Load harmonized cell type basis genes
harmonized_basis <- readRDS(file.path(data_dir, "harmonized_cell_type_basis_genes.rds"))
cat(sprintf("\nHarmonized cell types with basis genes: %d\n", length(harmonized_basis)))

# Load all genes for background
all_basis_genes <- readRDS(file.path(data_dir, "all_scdc_basis_genes.rds"))
all_genes <- unique(unlist(all_basis_genes))
cat(sprintf("Background gene universe: %d genes\n", length(all_genes)))

# Debug: Check gene identifier overlap
cat("\n----------------------------------------\n")
cat("GENE IDENTIFIER DEBUGGING\n")
cat("----------------------------------------\n")

# Sample some genes from each set to check format
effect_sample <- head(effect_genes$age, 10)
basis_sample <- head(all_genes, 10)

cat("Sample effect genes (age):\n")
cat(paste(effect_sample, collapse = ", "), "\n\n")

cat("Sample basis matrix genes:\n") 
cat(paste(basis_sample, collapse = ", "), "\n\n")

# Check overlap between effect genes and basis matrix genes
age_overlap <- length(intersect(effect_genes$age, all_genes))
region_overlap <- length(intersect(effect_genes$region, all_genes))
age_region_overlap <- length(intersect(effect_genes$age_region, all_genes))

cat("Gene identifier overlap with basis matrices:\n")
cat(sprintf("  Age genes: %d/%d (%.1f%%) overlap with basis genes\n", 
            age_overlap, length(effect_genes$age), 
            100 * age_overlap / max(1, length(effect_genes$age))))
cat(sprintf("  Region genes: %d/%d (%.1f%%) overlap with basis genes\n", 
            region_overlap, length(effect_genes$region), 
            100 * region_overlap / max(1, length(effect_genes$region))))
cat(sprintf("  Age×Region genes: %d/%d (%.1f%%) overlap with basis genes\n", 
            age_region_overlap, length(effect_genes$age_region), 
            100 * age_region_overlap / max(1, length(effect_genes$age_region))))

if (age_overlap < length(effect_genes$age) * 0.5) {
  cat("\n⚠️ WARNING: Low overlap detected! Gene identifier mismatch likely.\n")
  cat("   Effect genes may use different identifiers than basis matrices.\n")
}

# Function to perform hypergeometric test
perform_enrichment_test <- function(query_genes, target_genes, background_genes) {
  # Ensure all gene sets are in the same space
  query_genes <- intersect(query_genes, background_genes)
  target_genes <- intersect(target_genes, background_genes)
  
  overlap <- intersect(query_genes, target_genes)
  overlap_count <- length(overlap)
  
  # Hypergeometric test parameters
  pop_size <- length(background_genes)  # Total genes
  success_states <- length(query_genes)  # Query genes in background
  sample_size <- length(target_genes)   # Target genes
  
  # P(X >= overlap_count)
  p_value <- phyper(overlap_count - 1, success_states, pop_size - success_states, sample_size, lower.tail = FALSE)
  
  # Calculate expected overlap and fold enrichment
  expected <- (success_states * sample_size) / pop_size
  fold_enrichment <- ifelse(expected > 0, overlap_count / expected, 0)
  
  # Calculate Jaccard index
  union_size <- length(union(query_genes, target_genes))
  jaccard_index <- ifelse(union_size > 0, overlap_count / union_size, 0)
  
  return(list(
    overlap_count = overlap_count,
    query_size = length(query_genes),
    target_size = length(target_genes),
    expected = expected,
    fold_enrichment = fold_enrichment,
    jaccard_index = jaccard_index,
    p_value = p_value,
    overlap_genes = overlap
  ))
}

# Define effects to test
effects_to_test <- list(
  age = "Age Effect",
  region = "Region Effect", 
  age_region = "Age×Region Interaction"
)

# Run enrichment analysis
cat("\n===========================================\n")
cat("ENRICHMENT ANALYSIS\n")
cat("===========================================\n")

enrichment_results <- list()

for (effect_name in names(effects_to_test)) {
  cat(sprintf("\n%s:\n", effects_to_test[[effect_name]]))
  cat("----------------------------------------", "\n")
  
  query_genes <- effect_genes[[effect_name]]
  
  if (length(query_genes) == 0) {
    cat("  No query genes available\n")
    next
  }
  
  for (cell_type in names(harmonized_basis)) {
    target_genes <- harmonized_basis[[cell_type]]
    
    result <- perform_enrichment_test(query_genes, target_genes, all_genes)
    
    # Debug output
    cat(sprintf("  %s: %d query genes, %d target genes, %d overlap, FC=%.2f, p=%.3e\n",
                cell_type, result$query_size, result$target_size, 
                result$overlap_count, result$fold_enrichment, result$p_value))
    
    # Store result
    key <- paste(effect_name, cell_type, sep = "_vs_")
    enrichment_results[[key]] <- c(
      list(effect = effect_name, cell_type = cell_type),
      result
    )
  }
}

# Convert results to data frame
cat("\nConverting results to data frame...\n")

if (length(enrichment_results) == 0) {
  cat("❌ No enrichment results to analyze - no harmonized cell types with basis genes found!\n")
  cat("   This indicates an issue with the SCDC basis gene extraction step.\n")
  cat("   Please check the harmonized cell type mapping in the previous step.\n")
  stop("No enrichment results available")
}

enrichment_df <- do.call(rbind, lapply(enrichment_results, function(x) {
  data.frame(
    effect = x$effect,
    cell_type = x$cell_type,
    overlap_count = x$overlap_count,
    query_size = x$query_size,
    target_size = x$target_size,
    expected = x$expected,
    fold_enrichment = x$fold_enrichment,
    jaccard_index = x$jaccard_index,
    p_value = x$p_value,
    stringsAsFactors = FALSE
  )
}))

# Add FDR correction
enrichment_df$fdr <- p.adjust(enrichment_df$p_value, method = "fdr")
enrichment_df$significant <- enrichment_df$fdr < 0.05

# Sort by p-value
enrichment_df <- enrichment_df %>% arrange(p_value)

# Save results
write.csv(enrichment_df, file.path(output_dir, "harmonized_enrichment_results.csv"), row.names = FALSE)

# Print summary
cat("\n===========================================\n")
cat("ENRICHMENT SUMMARY\n")
cat("===========================================\n")

# Summary by effect
for (eff in unique(enrichment_df$effect)) {
  eff_data <- enrichment_df %>% filter(effect == eff)
  n_sig <- sum(eff_data$significant)
  
  cat(sprintf("\n%s:\n", effects_to_test[[eff]]))
  cat(sprintf("  Total tests: %d\n", nrow(eff_data)))
  cat(sprintf("  Significant (FDR < 0.05): %d\n", n_sig))
  
  if (n_sig > 0) {
    top_hits <- eff_data %>% 
      filter(significant) %>%
      head(5)
    
    cat("  Top enriched cell types:\n")
    for (i in 1:nrow(top_hits)) {
      cat(sprintf("    - %s: FC=%.2f, p=%.3e\n",
                  top_hits$cell_type[i],
                  top_hits$fold_enrichment[i],
                  top_hits$p_value[i]))
    }
  }
}

# Create visualizations
cat("\n===========================================\n")
cat("CREATING VISUALIZATIONS\n")
cat("===========================================\n")

# 1. Heatmap of enrichment
cat("Creating enrichment heatmap...\n")

# Prepare matrix for heatmap
enrichment_matrix <- enrichment_df %>%
  select(effect, cell_type, fold_enrichment) %>%
  tidyr::pivot_wider(names_from = effect, values_from = fold_enrichment, values_fill = 0) %>%
  tibble::column_to_rownames("cell_type") %>%
  as.matrix()

# Create significance matrix
sig_matrix <- enrichment_df %>%
  mutate(sig_mark = ifelse(significant, "*", "")) %>%
  select(effect, cell_type, sig_mark) %>%
  tidyr::pivot_wider(names_from = effect, values_from = sig_mark, values_fill = "") %>%
  tibble::column_to_rownames("cell_type") %>%
  as.matrix()

# Create heatmap
pdf(file.path(output_dir, "harmonized_enrichment_heatmap.pdf"), width = 10, height = 8)

# Handle case where all values are identical
max_val <- max(enrichment_matrix)
min_val <- min(enrichment_matrix)

if (max_val == min_val) {
  # All values are the same, create custom breaks
  breaks <- seq(max(0, min_val - 0.1), max_val + 0.1, length.out = 101)
} else {
  breaks <- seq(min_val, max_val, length.out = 101)
}

pheatmap(enrichment_matrix,
         main = "Classical Effects Enrichment in Harmonized Cell Types\nFold Enrichment (* = FDR < 0.05)",
         color = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100),
         breaks = breaks,
         cluster_rows = nrow(enrichment_matrix) > 1,
         cluster_cols = FALSE,
         display_numbers = sig_matrix,
         fontsize_number = 12,
         fontsize_row = 10,
         fontsize_col = 11,
         angle_col = 45)

dev.off()

cat("✅ Classical effects enrichment analysis completed!\n")
cat(sprintf("Results saved in: %s\n", output_dir))