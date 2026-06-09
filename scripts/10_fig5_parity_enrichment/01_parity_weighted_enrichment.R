#!/usr/bin/env Rscript

# PARITY-SPECIFIC WEIGHTED ENRICHMENT ANALYSIS - FOUR METHODS
# Adapted from Paper_Level/6 framework for parity effects analysis
# Uses 4 different weighting strategies on parity genes and cell types

# Windows-safe package loading
safe_library <- function(package) {
  if(!require(package, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Warning: Package %s not available. Some analyses may not work.\n", package))
    return(FALSE)
  }
  return(TRUE)
}

# Load required packages
cat("Loading packages for parity weighted enrichment...\n")
packages_loaded <- list(
  tidyverse = safe_library("tidyverse"),
  clusterProfiler = safe_library("clusterProfiler"),
  org.Mm.eg.db = safe_library("org.Mm.eg.db"),
  DOSE = safe_library("DOSE"),
  ReactomePA = safe_library("ReactomePA"),
  msigdbr = safe_library("msigdbr"),
  enrichplot = safe_library("enrichplot")
)

cat("================================================================\n")
cat("PARITY-SPECIFIC WEIGHTED ENRICHMENT ANALYSIS\n")
cat("================================================================\n\n")

# Paths via curated config. The cell_matrices inputs are not staged in the
# Repository checkpoints (they live in Set 1 working directories). When this
# script's inputs are unavailable the curated pipeline falls back to the
# staged checkpoint at checkpoints/enrichment_parity/all_parity_weighted_enrichment_results.csv,
# which 02_meta_analysis_integration.R reads directly. This script is
# preserved for full upstream reproducibility but short-circuits if inputs
# are missing.
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))
cell_matrices_dir <- file.path(CHECKPOINT_DIR, "enrichment_parity", "cell_matrices")
output_dir <- file.path(OUT$fig5ce, "weighted_enrichment_intermediate")
if (!dir.exists(cell_matrices_dir)) {
  message(sprintf("[skip] %s missing; using staged checkpoint instead.",
                  cell_matrices_dir))
  message("[skip] 02_meta_analysis_integration.R will read all_parity_weighted_enrichment_results.csv directly.")
  if (sys.nframe() == 0L) quit(save = "no", status = 0)
}

cat("Directory setup:\n")
cat(sprintf("  Base directory: %s\n", base_dir))
cat(sprintf("  Cell matrices directory: %s\n", cell_matrices_dir))
cat(sprintf("  Output directory: %s\n", output_dir))

# Create output directories
dir.create(file.path(output_dir, "results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. LOAD PARITY DATA
# ------------------------------------------------------------------------------
cat("\nLoading parity association data...\n")

# Load parity gene-cell associations
parity_associations <- read.csv(file.path(cell_matrices_dir, "parity_gene_cell_associations.csv"), 
                               stringsAsFactors = FALSE)
cat(sprintf("  Parity associations: %d gene-cell combinations\n", nrow(parity_associations)))

# Define the 4 parity-affected cell types
parity_cell_types <- c("358_CA3-do", "376_Astro", "78_Sst_HPF", "356_CA3-do")

# Load full expression matrices for composite weighting
expression_matrices <- list()
for(cell_type in parity_cell_types) {
  cell_clean <- gsub("[^A-Za-z0-9]", "_", cell_type)
  expr_file <- file.path(cell_matrices_dir, cell_clean, "full_expression_matrix.csv")
  
  if(file.exists(expr_file)) {
    expr_data <- read.csv(expr_file, stringsAsFactors = FALSE)
    expression_matrices[[cell_type]] <- expr_data
    cat(sprintf("  %s: %d genes loaded\n", cell_type, nrow(expr_data)))
  } else {
    cat(sprintf("  Warning: Expression matrix not found for %s\n", cell_type))
  }
}

# ------------------------------------------------------------------------------
# 2. PARITY WEIGHTING STRATEGY FUNCTIONS
# ------------------------------------------------------------------------------
cat("\nDefining parity-specific weighting strategies...\n")

# Function to create weighted gene lists for parity analysis
create_parity_weighted_gene_list <- function(cell_type, method = "association") {
  
  # Get parity associations for this cell type
  cell_data <- parity_associations[parity_associations$cell_type == cell_type, ]
  
  if(nrow(cell_data) == 0) {
    cat(sprintf("  Warning: No parity data found for %s\n", cell_type))
    return(NULL)
  }
  
  # Get full expression matrix for composite weighting
  expr_matrix <- expression_matrices[[cell_type]]
  
  if(is.null(expr_matrix)) {
    cat(sprintf("  Warning: No expression matrix found for %s\n", cell_type))
    return(NULL)
  }
  
  # Merge parity genes with expression data
  merged_data <- cell_data %>%
    left_join(expr_matrix, by = c("gene_symbol" = "gene")) %>%
    filter(!is.na(composite_score))  # Keep only genes with expression data
  
  if(nrow(merged_data) == 0) {
    cat(sprintf("  Warning: No merged data for %s\n", cell_type))
    return(NULL)
  }
  
  # Calculate weights based on method
  weights <- switch(method,
    "association" = {
      # Pure parity association score weighting (0.6-1.0 based on method count)
      merged_data$association_score
    },
    
    "expression" = {
      # Association × parity fold change magnitude
      if("abs_log2fc" %in% colnames(merged_data) && any(!is.na(merged_data$abs_log2fc))) {
        # Use parity-specific fold change
        norm_fc <- merged_data$abs_log2fc / max(merged_data$abs_log2fc, na.rm = TRUE)
        merged_data$association_score * (0.6 + 0.4 * norm_fc)
      } else {
        merged_data$association_score
      }
    },
    
    "statistical" = {
      # Association × parity statistical significance
      if("fdr" %in% colnames(merged_data) && any(!is.na(merged_data$fdr))) {
        neg_log_fdr <- -log10(merged_data$fdr + 1e-100)
        norm_stat <- neg_log_fdr / max(neg_log_fdr, na.rm = TRUE)
        merged_data$association_score * (0.6 + 0.4 * norm_stat)
      } else {
        merged_data$association_score
      }
    },
    
    "composite" = {
      # Parity association + cell-type specificity + expression magnitude
      base_weight <- merged_data$association_score * 0.4  # 40% parity association
      
      # 40% cell-type specificity from expression matrix
      cell_specificity <- if("composite_score" %in% colnames(merged_data)) {
        0.4 * merged_data$composite_score
      } else { 0 }
      
      # 20% parity magnitude
      expr_component <- if("abs_log2fc" %in% colnames(merged_data) && any(!is.na(merged_data$abs_log2fc))) {
        0.2 * (merged_data$abs_log2fc / max(merged_data$abs_log2fc, na.rm = TRUE))
      } else { 0 }
      
      # Final composite weight
      base_weight + cell_specificity + expr_component
    }
  )
  
  # Add direction from parity fold change
  if("log2fc" %in% colnames(merged_data) && any(!is.na(merged_data$log2fc))) {
    direction <- sign(merged_data$log2fc)
    direction[is.na(direction)] <- 1  # Default to positive
    weights <- weights * direction
  }
  
  # Create named vector and sort by weight (descending)
  gene_list <- setNames(weights, merged_data$gene_symbol)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}

# Test weighting strategies
cat("Testing parity weighting strategies...\n")
for(cell_type in parity_cell_types) {
  cat(sprintf("  Testing %s:\n", cell_type))
  for(method in c("association", "expression", "statistical", "composite")) {
    test_list <- create_parity_weighted_gene_list(cell_type, method)
    if(!is.null(test_list)) {
      cat(sprintf("    %s: %d genes, range [%.3f, %.3f]\n", 
                  method, length(test_list), min(test_list, na.rm = TRUE), max(test_list, na.rm = TRUE)))
    }
  }
}

# ------------------------------------------------------------------------------
# 3. PATHWAY DATABASE SETUP
# ------------------------------------------------------------------------------
cat("\nSetting up pathway databases for parity analysis...\n")

# Initialize pathway databases
pathway_databases <- list()

# GO databases
if(packages_loaded$org.Mm.eg.db) {
  pathway_databases$GO_BP <- list(type = "GO", ont = "BP", name = "GO Biological Process")
  pathway_databases$GO_MF <- list(type = "GO", ont = "MF", name = "GO Molecular Function")
  pathway_databases$GO_CC <- list(type = "GO", ont = "CC", name = "GO Cellular Component")
  cat("  ✓ GO databases configured\n")
}

# KEGG database
if(packages_loaded$clusterProfiler) {
  pathway_databases$KEGG <- list(type = "KEGG", name = "KEGG Pathways")
  cat("  ✓ KEGG database configured\n")
}

# Reactome database
if(packages_loaded$ReactomePA) {
  pathway_databases$Reactome <- list(type = "Reactome", name = "Reactome Pathways")
  cat("  ✓ Reactome database configured\n")
}

# MSigDB databases
if(packages_loaded$msigdbr) {
  pathway_databases$MSigDB_H <- list(type = "MSigDB", collection = "H", name = "MSigDB Hallmark")
  pathway_databases$MSigDB_C2 <- list(type = "MSigDB", collection = "C2", name = "MSigDB Curated")
  cat("  ✓ MSigDB databases configured\n")
}

cat(sprintf("Total databases configured: %d\n", length(pathway_databases)))

# ------------------------------------------------------------------------------
# 4. ENRICHMENT ANALYSIS FUNCTION
# ------------------------------------------------------------------------------
cat("\nDefining enrichment analysis functions...\n")

# Function to run enrichment for a specific method and database
run_parity_enrichment <- function(gene_list, database_info, cell_type, method) {
  
  if(is.null(gene_list) || length(gene_list) == 0) {
    return(NULL)
  }
  
  tryCatch({
    result <- switch(database_info$type,
      "GO" = {
        if(packages_loaded$clusterProfiler) {
          gseGO(geneList = gene_list,
                OrgDb = org.Mm.eg.db,
                ont = database_info$ont,
                keyType = "SYMBOL",
                minGSSize = 15,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                verbose = FALSE)
        } else { NULL }
      },
      
      "KEGG" = {
        if(packages_loaded$clusterProfiler) {
          # Convert symbols to ENTREZ IDs
          entrez_genes <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", 
                              OrgDb = org.Mm.eg.db, drop = TRUE)
          if(nrow(entrez_genes) > 0) {
            entrez_list <- gene_list[entrez_genes$SYMBOL]
            names(entrez_list) <- entrez_genes$ENTREZID
            gseKEGG(geneList = sort(entrez_list, decreasing = TRUE),
                    organism = "mmu",
                    minGSSize = 15,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)
          } else { NULL }
        } else { NULL }
      },
      
      "Reactome" = {
        if(packages_loaded$ReactomePA) {
          # Convert symbols to ENTREZ IDs
          entrez_genes <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", 
                              OrgDb = org.Mm.eg.db, drop = TRUE)
          if(nrow(entrez_genes) > 0) {
            entrez_list <- gene_list[entrez_genes$SYMBOL]
            names(entrez_list) <- entrez_genes$ENTREZID
            gsePathway(geneList = sort(entrez_list, decreasing = TRUE),
                      organism = "mouse",
                      minGSSize = 15,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      verbose = FALSE)
          } else { NULL }
        } else { NULL }
      },
      
      "MSigDB" = {
        if(packages_loaded$msigdbr) {
          # Get MSigDB gene sets
          msig_sets <- msigdbr(species = "Mus musculus", category = database_info$collection)
          msig_list <- split(msig_sets$gene_symbol, msig_sets$gs_name)
          
          GSEA(geneList = gene_list,
               TERM2GENE = msig_sets[,c("gs_name", "gene_symbol")],
               minGSSize = 15,
               maxGSSize = 500,
               pvalueCutoff = 0.05,
               verbose = FALSE)
        } else { NULL }
      }
    )
    
    if(!is.null(result) && nrow(result@result) > 0) {
      # Add metadata
      result@result$cell_type <- cell_type
      result@result$method <- method
      result@result$database <- database_info$name
      result@result$effect <- "parity"
      
      return(result@result)
    }
    
  }, error = function(e) {
    cat(sprintf("    Error in %s %s: %s\n", database_info$name, method, e$message))
    return(NULL)
  })
  
  return(NULL)
}

# ------------------------------------------------------------------------------
# 5. RUN COMPREHENSIVE ENRICHMENT ANALYSIS
# ------------------------------------------------------------------------------
cat("\nRunning comprehensive parity enrichment analysis...\n")

# Initialize results storage
all_enrichment_results <- list()
method_names <- c("association", "expression", "statistical", "composite")

# Run enrichment for each combination
total_analyses <- length(parity_cell_types) * length(method_names) * length(pathway_databases)
current_analysis <- 0

for(cell_type in parity_cell_types) {
  cat(sprintf("\nProcessing cell type: %s\n", cell_type))
  
  for(method in method_names) {
    cat(sprintf("  Method: %s\n", method))
    
    # Create weighted gene list for this method
    gene_list <- create_parity_weighted_gene_list(cell_type, method)
    
    if(is.null(gene_list)) {
      cat(sprintf("    Skipping %s - no gene list\n", method))
      next
    }
    
    for(db_name in names(pathway_databases)) {
      current_analysis <- current_analysis + 1
      cat(sprintf("    Database: %s (%d/%d)\n", db_name, current_analysis, total_analyses))
      
      # Run enrichment
      enrichment_result <- run_parity_enrichment(gene_list, pathway_databases[[db_name]], 
                                                 cell_type, method)
      
      if(!is.null(enrichment_result)) {
        result_key <- paste(cell_type, method, db_name, sep = "_")
        all_enrichment_results[[result_key]] <- enrichment_result
        cat(sprintf("      ✓ %d pathways found\n", nrow(enrichment_result)))
      } else {
        cat(sprintf("      - No significant pathways\n"))
      }
    }
  }
}

# ------------------------------------------------------------------------------
# 6. SAVE RESULTS
# ------------------------------------------------------------------------------
cat("\nSaving enrichment results...\n")

# Combine all results
if(length(all_enrichment_results) > 0) {
  combined_results <- do.call(rbind, all_enrichment_results)
  
  # Save by method
  for(method in method_names) {
    method_results <- combined_results[combined_results$method == method, ]
    if(nrow(method_results) > 0) {
      write.csv(method_results, 
                file.path(output_dir, "results", paste0(method, "_weighted_enrichment.csv")), 
                row.names = FALSE)
      cat(sprintf("  ✓ %s: %d pathways saved\n", method, nrow(method_results)))
    }
  }
  
  # Save complete results
  write.csv(combined_results, 
            file.path(output_dir, "results", "all_parity_weighted_enrichment_results.csv"), 
            row.names = FALSE)
  
  cat(sprintf("✓ Complete results: %d total pathway enrichments saved\n", nrow(combined_results)))
  
} else {
  cat("Warning: No enrichment results to save\n")
}

cat("\n================================================================\n")
cat("PARITY WEIGHTED ENRICHMENT ANALYSIS COMPLETE\n")
cat("================================================================\n")

# Summary statistics
if(exists("combined_results") && nrow(combined_results) > 0) {
  summary_stats <- combined_results %>%
    group_by(cell_type, method, database) %>%
    summarise(
      pathways = n(),
      significant = sum(p.adjust < 0.05, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\n📊 SUMMARY STATISTICS:\n")
  print(summary_stats)
  
  write.csv(summary_stats, 
            file.path(output_dir, "results", "parity_enrichment_summary.csv"), 
            row.names = FALSE)
  
  cat("\n🎯 FILES SAVED:\n")
  cat("  • Individual method results: weighted_enrichment/results/[method]_weighted_enrichment.csv\n")
  cat("  • Complete results: all_parity_weighted_enrichment_results.csv\n")
  cat("  • Summary statistics: parity_enrichment_summary.csv\n")
  
  cat("\n✅ Ready for Step 5: Meta-analysis integration\n")
}