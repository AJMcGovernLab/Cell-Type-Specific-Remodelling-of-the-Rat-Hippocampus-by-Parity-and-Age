# 03_pc_driver_analysis.R — gene-removal driver analysis using the paper's
# pre-computed PCA loadings + PC coordinates.
#
# The paper's `improved_comprehensive_driver_analysis.R` reads three inputs
# produced by upstream scripts (pca_analysis_clean.R + create_binary_parity_pca.R
# + create_significance_outlined_heatmap.R):
#   - pca_loadings_top10.csv
#   - binary_parity_pca_coordinates.csv
#   - factor_association_significant_models_summary.csv
#   - clean_normalized_counts.csv  (expression matrix with Gene_1..Gene_N row names)
# All four are now staged in Repository/checkpoints/pca_drivers/ and
# Repository/data/bulk_rnaseq/.

suppressPackageStartupMessages({ library(dplyr); library(data.table) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

DRIVER_THRESH <- 0.1   # paper's significant_threshold
PCA_CKPT_DIR  <- file.path(CHECKPOINT_DIR, "pca_drivers")

main <- function() {
  out_dir <- file.path(OUT$suppfig234, "PC_Drivers")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Expression matrix — paper's script uses clean_normalized_counts.csv with
  # gene_name = "Gene_1", "Gene_2", ..., "Gene_N" (just row indices).
  expr_df <- read.csv(file.path(BULK_DIR, "clean_normalized_counts.csv"),
                      check.names = FALSE)
  expr_mat <- as.matrix(expr_df)
  rownames(expr_mat) <- paste0("Gene_", seq_len(nrow(expr_mat)))

  loadings  <- read.csv(file.path(PCA_CKPT_DIR, "pca_loadings_top10.csv"),
                        stringsAsFactors = FALSE)
  pc_coords <- read.csv(file.path(PCA_CKPT_DIR, "binary_parity_pca_coordinates.csv"),
                        stringsAsFactors = FALSE) |>
    mutate(
      age           = as.factor(age),
      parity_binary = as.factor(parity_binary),
      region        = as.factor(region)
    )
  sig_models <- read.csv(file.path(PCA_CKPT_DIR,
                                    "factor_association_significant_models_summary.csv"),
                         stringsAsFactors = FALSE)

  # Only the Parity_Binary associations for Supp Fig 4.
  sig_models <- sig_models[sig_models$Model == "Parity_Binary", ]
  message(sprintf("[pc] parity-associated PCs: %s",
                  paste(sig_models$PC, collapse = ", ")))

  gene_meta <- fread(GENE_META, data.table = FALSE)

  # Build expr_to_ensembl_id map: paper's loadings use Gene_N, we need the Ensembl.
  # Gene_N corresponds to row N in the original counts CSV (1-indexed).
  gene_to_ensembl <- if ("genes" %in% colnames(gene_meta)) gene_meta$genes[seq_len(nrow(expr_mat))]
                     else rownames(expr_mat)

  for (i in seq_len(nrow(sig_models))) {
    pc_name   <- sig_models$PC[i]
    pc_num    <- as.numeric(sub("PC", "", pc_name))

    pc_col    <- pc_name
    pc_scores <- pc_coords[[pc_col]]

    loading_vec <- loadings[[pc_col]]
    names(loading_vec) <- loadings$gene

    message(sprintf("[pc%d] fitting base model (Parity_Binary)", pc_num))
    formula_full <- as.formula(paste(pc_col, "~ parity_binary"))
    base_fit  <- lm(formula_full, data = pc_coords)
    base_sum  <- summary(base_fit)
    base_r2   <- base_sum$r.squared
    # Extract coefficient row for parity_binary
    term_row <- grep("^parity_binary", rownames(base_sum$coefficients), value = TRUE)[1]
    base_p   <- base_sum$coefficients[term_row, "Pr(>|t|)"]
    base_coef <- base_sum$coefficients[term_row, "Estimate"]

    # Test each gene
    n_genes <- length(loading_vec)
    reduced_r2 <- numeric(n_genes)
    reduced_p  <- numeric(n_genes)
    reduced_coef <- numeric(n_genes)

    for (g in seq_len(n_genes)) {
      gene_name <- names(loading_vec)[g]
      gene_load <- loading_vec[g]

      # Row index from Gene_<N>
      row_idx <- as.integer(sub("Gene_", "", gene_name))
      if (is.na(row_idx) || row_idx > nrow(expr_mat)) next
      gene_expr <- expr_mat[row_idx, ]

      contrib  <- gene_expr * gene_load
      # Map column names to sample rows in pc_coords
      reduced_scores <- pc_scores - contrib[match(pc_coords$sample %||%
                                                    rownames(pc_coords), colnames(expr_mat))]
      tmp <- pc_coords; tmp[[pc_col]] <- reduced_scores
      fit <- tryCatch(lm(formula_full, data = tmp), error = function(e) NULL)
      if (is.null(fit)) { reduced_p[g] <- 1; reduced_r2[g] <- 0; next }
      s <- summary(fit)
      reduced_r2[g]   <- s$r.squared
      if (term_row %in% rownames(s$coefficients)) {
        reduced_p[g]    <- s$coefficients[term_row, "Pr(>|t|)"]
        reduced_coef[g] <- s$coefficients[term_row, "Estimate"]
      } else {
        reduced_p[g]    <- 1.0
        reduced_coef[g] <- 0
      }

      if (g %% 2000 == 0) message(sprintf("[pc%d]   %d/%d", pc_num, g, n_genes))
    }

    log_p_deg  <- log10(pmax(reduced_p, 1e-300) / pmax(base_p, 1e-300))
    r2_loss    <- base_r2 - reduced_r2
    contribution_score <- log_p_deg + r2_loss * 100

    df <- data.frame(
      ensembl_id         = gene_to_ensembl[seq_len(n_genes)],
      gene_name          = names(loading_vec),
      loading            = loading_vec,
      abs_loading        = abs(loading_vec),
      original_pval      = base_p,
      reduced_pval       = reduced_p,
      log_pval_degradation = log_p_deg,
      original_coef      = base_coef,
      reduced_coef       = reduced_coef,
      coef_degradation   = abs(base_coef - reduced_coef),
      original_r2        = base_r2,
      reduced_r2         = reduced_r2,
      r2_loss            = r2_loss,
      contribution_score = contribution_score,
      PC                 = pc_num,
      model              = "Parity_Binary",
      test_term          = term_row,
      stringsAsFactors   = FALSE
    )
    df$external_gene_name <- gene_meta$external_gene_name[match(df$ensembl_id,
                                                                 gene_meta$genes)]
    df <- df[df$contribution_score > DRIVER_THRESH, ]
    df <- df[order(-df$contribution_score), ]
    df$final_rank <- seq_len(nrow(df))

    out_file <- file.path(out_dir,
                          sprintf("PC%d_Parity_Binary_improved_drivers.csv", pc_num))
    write.csv(df, out_file, row.names = FALSE)
    message(sprintf("[pc%d] %d driver genes (threshold = %.2f)",
                    pc_num, nrow(df), DRIVER_THRESH))
  }
}

if (sys.nframe() == 0L) main()
