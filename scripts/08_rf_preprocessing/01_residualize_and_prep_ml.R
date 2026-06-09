# 06c/01_residualize_and_prep_ml.R — regenerate the RF ML-data checkpoints.
#
# Ports two legacy scripts into a single config-driven step:
#   - residualized_parity_analysis.R     (regress out age + region per gene)
#   - random_forest_parity_analysis.R    (lines 22-48: transpose + metadata for ML)
#
# Outputs (overwrites the staged checkpoint so reruns regenerate it):
#   checkpoints/rf_inputs/ml_expression_data.csv
#   checkpoints/rf_inputs/ml_metadata.csv
#
# The §07 random-forest analysis (`09_suppfig2to5_parity_genes/02_random_forest.R`)
# reads from this checkpoint, so this section closes the regeneration gap for
# Supp Fig 5 random-forest top-23 gene list.

suppressPackageStartupMessages({
  library(tidyverse)
})
source(Sys.getenv(
  "REPRO_CONFIG",
  "f:/Parity/Final/Repository/scripts/config.R"
))

RF_CKPT_DIR <- file.path(CHECKPOINT_DIR, "rf_inputs")

# ---- 1. Residualize age + region per gene --------------------------------
#
# For each gene g: fit lm(g ~ age + region), take residuals + mean of g.
# The "+ mean" preserves the gene's baseline expression so RF input remains
# on the same scale as the original; only age- and region-driven variance is
# removed. Parity effects (and parity:age, parity:region interactions) are
# fully retained.

residualize_one_gene <- function(gene_expr, age, region) {
  fit <- lm(gene_expr ~ age + region)
  residuals(fit) + mean(gene_expr)
}

residualize_all <- function(expr_matrix, metadata) {
  # expr_matrix: genes (rows) x samples (cols)
  # metadata$sample must align to colnames(expr_matrix)
  stopifnot(all(metadata$sample == colnames(expr_matrix)))

  age    <- factor(metadata$age)
  region <- factor(metadata$region)

  out <- matrix(NA_real_,
                nrow = nrow(expr_matrix),
                ncol = ncol(expr_matrix),
                dimnames = dimnames(expr_matrix))

  for (i in seq_len(nrow(expr_matrix))) {
    g <- as.numeric(expr_matrix[i, ])
    out[i, ] <- residualize_one_gene(g, age, region)
    if (i %% 2000 == 0) message(sprintf("[resid] %d / %d genes",
                                        i, nrow(expr_matrix)))
  }
  out
}

# ---- 2. Transpose to ML format + write metadata --------------------------

prep_ml_data <- function(residualized, metadata, out_dir) {
  ml_data <- t(residualized)  # samples (rows) x genes (cols)

  # Save expression: rownames = sample, colnames = ENSRNOG gene IDs.
  # write.csv default emits "" for the row-name header, which matches
  # the staged checkpoint format.
  out_expr <- file.path(out_dir, "ml_expression_data.csv")
  write.csv(ml_data, out_expr)
  message(sprintf("[ml] wrote %s (%d samples x %d genes)",
                  out_expr, nrow(ml_data), ncol(ml_data)))

  # ML metadata: per-sample factors RF needs.
  ml_meta <- data.frame(
    sample        = metadata$sample,
    parity        = metadata$parity,
    parity_binary = ifelse(metadata$parity == "Nulliparous",
                           "Nulliparous", "Parous"),
    parity_3group = metadata$parity,
    stringsAsFactors = FALSE
  )
  out_meta <- file.path(out_dir, "ml_metadata.csv")
  write.csv(ml_meta, out_meta, row.names = FALSE)
  message(sprintf("[ml] wrote %s (%d samples)", out_meta, nrow(ml_meta)))

  invisible(list(ml_data = ml_data, ml_meta = ml_meta))
}

# ---- main -----------------------------------------------------------------

main <- function() {
  dir.create(RF_CKPT_DIR, recursive = TRUE, showWarnings = FALSE)

  expr_path <- file.path(BULK_DIR, "filtered_expression_for_anova.csv")
  meta_path <- file.path(BULK_DIR, "aligned_sample_metadata.csv")

  message(sprintf("[ml] loading expression: %s", expr_path))
  expr <- read.csv(expr_path, row.names = 1, check.names = FALSE)
  expr_matrix <- as.matrix(expr)

  message(sprintf("[ml] loading metadata:    %s", meta_path))
  metadata <- read.csv(meta_path, stringsAsFactors = FALSE)
  metadata <- metadata[match(colnames(expr_matrix), metadata$sample), ,
                       drop = FALSE]

  message(sprintf("[ml] %d genes x %d samples", nrow(expr_matrix),
                  ncol(expr_matrix)))

  message("[resid] residualizing age + region per gene")
  residualized <- residualize_all(expr_matrix, metadata)

  prep_ml_data(residualized, metadata, RF_CKPT_DIR)

  message("[06c] RF preprocessing complete")
  invisible(NULL)
}

if (sys.nframe() == 0L) main()
