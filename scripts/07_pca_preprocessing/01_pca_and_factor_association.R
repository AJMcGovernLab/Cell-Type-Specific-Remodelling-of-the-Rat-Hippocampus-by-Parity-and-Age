# 06b/01_pca_and_factor_association.R — regenerate the PCA driver checkpoints.
#
# Ports three legacy scripts into a single config-driven step:
#   - pca_analysis_clean.R                (PCA on clean_normalized_counts.csv)
#   - create_binary_parity_pca.R          (PC coordinates + binary parity metadata)
#   - factor_association_analysis.R       (lm models per PC × factor combination)
#
# Outputs (overwrites the staged checkpoint so reruns regenerate it):
#   checkpoints/pca_drivers/pca_loadings_top10.csv
#   checkpoints/pca_drivers/binary_parity_pca_coordinates.csv
#   checkpoints/pca_drivers/factor_association_significant_models_summary.csv
#
# The §07 PCA driver-gene analysis (`09_suppfig2to5_parity_genes/03_pc_driver_analysis.R`)
# reads from this checkpoint, so this section closes the regeneration gap for
# Supp Fig 4 driver-gene counts.

suppressPackageStartupMessages({
  library(tidyverse)
})
source(Sys.getenv(
  "REPRO_CONFIG",
  "f:/Parity/Final/Repository/scripts/config.R"
))

PCA_CKPT_DIR <- file.path(CHECKPOINT_DIR, "pca_drivers")
N_PCS_SAVE   <- 10

# ---- 1. Load expression data and metadata ---------------------------------

run_pca <- function(counts_path, meta_path) {
  message(sprintf("[pca] loading expression: %s", counts_path))
  counts <- read.csv(counts_path, check.names = FALSE)
  message(sprintf("[pca] loading metadata:    %s", meta_path))
  metadata <- read.csv(meta_path, stringsAsFactors = FALSE)

  # Gene names: use Gene_1, Gene_2, ... (matches legacy convention used by
  # downstream §07 driver-gene analysis).
  rownames(counts) <- paste0("Gene_", seq_len(nrow(counts)))
  message(sprintf("[pca] %d genes x %d samples", nrow(counts), ncol(counts)))

  # Align metadata to count columns
  metadata <- metadata[match(colnames(counts), metadata$sample), , drop = FALSE]
  stopifnot(all(metadata$sample == colnames(counts)))

  # Add binary parity
  metadata$parity_binary <- ifelse(metadata$parity == "Nulliparous",
                                   "Nulliparous", "Parous")

  # PCA: samples as rows, genes as columns; center + scale
  pca_input <- t(as.matrix(counts))
  message("[pca] running prcomp (center + scale)")
  pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)

  list(pca = pca_result, metadata = metadata)
}

# ---- 2. Save loadings for top N PCs ---------------------------------------

save_loadings <- function(pca, out_dir, n_pcs = N_PCS_SAVE) {
  loadings <- as.data.frame(pca$rotation[, seq_len(n_pcs)])
  loadings$gene <- rownames(loadings)
  loadings$gene_index <- seq_len(nrow(loadings))
  out <- file.path(out_dir, "pca_loadings_top10.csv")
  write.csv(loadings, out, row.names = FALSE)
  message(sprintf("[pca] wrote %s (%d genes x %d PCs)",
                  out, nrow(loadings), n_pcs))
}

# ---- 3. Save PC coordinates with metadata + binary parity -----------------

save_coordinates <- function(pca, metadata, out_dir, n_pcs = N_PCS_SAVE) {
  coords <- as.data.frame(pca$x[, seq_len(n_pcs)])
  coords$sample <- rownames(coords)
  coords <- coords %>% left_join(metadata, by = "sample")
  out <- file.path(out_dir, "binary_parity_pca_coordinates.csv")
  write.csv(coords, out, row.names = FALSE)
  message(sprintf("[pca] wrote %s (%d samples x %d PCs)",
                  out, nrow(coords), n_pcs))
  invisible(coords)
}

# ---- 4. Factor-association lm models --------------------------------------
#
# For each of the top 10 PCs, fit 7 lm models testing different factor
# combinations. Filter to significant model fits and write the summary
# the §07 driver-gene script consumes.

MODELS <- list(
  Age           = "~ age",
  Region        = "~ region",
  Parity_Binary = "~ parity_binary",
  Age_Region    = "~ age + region + age:region",
  Age_Parity    = "~ age + parity_binary + age:parity_binary",
  Region_Parity = "~ region + parity_binary + region:parity_binary",
  Three_Way     = paste("~ age + region + parity_binary +",
                        "age:region + age:parity_binary +",
                        "region:parity_binary + age:region:parity_binary")
)

fit_one_pc <- function(pc_name, pc_data) {
  out <- data.frame()
  for (model_name in names(MODELS)) {
    f <- as.formula(paste(pc_name, MODELS[[model_name]]))
    fit <- tryCatch(lm(f, data = pc_data), error = function(e) NULL)
    if (is.null(fit)) next
    s <- summary(fit)
    f_stat <- s$fstatistic
    if (is.null(f_stat)) next
    p_val <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
    out <- rbind(out, data.frame(
      PC            = pc_name,
      Model         = model_name,
      R_squared     = s$r.squared,
      Model_p_value = unname(p_val)
    ))
  }
  out
}

run_factor_association <- function(coords, out_dir, n_pcs = N_PCS_SAVE) {
  pc_data <- coords %>%
    mutate(
      age           = as.factor(age),
      region        = as.factor(region),
      parity_binary = as.factor(parity_binary)
    )

  message(sprintf("[fa] fitting 7 models x %d PCs", n_pcs))
  pc_names <- paste0("PC", seq_len(n_pcs))
  all_fits <- bind_rows(lapply(pc_names, fit_one_pc, pc_data = pc_data))

  # Filter to model fits with p < 0.05; format columns to match legacy schema.
  sig <- all_fits %>%
    filter(Model_p_value < 0.05) %>%
    arrange(PC, desc(R_squared)) %>%
    mutate(
      R_squared_percent = sprintf("%.1f%%", 100 * R_squared),
      p_value_formatted = formatC(Model_p_value, format = "e", digits = 2)
    ) %>%
    select(PC, Model, R_squared, Model_p_value,
           R_squared_percent, p_value_formatted)

  out <- file.path(out_dir, "factor_association_significant_models_summary.csv")
  write.csv(sig, out, row.names = FALSE)
  message(sprintf("[fa] wrote %s (%d significant model fits)",
                  out, nrow(sig)))

  invisible(sig)
}

# ---- main -----------------------------------------------------------------

main <- function() {
  dir.create(PCA_CKPT_DIR, recursive = TRUE, showWarnings = FALSE)

  counts_path <- file.path(BULK_DIR, "clean_normalized_counts.csv")
  meta_path   <- file.path(BULK_DIR, "sample_metadata.csv")

  res <- run_pca(counts_path, meta_path)
  save_loadings(res$pca, PCA_CKPT_DIR)
  coords <- save_coordinates(res$pca, res$metadata, PCA_CKPT_DIR)
  run_factor_association(coords, PCA_CKPT_DIR)

  message("[06b] PCA preprocessing complete")
  invisible(NULL)
}

if (sys.nframe() == 0L) main()
