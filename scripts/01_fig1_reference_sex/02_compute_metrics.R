# 02_compute_metrics.R — normalized_metrics_table.csv (source of Figure 1 numbers).
# Condensed from Set 1/Transfer/Final/generate_nature_biotech_figures_final_v4.R (lines 30-90).

suppressPackageStartupMessages(library(tidyverse))
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

DATASET_LABEL <- c(
  mouse10x_2020       = "Mouse 10x 2020",
  yao_hippo_10x       = "Yao Hippo 10x",
  mouse_smartseq_2019 = "Mouse Smart-seq"
)

compute_metrics <- function(all_results, summary_data) {
  rows <- list()

  for (ref_name in names(all_results)) {
    res  <- all_results[[ref_name]]
    props <- res$proportions
    prop_cols <- setdiff(colnames(props), c("sample", "sample_id"))
    mat <- as.matrix(props[, prop_cols])

    # Parse "dataset_sex" from ref_name (dataset may contain underscores)
    parts   <- strsplit(ref_name, "_")[[1]]
    sex     <- parts[length(parts)]
    dataset <- paste(parts[-length(parts)], collapse = "_")

    shannon <- apply(mat, 1, function(x) {
      x <- x[x > 0]
      if (length(x) <= 1) 0 else -sum(x * log2(x))
    })

    n_cells <- summary_data$n_cells[summary_data$reference == ref_name]

    rows[[ref_name]] <- data.frame(
      reference                 = ref_name,
      dataset                   = dataset,
      dataset_label             = unname(DATASET_LABEL[dataset]),
      sex                       = factor(sex, levels = c("male","female","mixed")),
      technology                = ifelse(grepl("smartseq", dataset), "Smart-seq2", "10X Genomics"),
      n_ref_cells               = n_cells,
      n_detected_celltypes      = ncol(mat),
      types_per_k_cells         = (ncol(mat) / n_cells) * 1000,
      entropy_per_k_cells       = (mean(shannon) / n_cells) * 1000,
      sparsity_per_k_cells      = (mean(mat == 0) * 100 / n_cells) * 1000,
      mean_prop_per_k_cells     = (mean(mat[mat > 0]) / n_cells) * 1000,
      detected_per_sample_per_k = (mean(apply(mat > 0.01, 1, sum)) / n_cells) * 1000,
      diversity_index           = mean(shannon) / log2(ncol(mat)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

main <- function() {
  out_dir <- OUT$fig1
  # Prefer the outputs/ copy (written by 01_run_deconvolution.R);
  # otherwise fall back to the pre-computed SCDC checkpoints so we can
  # reproduce Figure 1 without re-running SCDC.
  rds <- file.path(out_dir, "all_results_scdc.rds")
  csv <- file.path(out_dir, "deconvolution_summary.csv")
  if (!file.exists(rds)) rds <- file.path(SCDC_CHECKPOINTS, "all_results_scdc.rds")
  if (!file.exists(csv)) csv <- file.path(SCDC_CHECKPOINTS, "deconvolution_summary.csv")
  all_results  <- readRDS(rds)
  summary_data <- read.csv(csv)

  metrics_df <- compute_metrics(all_results, summary_data)
  write.csv(metrics_df, file.path(out_dir, "normalized_metrics_table.csv"), row.names = FALSE)
  message(sprintf("[metrics] Wrote %s", file.path(out_dir, "normalized_metrics_table.csv")))
  invisible(metrics_df)
}

if (sys.nframe() == 0L) main()
