# 05_create_umap.R — UMAP of harmonized cell types (Figure 2d-e).
# Stub: the full-reference UMAP requires access to raw single-cell data beyond
# the SCDC checkpoint. Complete reproduction requires the original reference
# H5s in Repository/data/references/. See Set 1/Final_Results_Summary/
# 2_Cell_Harmonization/UMAP_Visualization/create_final_harmonized_umap.R
# for the full (longer) implementation.

suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2)
  if (!requireNamespace("uwot", quietly = TRUE)) install.packages("uwot")
  library(uwot)
})
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

N_NEIGHBORS <- 30
MIN_DIST    <- 0.3

main <- function() {
  out_dir <- OUT$fig2
  aligned <- readRDS(file.path(out_dir, "aligned_signatures.rds"))
  mapping <- read.csv(file.path(out_dir, "refined_cell_type_mapping.csv"))

  sig_mat <- do.call(cbind, aligned)
  sig_mat <- sig_mat[, mapping$original_name[match(colnames(sig_mat), mapping$original_name)], drop = FALSE]

  set.seed(42)
  emb <- uwot::umap(t(sig_mat), n_neighbors = N_NEIGHBORS, min_dist = MIN_DIST)

  df <- data.frame(
    UMAP1         = emb[, 1],
    UMAP2         = emb[, 2],
    original_name = colnames(sig_mat),
    dataset       = mapping$dataset[match(colnames(sig_mat), mapping$original_name)],
    harmonized    = mapping$unified_name[match(colnames(sig_mat), mapping$original_name)],
    cell_class    = mapping$cell_class[match(colnames(sig_mat), mapping$original_name)]
  )

  p_harm <- ggplot(df, aes(UMAP1, UMAP2, color = harmonized)) +
    geom_point(size = 2, alpha = 0.9) +
    theme_classic(base_size = 10) +
    theme(legend.position = "right") +
    labs(title = "Figure 2e | Harmonized cell-type UMAP", color = NULL)

  p_class <- ggplot(df, aes(UMAP1, UMAP2, color = cell_class)) +
    geom_point(size = 2, alpha = 0.9) +
    theme_classic(base_size = 10) +
    labs(title = "Figure 2d | Cell class UMAP", color = NULL)

  ggsave(file.path(out_dir, "final_harmonized_umap_by_class.pdf"), p_class, width = 8, height = 6)
  ggsave(file.path(out_dir, "final_harmonized_umap_by_type.pdf"),  p_harm,  width = 10, height = 6)
  write.csv(df, file.path(out_dir, "umap_coordinates.csv"), row.names = FALSE)

  message(sprintf("[umap] Wrote UMAP PDFs for %d cell types", nrow(df)))
}

if (sys.nframe() == 0L) main()
