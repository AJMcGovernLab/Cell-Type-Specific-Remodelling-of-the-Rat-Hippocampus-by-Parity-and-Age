# 06_parity_focused_heatmap.R — Figure 5a (parity-focused) heatmap.
#
# A subset of `Figure5a_comprehensive_heatmap.png`:
#   - columns: Parity, Age × Parity, Parity × Region, Three-way interaction
#   - rows   : only cells where at least one of those four p-values < 0.10
#   - colour : -log10(p), 0 (blue) → 1 (white) → 2 (red), values above 2 saturate
#
# Writes:
#   - outputs/06_fig5ab_parity_proportions/Figure5a_parity_focused_heatmap.{pdf,png}

suppressPackageStartupMessages({
  library(tidyverse)
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE)
  library(ComplexHeatmap); library(circlize)
})
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

P_THRESH         <- 0.10   # row filter (cell must have any p below this)
P_SIG            <- 0.05   # cell outline + white text threshold

# Five-stop -log10(p) colour ramp:
#   0 (p=1)    deep blue   #303595
#   0.5        light blue  #73ABD0
#   1 (p=0.10) white
#   1.5        orange      #FEB365
#   2 (p=0.01) deep red    #A50025
COL_BREAKS  <- c(0, 0.5, 1, 1.5, 2)
COL_RAMP    <- c("#303595", "#73ABD0", "#FFFFFF", "#FEB365", "#A50025")

main <- function() {
  out_dir <- OUT$fig5ab
  df <- read.csv(file.path(out_dir, "fixed_individual_anova_results.csv"))

  # Keep only the four parity-related effect columns.
  # Display order: interaction terms first, Parity main effect on the right.
  parity_effects <- c("Age × Parity"         = "age_parity_pval",
                      "Parity × Region"      = "parity_region_pval",
                      "Three-way"            = "three_way_pval",
                      "Parity"               = "parity_pval")

  long <- df %>%
    mutate(row = sprintf("%s | %s", dataset, cell_type)) %>%
    select(row, all_of(parity_effects)) %>%
    pivot_longer(-row, names_to = "effect", values_to = "pval") %>%
    mutate(effect      = recode(effect, !!!setNames(names(parity_effects),
                                                    unname(parity_effects))),
           neg_log10_p = -log10(pmax(pval, 1e-300)))

  # Cells passing the threshold in at least one of the four parity effects
  keep_rows <- long %>%
    group_by(row) %>%
    summarise(min_p = min(pval, na.rm = TRUE), .groups = "drop") %>%
    filter(min_p < P_THRESH) %>%
    pull(row)

  message(sprintf("[fig5a-parity] keeping %d of %d cells with min(parity p) < %.2f",
                  length(keep_rows), n_distinct(long$row), P_THRESH))

  if (length(keep_rows) == 0) {
    message("[fig5a-parity] no cells passed threshold; aborting")
    return(invisible(NULL))
  }

  # Wide -log10(p) matrix (used for cell colour)
  mat <- long %>%
    filter(row %in% keep_rows) %>%
    select(row, effect, neg_log10_p) %>%
    pivot_wider(names_from = effect, values_from = neg_log10_p) %>%
    column_to_rownames("row") %>%
    as.matrix()
  mat[!is.finite(mat)] <- 0

  # Parallel raw-p matrix (used for cell text labels)
  pmat <- long %>%
    filter(row %in% keep_rows) %>%
    select(row, effect, pval) %>%
    pivot_wider(names_from = effect, values_from = pval) %>%
    column_to_rownames("row") %>%
    as.matrix()

  # Order columns explicitly so left→right reads
  # interactions → Parity main effect (right-most). Rows are clustered by
  # ComplexHeatmap (dendrogram appears on the left).
  mat  <- mat[, names(parity_effects)]
  pmat <- pmat[, names(parity_effects)]

  col_fun <- colorRamp2(COL_BREAKS, COL_RAMP)

  # Column split: put the Parity main effect in its own group so there is a
  # visual gap between it and the three interaction columns on the left.
  col_groups <- factor(ifelse(colnames(mat) == "Parity",
                              "main", "interaction"),
                       levels = c("interaction", "main"))

  # Format p-values for in-cell labels: 3 decimal places for p >= 0.001,
  # otherwise scientific notation with 2 sig figs.
  format_p <- function(p) {
    if (is.na(p))                     return("")
    if (p < 0.001)                    return(formatC(p, format = "e", digits = 2))
    formatC(p, format = "f", digits = 3)
  }

  ht <- Heatmap(mat,
                name             = "-log10(p)",
                col              = col_fun,
                cluster_rows     = TRUE,
                cluster_columns  = FALSE,
                column_split     = col_groups,
                column_gap       = unit(3, "mm"),
                row_dend_side    = "left",
                row_dend_width   = unit(15, "mm"),
                row_names_side   = "right",
                row_names_gp     = gpar(fontsize = 10, col = "black"),
                column_names_gp  = gpar(fontsize = 11, col = "black"),
                column_names_rot = 45,
                rect_gp          = gpar(col = "white", lwd = 0.5),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  p <- pmat[i, j]
                  # White text for p < 0.05, black otherwise
                  txt_col <- if (!is.na(p) && p < P_SIG) "white" else "black"
                  grid.text(format_p(p), x, y,
                            gp = gpar(fontsize = 10, col = txt_col))
                  # Black outline around p < 0.05 cells (overdraws the white
                  # default border for that cell only).
                  if (!is.na(p) && p < P_SIG) {
                    grid.rect(x, y, width, height,
                              gp = gpar(col = "black", fill = NA, lwd = 1.5))
                  }
                },
                heatmap_legend_param = list(
                  at         = c(0, 0.5, 1, 1.5, 2),
                  labels     = c("0", "0.5", "1", "1.5", "≥2"),
                  title_gp   = gpar(fontsize = 11, fontface = "bold",
                                    col = "black"),
                  labels_gp  = gpar(fontsize = 10, col = "black")
                ),
                column_title     = NULL,
                row_title        = NULL)

  pdf_path <- file.path(out_dir, "Figure5a_parity_focused_heatmap.pdf")
  png_path <- file.path(out_dir, "Figure5a_parity_focused_heatmap.png")
  fig_h    <- max(4, nrow(mat) * 0.32 + 2)
  fig_w    <- 7

  pdf(pdf_path, width = fig_w, height = fig_h)
  draw(ht, padding = unit(c(8, 4, 4, 4), "mm"))
  dev.off()

  png(png_path, width = fig_w, height = fig_h, units = "in", res = 300)
  draw(ht, padding = unit(c(8, 4, 4, 4), "mm"))
  dev.off()

  message(sprintf("[fig5a-parity] wrote %s (%d × %d)", pdf_path,
                  nrow(mat), ncol(mat)))
}

if (sys.nframe() == 0L) main()
