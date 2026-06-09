# 06_dotplot_figure4.R — Figure 4 dot plots for Age, Region, and Age×Region.
#
# Replaces the per-cell volcano panels (unreadable when more cells reach ANOVA
# significance) with three cell-type × pathway dot plots:
#   - Figure4a_age_dotplot.{pdf,png}
#   - Figure4b_region_dotplot.{pdf,png}
#   - Figure4c_interaction_dotplot.{pdf,png}
#   - Figure4_combined_dotplots.{pdf,png}   (composite: shared legend, balanced)
#
# Layout per panel:
#   x-axis: cell types (rotated 45°)
#   y-axis: pathway description (sorted by meta-score, highest at top)
#   colour: meta NES (blue → red)
#   size  : -log10(meta-FDR)
#
# Cell-level significance filter: only cells with at least one pathway at
# -log10(meta-FDR) >= MIN_NEG_LOG10_FDR_CELL are shown.
#
# Composite figure: all three panels share a single colour + size legend, with
# panel heights proportional to pathway count so dots stay visually consistent.

suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(patchwork)
})
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

MAX_DESC_CHARS         <- 45   # truncate long pathway descriptions
HIGH_CONF_TIERS        <- c("Ultra-High", "High", "Moderate")
MIN_NEG_LOG10_FDR_CELL <- 5    # cells must have at least one pathway above this
DOT_SIZE_RANGE         <- c(1.5, 7)
COLOUR_LOW             <- "#2c7bb6"   # low NES (blue)
COLOUR_HIGH            <- "#d7191c"   # high NES (red)

# ---------------------------------------------------------------------------
# Per-panel data preparation: returns the plot_data + bookkeeping needed by
# both standalone-plot and composite-plot rendering.
# ---------------------------------------------------------------------------
prepare_panel_data <- function(meta_results, effect_key, label_for_log) {
  eff_data <- meta_results %>% filter(effect == effect_key)
  if (nrow(eff_data) == 0) return(NULL)

  high_conf_pathways <- eff_data %>%
    filter(confidence_tier %in% HIGH_CONF_TIERS) %>%
    pull(pathway_id) %>%
    unique()
  if (length(high_conf_pathways) == 0) return(NULL)

  plot_data <- eff_data %>%
    filter(pathway_id %in% high_conf_pathways,
           !is.na(cell_type), nzchar(cell_type), cell_type != "NA",
           !is.na(description_short), nzchar(description_short),
           description_short != "NA",
           !is.na(meta_nes), !is.na(meta_fdr), !is.na(meta_score)) %>%
    mutate(neg_log10_fdr = -log10(pmax(meta_fdr, 1e-300)))

  # Cell-level significance filter
  qualifying_cells <- plot_data %>%
    group_by(cell_type) %>%
    summarise(max_sig = max(neg_log10_fdr, na.rm = TRUE), .groups = "drop") %>%
    filter(max_sig >= MIN_NEG_LOG10_FDR_CELL) %>%
    pull(cell_type)

  dropped <- setdiff(unique(plot_data$cell_type), qualifying_cells)
  if (length(dropped) > 0) {
    message(sprintf("[dotplot] %s: dropped %d cell(s) with max -log10(FDR) < %g: %s",
                    label_for_log, length(dropped),
                    MIN_NEG_LOG10_FDR_CELL,
                    paste(dropped, collapse = ", ")))
  }
  plot_data <- plot_data %>% filter(cell_type %in% qualifying_cells)
  if (nrow(plot_data) == 0) return(NULL)

  # Order cells by max meta_score (strongest signal cell on the left).
  cell_order <- plot_data %>%
    group_by(cell_type) %>%
    summarise(max_score = max(meta_score, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_score)) %>%
    pull(cell_type)
  plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_order)

  # Order pathways: highest max meta_score at top of y-axis.
  pathway_order <- plot_data %>%
    group_by(description_short) %>%
    summarise(max_score = max(meta_score, na.rm = TRUE), .groups = "drop") %>%
    arrange(max_score) %>%
    pull(description_short)
  plot_data$description_short <- factor(plot_data$description_short,
                                        levels = pathway_order)

  list(data = plot_data,
       n_pathways = length(pathway_order),
       n_cells = length(cell_order))
}

# ---------------------------------------------------------------------------
# Per-panel ggplot. `scales_shared` controls whether colour/size scales use
# global limits (TRUE → for composite figure) or per-panel limits (FALSE → for
# standalone files).
# ---------------------------------------------------------------------------
build_panel <- function(plot_data, title, nes_limits, sig_limits,
                        scales_shared = FALSE) {
  p <- ggplot(plot_data, aes(x = cell_type, y = description_short,
                             colour = meta_nes, size = neg_log10_fdr)) +
    geom_point(alpha = 0.85)

  if (scales_shared) {
    p <- p +
      scale_colour_gradient(name = "NES",
                            low = COLOUR_LOW, high = COLOUR_HIGH,
                            limits = nes_limits) +
      scale_size_continuous(name = expression(-log[10] * "(meta-FDR)"),
                            range = DOT_SIZE_RANGE,
                            limits = sig_limits)
  } else {
    p <- p +
      scale_colour_gradient(name = "NES",
                            low = COLOUR_LOW, high = COLOUR_HIGH) +
      scale_size_continuous(name = expression(-log[10] * "(meta-FDR)"),
                            range = DOT_SIZE_RANGE)
  }

  p +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(text             = element_text(colour = "black"),
          plot.title       = element_text(size = 14, face = "bold",
                                          colour = "black"),
          axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1,
                                          colour = "black", size = 11),
          axis.text.y      = element_text(size = 11, colour = "black"),
          axis.title       = element_text(colour = "black"),
          legend.title     = element_text(colour = "black"),
          legend.text      = element_text(colour = "black"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "gray92", linewidth = 0.3),
          legend.position  = "right")
}

main <- function() {
  out_dir <- OUT$fig4
  meta_results <- read.csv(file.path(REPO_ROOT, "checkpoints", "enrichment_fig4",
                                     "all_meta_analysis_results.csv"))

  meta_results <- meta_results %>%
    mutate(description_short = ifelse(nchar(description) > MAX_DESC_CHARS,
                                      paste0(substr(description, 1, MAX_DESC_CHARS - 3), "..."),
                                      description))

  effects <- list(
    age         = list(effect_key = "age",         title = "Age effects",
                       file_stem  = "Figure4a_age_dotplot"),
    region      = list(effect_key = "region",      title = "Region effects",
                       file_stem  = "Figure4b_region_dotplot"),
    interaction = list(effect_key = "age_region",  title = "Age × Region interaction",
                       file_stem  = "Figure4c_interaction_dotplot")
  )

  # Prepare data for every panel up front so we can compute global scale limits.
  panels <- list()
  for (eff_name in names(effects)) {
    cfg <- effects[[eff_name]]
    prep <- prepare_panel_data(meta_results, cfg$effect_key, eff_name)
    if (!is.null(prep)) {
      panels[[eff_name]] <- list(data = prep$data, cfg = cfg,
                                 n_pathways = prep$n_pathways,
                                 n_cells = prep$n_cells)
    }
  }

  if (length(panels) == 0) {
    message("[dotplot] no panels produced — nothing to plot")
    return(invisible(NULL))
  }

  # Global colour/size limits across all panels — used by the composite plot
  # so all three share one legend.
  all_data <- bind_rows(lapply(panels, `[[`, "data"))
  nes_limits <- range(all_data$meta_nes, na.rm = TRUE)
  sig_limits <- range(all_data$neg_log10_fdr, na.rm = TRUE)

  # ---- Standalone panels (per-panel auto-scaled, legends shown) -----------
  standalone_plots <- list()
  for (eff_name in names(panels)) {
    pan <- panels[[eff_name]]
    p <- build_panel(pan$data, pan$cfg$title,
                     nes_limits = NULL, sig_limits = NULL,
                     scales_shared = FALSE)

    h <- max(4, pan$n_pathways * 0.25 + 2)
    w <- max(6, pan$n_cells * 0.8 + 4)
    ggsave(file.path(out_dir, paste0(pan$cfg$file_stem, ".pdf")), p,
           width = w, height = h, limitsize = FALSE)
    ggsave(file.path(out_dir, paste0(pan$cfg$file_stem, ".png")), p,
           width = w, height = h, dpi = 300, limitsize = FALSE)

    standalone_plots[[eff_name]] <- p
    message(sprintf("[dotplot] %-12s : %d pathways x %d cells (standalone)",
                    eff_name, pan$n_pathways, pan$n_cells))
  }

  # ---- Composite figure (shared scales, single legend) --------------------
  composite_panels <- lapply(names(panels), function(eff_name) {
    pan <- panels[[eff_name]]
    build_panel(pan$data, pan$cfg$title,
                nes_limits = nes_limits, sig_limits = sig_limits,
                scales_shared = TRUE)
  })

  # Vertical stack with heights proportional to pathway counts so dots stay
  # visually consistent across panels. Collect legends into a single shared
  # legend on the right.
  heights_for_layout <- sapply(panels, function(p) max(2, p$n_pathways))
  composite <- wrap_plots(composite_panels, ncol = 1) +
    plot_layout(guides  = "collect",
                heights = heights_for_layout) &
    theme(legend.position = "right",
          legend.box      = "vertical")

  # Total composite height: ~0.25 inch per pathway row + ~3 inch overhead for
  # the three panel titles + per-panel x-axis label spacing.
  total_pathways <- sum(heights_for_layout)
  composite_h <- max(8, total_pathways * 0.25 + 4)
  composite_w <- max(7, max(sapply(panels, `[[`, "n_cells")) * 0.9 + 5)

  ggsave(file.path(out_dir, "Figure4_combined_dotplots.pdf"), composite,
         width = composite_w, height = composite_h, limitsize = FALSE)
  ggsave(file.path(out_dir, "Figure4_combined_dotplots.png"), composite,
         width = composite_w, height = composite_h, dpi = 300,
         limitsize = FALSE)

  message(sprintf("[dotplot] composite figure: %d panels, %d total pathway rows, %.1f x %.1f inches",
                  length(panels), total_pathways, composite_w, composite_h))
  message("[dotplot] wrote Figure 4 dot plots to ", out_dir)
}

if (sys.nframe() == 0L) main()
