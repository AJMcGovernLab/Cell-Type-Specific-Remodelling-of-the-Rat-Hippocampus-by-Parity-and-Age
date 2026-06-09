# 03_plot_parity_volcano.R — Figure 5 parity volcano rebuilt from REAL log2 fold changes.
#
# Replaces the fabricated `01_create_parity_volcano_plots.R` (which used runif()
# bins + grepl(cell_type) for direction). Two outputs:
#   - parity_volcano_per_dataset.pdf  (faceted by mouse10x_2020 / smartseq_2019 / yao_hippo_10x)
#   - parity_volcano_unified.pdf      (all 27 cell × dataset combinations together)

suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(ggrepel)
})
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

P_THRESH  <- 0.05
FC_THRESH <- 0.05   # smaller than fig3 — parity proportional changes are typically subtle
FC_CAP    <- 1.5

cap_fc <- function(x, cap = FC_CAP) pmin(pmax(x, -cap), cap)

theme_volcano <- theme_minimal(base_size = 12) +
  theme(plot.title       = element_text(size = 14, face = "bold"),
        strip.text       = element_text(size = 11, face = "bold"),
        legend.position  = "bottom",
        plot.margin      = unit(c(1, 1, 1.4, 1), "cm"))

faceted_panel <- function(d) {
  ggplot(d, aes(x = log2fc_parity, y = neg_log10_p)) +
    geom_hline(yintercept = -log10(P_THRESH), linetype = "dashed", colour = "gray50") +
    geom_vline(xintercept = 0,                 linetype = "solid",  colour = "black", alpha = 0.3) +
    geom_vline(xintercept = c(-FC_THRESH, FC_THRESH), linetype = "dashed", colour = "gray50") +
    geom_point(aes(colour = direction, size = abs(log2fc_parity)), alpha = 0.8) +
    geom_text_repel(data = . %>% filter(parity_significant),
                    aes(label = cell_type), size = 2.7,
                    box.padding = 0.35, point.padding = 0.3,
                    max.overlaps = 30, seed = 42) +
    scale_colour_manual(values = c("Higher in parous" = "#d73027",
                                   "Higher in nullip" = "#4575b4",
                                   "Not Significant"  = "#969696"),
                        name = "Parity effect") +
    scale_size_continuous(name = "Effect size", range = c(2.2, 7.5)) +
    facet_wrap(~ dataset, scales = "free", ncol = 3) +
    labs(title = "Parity effects on cell-type proportions (per dataset)",
         x = "log2 FC (Parous vs Nulliparous)",
         y = "-log10(p-value)") +
    theme_volcano
}

unified_panel <- function(d) {
  d_lab <- d %>% filter(parity_significant)

  # Piecewise y-axis transformation: compresses the 0-1.5 region into the
  # bottom ~15% of the plot and expands the 1.5+ region into the upper ~85%,
  # so the four parity-significant cells (-log10 p ≈ 1.55–2.07) are clearly
  # separated.
  #
  #   y in [0, 1.5]  → display 0      to 0.15  (factor 0.1)
  #   y in [1.5, ∞)  → display 0.15   to 0.15 + (y - 1.5)
  sq_trans <- scales::trans_new(
    name      = "piecewise_top",
    transform = function(x) ifelse(x <= 1.5, 0.1 * x, 0.15 + (x - 1.5)),
    inverse   = function(y) ifelse(y <= 0.15, 10 * y, 1.5 + y - 0.15)
  )

  ymax_data <- max(d$neg_log10_p, na.rm = TRUE)

  ggplot(d, aes(x = log2fc_parity, y = neg_log10_p)) +
    geom_vline(xintercept = 0, colour = "black", alpha = 0.3) +
    geom_point(aes(colour = direction,
                   size   = parity_significant),
               alpha = 0.85) +
    geom_text_repel(data = d_lab,
                    aes(label = paste0(cell_type, " (", dataset, ")")),
                    size            = 5,
                    fontface        = "bold",
                    colour          = "black",
                    box.padding     = 0.45,
                    point.padding   = 0.35,
                    max.overlaps    = 30,
                    min.segment.length = 0,
                    seed            = 42) +
    scale_colour_manual(values = c("Higher in parous" = "#d73027",
                                   "Higher in nullip" = "#4575b4",
                                   "Not Significant"  = "#bdbdbd")) +
    scale_size_manual(values = c(`TRUE` = 6, `FALSE` = 2)) +
    scale_y_continuous(transform = sq_trans,
                       breaks    = c(0, 0.5, 1, 1.5, 2),
                       limits    = c(0, max(2.1, ymax_data) * 1.02)) +
    guides(colour = "none", size = "none") +
    coord_cartesian(xlim = c(-0.125, 0.125), clip = "off") +
    labs(title = "Parity effects on cell-type proportions",
         x = "log2 FC (Parous - Nulliparous)",
         y = "-log10(p-value)") +
    theme_volcano +
    theme(text             = element_text(colour = "black"),
          plot.title       = element_text(size = 14, face = "bold",
                                          colour = "black"),
          axis.title       = element_text(size = 12, face = "bold",
                                          colour = "black"),
          axis.text        = element_text(size = 11, colour = "black"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "gray92",
                                          linewidth = 0.3))
}

main <- function() {
  out_dir <- OUT$fig5ab
  d <- read.csv(file.path(out_dir, "per_dataset_log2fc.csv")) %>%
    mutate(
      log2fc_parity        = cap_fc(log2fc_parity),
      neg_log10_p          = -log10(pmax(parity_pval, 1e-300)),
      parity_significant   = parity_pval < P_THRESH,
      # Direction now uses just the p-value and sign — colour every
      # parity-significant cell regardless of effect-size magnitude so the
      # plot is visually consistent with the §3.4.1 hit list (4 cells).
      direction = case_when(
        parity_significant & log2fc_parity >= 0 ~ "Higher in parous",
        parity_significant & log2fc_parity <  0 ~ "Higher in nullip",
        TRUE                                    ~ "Not Significant"
      )
    )

  facet <- faceted_panel(d)
  uni   <- unified_panel(d)

  ggsave(file.path(out_dir, "parity_volcano_per_dataset.pdf"), facet,
         width = 13, height = 6)
  ggsave(file.path(out_dir, "parity_volcano_unified.pdf"),     uni,
         width = 9, height = 7.5)
  ggsave(file.path(out_dir, "parity_volcano_per_dataset.png"), facet,
         width = 13, height = 6, dpi = 300)
  ggsave(file.path(out_dir, "parity_volcano_unified.png"),     uni,
         width = 9, height = 7.5, dpi = 300)

  message("[plot-fig5] wrote parity volcanos to ", out_dir)
  invisible(list(facet = facet, uni = uni))
}

if (sys.nframe() == 0L) main()
