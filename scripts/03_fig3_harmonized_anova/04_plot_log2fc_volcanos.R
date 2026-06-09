# 04_plot_log2fc_volcanos.R — Figure 3 panels rebuilt from REAL log2 fold changes.
#
# Replaces the fabricated `08_create_volcano_plots.R` with three panels driven
# by harmonized_log2fc.csv:
#   - age_effects_volcano.pdf       (Figure 3a)
#   - region_effects_volcano.pdf    (Figure 3b)
#   - interaction_effects_plot.pdf  (Figure 3c)
#
# Refresh 2026-05-24:
#   * Axis convention kept as Old vs Young / Dorsal vs Ventral
#     (positive = Old higher / Dorsal higher). Axis labels now written as
#     "Old - Young" and "Dorsal - Ventral" so direction is unambiguous on
#     inspection: positive numbers are Old-enriched / Dorsal-enriched,
#     negative numbers are Young-enriched / Ventral-enriched.
#   * Direction colour legend made explicit: red (right) = "Old-enriched" /
#     "Dorsal-enriched", blue (left) = "Young-enriched" / "Ventral-enriched".
#   * Interaction plot tightened to ±1 on both axes (was ±2). Cells with
#     |log2fc_age| > 1 OR |log2fc_region| > 1 are dynamically excluded
#     (was a hardcoded outlier list).
#
# Refresh 2026-05-10:
#   * Larger labels on significant cells (ggrepel size 4 → 4.5)
#   * Arrows + YOUNG/OLD/VENTRAL/DORSAL labels removed (axis label is enough)
#   * Legend title shortened so it doesn't repeat the main title

suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(ggrepel)
})
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

P_THRESH        <- 0.05
FC_THRESH       <- 0.2
FC_CAP          <- 4    # squash extreme log2FC for the volcano panels
FC_CAP_INTER    <- 1    # tighter cap for the age:region interaction quadrant
LABEL_SIZE      <- 4.5  # ggrepel label size for significant cells

cap_fc <- function(x, cap = FC_CAP) pmin(pmax(x, -cap), cap)

base_volcano <- function(d, x_col, p_col, sig_col, dir_levels, title, x_lab,
                         right_col = "#d73027", left_col = "#4575b4",
                         ns_col = "#969696") {
  d <- d %>% mutate(
    !!x_col := cap_fc(.data[[x_col]]),
    neg_log10_p = -log10(pmax(.data[[p_col]], 1e-300)),
    direction   = case_when(
      .data[[sig_col]] & .data[[x_col]] >  FC_THRESH ~ dir_levels[1],
      .data[[sig_col]] & .data[[x_col]] < -FC_THRESH ~ dir_levels[2],
      TRUE                                          ~ "Not Significant"
    )
  )

  cols <- setNames(c(right_col, left_col, ns_col),
                   c(dir_levels[1], dir_levels[2], "Not Significant"))

  ymax <- max(d$neg_log10_p, na.rm = TRUE) * 1.05

  ggplot(d, aes(x = .data[[x_col]], y = neg_log10_p)) +
    geom_hline(yintercept = -log10(P_THRESH), linetype = "dashed",
               colour = "gray50") +
    geom_vline(xintercept = 0, linetype = "solid",  colour = "black",
               alpha = 0.3) +
    geom_vline(xintercept = c(-FC_THRESH, FC_THRESH),
               linetype = "dashed", colour = "gray50") +
    geom_point(aes(colour = direction, size = abs(.data[[x_col]])),
               alpha = 0.75) +
    geom_text_repel(data = d %>% filter(.data[[sig_col]]),
                    aes(label = cell_type),
                    size            = LABEL_SIZE,
                    fontface        = "bold",
                    box.padding     = 0.45,
                    point.padding   = 0.4,
                    max.overlaps    = 30,
                    seed            = 42) +
    scale_colour_manual(values = cols, name = "Direction") +
    scale_size_continuous(name = "Effect size", range = c(2.2, 7.5)) +
    coord_cartesian(clip = "off", ylim = c(0, ymax)) +
    labs(title = title, x = x_lab, y = "-log10(p-value)") +
    theme_minimal(base_size = 12) +
    theme(plot.title       = element_text(size = 14, face = "bold"),
          legend.position  = "bottom",
          plot.margin      = unit(c(1, 1, 1, 1), "cm"))
}

interaction_plot <- function(d) {
  # Dynamically exclude cells outside the ±1 plotting window on either axis.
  excluded <- d %>%
    filter(abs(log2fc_age) > FC_CAP_INTER | abs(log2fc_region) > FC_CAP_INTER) %>%
    pull(cell_type)
  if (length(excluded) > 0) {
    message(sprintf("[plot-fig3] interaction plot excludes %d cells outside ±%g: %s",
                    length(excluded), FC_CAP_INTER, paste(excluded, collapse = ", ")))
  }

  d <- d %>%
    filter(abs(log2fc_age) <= FC_CAP_INTER, abs(log2fc_region) <= FC_CAP_INTER) %>%
    mutate(
      interaction_significant = age_region_pval < P_THRESH,
      interaction_neg_log10_p = -log10(pmax(age_region_pval, 1e-300))
    )

  ggplot(d, aes(x = log2fc_age, y = log2fc_region)) +
    geom_hline(yintercept = 0, colour = "gray30", alpha = 0.5) +
    geom_vline(xintercept = 0, colour = "gray30", alpha = 0.5) +
    geom_point(aes(colour = interaction_neg_log10_p,
                   size   = interaction_significant), alpha = 0.85) +
    geom_text_repel(data = d %>% filter(interaction_significant),
                    aes(label = cell_type),
                    size            = LABEL_SIZE,
                    fontface        = "bold",
                    box.padding     = 0.45,
                    point.padding   = 0.4,
                    max.overlaps    = 30,
                    seed            = 42) +
    scale_colour_viridis_c(name = "-log10(age:region p)", option = "plasma") +
    scale_size_manual(values = c(`TRUE` = 4.5, `FALSE` = 2.4),
                      labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s."),
                      name = "age:region") +
    coord_cartesian(xlim = c(-FC_CAP_INTER, FC_CAP_INTER),
                    ylim = c(-FC_CAP_INTER, FC_CAP_INTER), clip = "off") +
    labs(title = "Age × Region interaction (real log2 fold changes)",
         x = "Age log2 FC (Old - Young)",
         y = "Region log2 FC (Dorsal - Ventral)") +
    theme_minimal(base_size = 12) +
    theme(plot.title      = element_text(size = 14, face = "bold"),
          legend.position = "bottom")
}

main <- function() {
  out_dir <- OUT$fig3
  d <- read.csv(file.path(out_dir, "harmonized_log2fc.csv")) %>%
    mutate(
      age_significant    = age_pval    < P_THRESH & abs(log2fc_age)    > FC_THRESH,
      region_significant = region_pval < P_THRESH & abs(log2fc_region) > FC_THRESH
    )

  age <- base_volcano(
    d, x_col = "log2fc_age", p_col = "age_pval", sig_col = "age_significant",
    dir_levels = c("Old-enriched", "Young-enriched"),
    title  = "Age effects on cell-type proportions",
    x_lab  = "log2 FC (Old - Young)"
  )
  region <- base_volcano(
    d, x_col = "log2fc_region", p_col = "region_pval",
    sig_col = "region_significant",
    dir_levels = c("Dorsal-enriched", "Ventral-enriched"),
    title  = "Region effects on cell-type proportions",
    x_lab  = "log2 FC (Dorsal - Ventral)"
  )
  inter <- interaction_plot(d)

  ggsave(file.path(out_dir, "age_effects_volcano.pdf"),
         age,    width = 9, height = 7)
  ggsave(file.path(out_dir, "region_effects_volcano.pdf"),
         region, width = 9, height = 7)
  ggsave(file.path(out_dir, "interaction_effects_plot.pdf"),
         inter,  width = 9, height = 7)
  ggsave(file.path(out_dir, "age_effects_volcano.png"),
         age,    width = 9, height = 7, dpi = 300)
  ggsave(file.path(out_dir, "region_effects_volcano.png"),
         region, width = 9, height = 7, dpi = 300)
  ggsave(file.path(out_dir, "interaction_effects_plot.png"),
         inter,  width = 9, height = 7, dpi = 300)

  message("[plot-fig3] wrote age/region/interaction volcanos to ", out_dir)
}

if (sys.nframe() == 0L) main()
