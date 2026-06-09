# 07_overview_panels_AB.R — Figure 4 Panels A and B (overview panels).
#
# Panel A: Pathway distribution across experimental factors (bar chart).
#          Bar height = percentage of total pathway tests attributable to each
#          factor; in-bar label = total pathway count; above-bar label = number
#          of high-confidence (Moderate+) pathways and their % of that factor.
# Panel B: Cell-type vulnerability ranking (horizontal bar chart).
#          Bars are cell-type × factor combinations sorted by vulnerability_score
#          (= mean_meta_score × (high_conf_total + 1)) and coloured by
#          factor. Labels show high-confidence pathway count per bar.

suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2)
})
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

# Factor colours — match the original Figure 4 colour scheme
factor_cols <- c(
  "Age"        = "#E31A1C",
  "Region"     = "#1F78B4",
  "Age×Region" = "#33A02C"
)

# Effect-key → display-name
factor_display <- c("age" = "Age", "region" = "Region", "age_region" = "Age×Region")

HIGH_CONF_TIERS <- c("Ultra-High", "High", "Moderate")

main <- function() {
  out_dir <- OUT$fig4
  meta_results <- read.csv(file.path(REPO_ROOT, "checkpoints", "enrichment_fig4",
                                     "all_meta_analysis_results.csv"))
  vuln <- read.csv(file.path(out_dir, "cell_type_vulnerability_ranking.csv"))

  # ----------------------------------------------------------------------
  # Panel A — Pathway distribution
  # ----------------------------------------------------------------------
  total_tests <- nrow(meta_results)

  panel_a_data <- meta_results %>%
    mutate(factor_display = factor_display[effect]) %>%
    group_by(factor_display) %>%
    summarise(
      total_pathways = n(),
      high_conf      = sum(confidence_tier %in% HIGH_CONF_TIERS),
      .groups        = "drop"
    ) %>%
    mutate(
      percentage    = 100 * total_pathways / total_tests,
      high_conf_pct = 100 * high_conf / total_pathways,
      factor_display = factor(factor_display, levels = names(factor_cols))
    )

  panel_a <- ggplot(panel_a_data,
                    aes(x = factor_display, y = percentage,
                        fill = factor_display)) +
    geom_col(alpha = 0.85, colour = "white", linewidth = 0.5) +
    geom_text(aes(label = paste0(format(total_pathways, big.mark = ","),
                                 "\npathways")),
              vjust = 1.5, colour = "black", size = 4, fontface = "bold") +
    geom_text(aes(label = sprintf("%d high-conf\n(%.3f%%)",
                                  high_conf, high_conf_pct),
                  y = percentage + 2),
              vjust = 0, size = 3.5, fontface = "bold", colour = "black") +
    scale_fill_manual(values = factor_cols, name = "Factor") +
    scale_y_continuous(limits = c(0, max(panel_a_data$percentage) + 10),
                       expand = c(0, 0)) +
    labs(title = "Pathway distribution across experimental factors",
         subtitle = "High-confidence counts shown above each bar",
         x = "Experimental factor",
         y = "Percentage of total pathway tests (%)") +
    theme_minimal(base_size = 12) +
    theme(text             = element_text(colour = "black"),
          plot.title       = element_text(size = 14, face = "bold",
                                          colour = "black"),
          plot.subtitle    = element_text(size = 11, colour = "black"),
          axis.title       = element_text(size = 11, face = "bold",
                                          colour = "black"),
          axis.text        = element_text(size = 11, colour = "black"),
          legend.title     = element_text(colour = "black"),
          legend.text      = element_text(colour = "black"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position  = "bottom")

  ggsave(file.path(out_dir, "Figure4_PanelA_pathway_distribution.pdf"),
         panel_a, width = 7, height = 6)
  ggsave(file.path(out_dir, "Figure4_PanelA_pathway_distribution.png"),
         panel_a, width = 7, height = 6, dpi = 300)

  message(sprintf("[fig4-A] %d total measurements; high-conf per factor: %s",
                  total_tests,
                  paste(sprintf("%s=%d", panel_a_data$factor_display,
                                panel_a_data$high_conf), collapse = ", ")))

  # ----------------------------------------------------------------------
  # Panel B — Cell-type vulnerability ranking
  # ----------------------------------------------------------------------
  panel_b_data <- vuln %>%
    mutate(factor_display = factor_display[effect],
           cell_factor    = paste0(cell_type, " - ", factor_display)) %>%
    arrange(desc(vulnerability_score)) %>%
    # Drop pure-zero entries (n=0 high-confidence, vulnerability_score=0); they
    # clutter the bottom of the plot without adding information.
    filter(vulnerability_score > 0) %>%
    mutate(cell_factor    = factor(cell_factor, levels = rev(cell_factor)),
           factor_display = factor(factor_display, levels = names(factor_cols)))

  panel_b <- ggplot(panel_b_data,
                    aes(x = cell_factor, y = vulnerability_score,
                        fill = factor_display)) +
    geom_col(alpha = 0.85, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%d pathway%s",
                                  high_conf_total,
                                  ifelse(high_conf_total == 1, "", "s"))),
              hjust = -0.1, size = 3.2, colour = "black") +
    scale_fill_manual(values = factor_cols, name = "Factor") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
    coord_flip() +
    labs(title = "Cell-type vulnerability ranking",
         subtitle = "Vulnerability = mean meta-score × (high-confidence pathway count + 1)",
         x = NULL,
         y = "Vulnerability score") +
    theme_minimal(base_size = 12) +
    theme(text             = element_text(colour = "black"),
          plot.title       = element_text(size = 14, face = "bold",
                                          colour = "black"),
          plot.subtitle    = element_text(size = 11, colour = "black"),
          axis.title       = element_text(size = 11, face = "bold",
                                          colour = "black"),
          axis.text.x      = element_text(size = 11, colour = "black"),
          axis.text.y      = element_text(size = 11, colour = "black"),
          legend.title     = element_text(colour = "black"),
          legend.text      = element_text(colour = "black"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position  = "bottom")

  ggsave(file.path(out_dir, "Figure4_PanelB_vulnerability_ranking.pdf"),
         panel_b, width = 9, height = max(4, nrow(panel_b_data) * 0.45 + 2))
  ggsave(file.path(out_dir, "Figure4_PanelB_vulnerability_ranking.png"),
         panel_b, width = 9, height = max(4, nrow(panel_b_data) * 0.45 + 2),
         dpi = 300)

  message(sprintf("[fig4-B] %d cell-type × factor entries plotted",
                  nrow(panel_b_data)))
}

if (sys.nframe() == 0L) main()
