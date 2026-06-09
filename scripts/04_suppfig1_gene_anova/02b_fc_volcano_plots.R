# 02b_fc_volcano_plots.R — Supplementary Figure 1.
#
# Age / Region / Parity: log2 fold-change volcanos. Points coloured by:
#   dark red   = FDR<0.05 and log2FC >  +0.3  (up, past FC barrier)
#   dark blue  = FDR<0.05 and log2FC <  -0.3  (down, past FC barrier)
#   orange     = FDR<0.05 but |log2FC| <= 0.3 (significant, not past FC barrier)
#   grey       = not significant
#
# Age×Region: an interaction scatter — log2FC on each axis (age effect in dorsal
# vs age effect in ventral), coloured by the interaction p-value. Points off the
# y=x diagonal have region-dependent age effects (the interaction).
#
# Per-gene log2FC is computed from the linear normalized-count matrix the ANOVA
# used; the p-value/FDR is the corresponding three-way ANOVA term.

suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(ggplot2); library(ggrepel)
})
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

FDR_THRESH <- 0.05
FC_THRESH  <- 0.3
LABEL_N    <- 10

DIR_COLOURS <- c("Up"           = "#a50f15",   # dark red
                 "Down"         = "#08519c",   # dark blue
                 "Sig (low FC)" = "#ff7f00",   # orange
                 "n.s."         = "#cccccc")
DIR_LEVELS  <- names(DIR_COLOURS)

main <- function() {
  out_dir <- OUT$suppfig1
  expr <- fread(file.path(BULK_DIR, "filtered_expression_for_anova.csv"), data.table = FALSE)
  m <- as.matrix(expr[, -1]); rownames(m) <- expr[[1]]
  md <- read.csv(file.path(BULK_DIR, "aligned_sample_metadata.csv"), stringsAsFactors = FALSE)
  stopifnot(all(colnames(m) == md$sample))

  gm <- function(sel) rowMeans(m[, sel, drop = FALSE])
  l2 <- function(a, b) log2((a + 1) / (b + 1))
  old <- md$age == "Old"; yng <- md$age == "Young"
  dor <- md$region == "Dorsal"; ven <- md$region == "Ventral"
  par <- md$parity %in% c("Primiparous", "Biparous"); nul <- md$parity == "Nulliparous"

  fc <- data.frame(
    gene_id     = rownames(m),
    Age         = l2(gm(old), gm(yng)),
    Region      = l2(gm(dor), gm(ven)),
    Parity      = l2(gm(par), gm(nul)),
    age_dorsal  = l2(gm(old & dor), gm(yng & dor)),   # age effect within dorsal
    age_ventral = l2(gm(old & ven), gm(yng & ven)),   # age effect within ventral
    age_nullip  = l2(gm(old & nul), gm(yng & nul)),   # age effect within nulliparous
    age_parous  = l2(gm(old & par), gm(yng & par))    # age effect within parous
  )
  an <- read.csv(file.path(out_dir, "threeway_anova_results_CORRECTED.csv"), stringsAsFactors = FALSE)
  an$label <- ifelse(is.na(an$external_gene_name) | an$external_gene_name == "",
                     an$gene_id, an$external_gene_name)

  # ---- Age / Region / Parity main-effect volcanos -------------------------
  main_volcano <- function(fc_col, pcol, fdrcol, ttl, file) {
    d <- an %>%
      transmute(gene_id, label, pval = .data[[pcol]], fdr = .data[[fdrcol]]) %>%
      left_join(fc %>% select(gene_id, log2fc = all_of(fc_col)), by = "gene_id") %>%
      mutate(neg_log10p = -log10(pmax(pval, 1e-300)),
             direction = case_when(
               fdr < FDR_THRESH & log2fc >  FC_THRESH ~ "Up",
               fdr < FDR_THRESH & log2fc < -FC_THRESH ~ "Down",
               fdr < FDR_THRESH                        ~ "Sig (low FC)",
               TRUE                                    ~ "n.s."))
    d$direction <- factor(d$direction, levels = DIR_LEVELS)
    lab <- d %>% filter(direction %in% c("Up", "Down")) %>%
      arrange(desc(neg_log10p)) %>% slice_head(n = LABEL_N)
    p <- ggplot(d, aes(log2fc, neg_log10p, colour = direction)) +
      geom_vline(xintercept = c(-FC_THRESH, FC_THRESH), linetype = "dashed", colour = "gray70") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "gray70") +
      geom_point(alpha = 0.7, size = 1.5) +
      geom_text_repel(data = lab, aes(label = label), size = 3, max.overlaps = 20,
                      seed = 42, colour = "black", show.legend = FALSE) +
      scale_colour_manual(values = DIR_COLOURS, drop = FALSE, name = NULL) +
      labs(title = ttl, x = "log2 fold change", y = "-log10(p-value)") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(size = 13, face = "bold", colour = "black"),
            legend.position = "bottom", panel.grid.minor = element_blank())
    ggsave(file.path(out_dir, paste0(file, ".png")), p, width = 6, height = 6, dpi = 600)
    message(sprintf("[volcano] %s -> %s.png", fc_col, file))
  }
  main_volcano("Age",    "age_pval",    "age_fdr",    "Age (Old vs Young)",             "supp_fig1_age_volcano")
  main_volcano("Region", "region_pval", "region_fdr", "Region (Dorsal vs Ventral)",     "supp_fig1_region_volcano")
  main_volcano("Parity", "parity_pval", "parity_fdr", "Parity (Parous vs Nulliparous)", "supp_fig1_parity_volcano")

  # ---- interaction scatters: FC vs FC on each axis, coloured by p ---------
  interaction_scatter <- function(xcol, ycol, pcol, fdrcol, xlab, ylab, ttl, sub, file) {
    di <- an %>%
      transmute(gene_id, label, pval = .data[[pcol]], fdr = .data[[fdrcol]]) %>%
      left_join(fc %>% select(gene_id, x = all_of(xcol), y = all_of(ycol)), by = "gene_id") %>%
      mutate(neg_log10p = -log10(pmax(pval, 1e-300)))
    rng <- max(abs(c(di$x, di$y)), na.rm = TRUE)
    labi <- di %>% filter(fdr < FDR_THRESH) %>% arrange(desc(neg_log10p)) %>% slice_head(n = LABEL_N)
    p <- ggplot(di %>% arrange(neg_log10p), aes(x, y, colour = neg_log10p)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "gray60") +
      geom_hline(yintercept = 0, colour = "gray85") + geom_vline(xintercept = 0, colour = "gray85") +
      geom_point(alpha = 0.8, size = 1.6) +
      geom_text_repel(data = labi, aes(label = label), size = 3, max.overlaps = 20,
                      seed = 42, colour = "black", show.legend = FALSE) +
      scale_colour_viridis_c(option = "C", name = "-log10(p)", direction = -1) +
      coord_equal(xlim = c(-rng, rng), ylim = c(-rng, rng)) +
      labs(title = ttl, subtitle = sub, x = xlab, y = ylab) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(size = 13, face = "bold", colour = "black"),
            plot.subtitle = element_text(size = 9, colour = "gray35"),
            panel.grid.minor = element_blank(), legend.position = "right")
    ggsave(file.path(out_dir, paste0(file, ".png")), p, width = 6.8, height = 6, dpi = 600)
    message(sprintf("[interaction] %s -> %s.png", ttl, file))
  }
  interaction_scatter("age_dorsal", "age_ventral", "age_region_pval", "age_region_fdr",
                      "Age effect in dorsal (log2FC)", "Age effect in ventral (log2FC)",
                      "Age × Region interaction",
                      "Age effect (log2FC, Old vs Young) within each region; off-diagonal = interaction",
                      "supp_fig1_age_region_volcano")
  interaction_scatter("age_nullip", "age_parous", "age_parity_pval", "age_parity_fdr",
                      "Age effect in nulliparous (log2FC)", "Age effect in parous (log2FC)",
                      "Age × Parity interaction",
                      "Age effect (log2FC, Old vs Young) within each parity group; off-diagonal = interaction",
                      "supp_fig1_age_parity_volcano")
}

if (sys.nframe() == 0L) main()
