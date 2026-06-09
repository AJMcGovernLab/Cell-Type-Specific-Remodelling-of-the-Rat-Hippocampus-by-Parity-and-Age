# 03e_themed_volcanos.R — Figure 5c-e volcanos in the simple (03d) layout, but
# with dots coloured by a FIXED functional-theme taxonomy shared across all
# three cells, so the panels use ONE common legend.
#
# Direction is still readable from the x-axis (NES < 0 = down, NES > 0 = up);
# colour now encodes the pathway's functional theme instead of direction.
# Themes are assigned by keyword (classify_theme) rather than per-cell rrvgo, so
# the same colour means the same theme in every panel and the legend is shared.
#
# Outputs (in outputs/10_fig5_parity_enrichment/directional_volcanos/):
#   volcano_themed_5c_CA3_358 / 5d_Astro / 5e_Sst  -- individual panels, no legend
#   volcano_themed_combined_5cde                    -- 5c|5d|5e faceted, ONE legend

suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(ggrepel)
})

source(Sys.getenv("REPRO_CONFIG",
                  "f:/Parity/Final/Repository/scripts/config.R"))

P_THRESH    <- 0.05
NES_THRESH  <- 1.0
LABEL_TOP_N <- 5

# Fixed theme taxonomy + colours (shared across all panels).
THEME_COLOURS <- c(
  "Synaptic / projection"   = "#1b9e77",  # teal
  "Translation / ribosome"  = "#d95f02",  # orange
  "Proteostasis"            = "#7570b3",  # purple
  "Cell cycle / DNA damage" = "#e7298a",  # magenta
  "Transcription / RNA"     = "#e6ab02",  # gold
  "Metabolism / energy"     = "#66a61e",  # green   (mito + metabolism merged)
  "Ion transport / channel" = "#386cb0",  # blue
  "Vesicle / trafficking"   = "#a6761d",  # brown
  "Cytoskeleton / ECM"      = "#e31a1c",  # red     (cytoskeleton + ECM merged)
  "Immune / disease"        = "#984ea3",  # violet
  "Development"             = "#00bcd4",  # cyan (moved out of the pink/magenta family)
  "Other"                   = "#c0c0c0",  # light grey
  "Not significant"         = "#ececec"   # faint grey
)

# Keyword classifier. Underscores are normalised to spaces first so the same
# keyword matches both GO descriptions ("electron transport") and raw pathway
# IDs ("REACTOME_..._ELECTRON_TRANSPORT"). Order matters: first match wins.
classify_theme <- function(desc) {
  d <- gsub("_", " ", toupper(desc))
  dplyr::case_when(
    grepl("SYNAP|NEURON|AXON|DENDRIT|MYELIN|NERVE|JUNCTION|CELL PROJECTION|PLASMA MEMBRANE BOUNDED", d) ~ "Synaptic / projection",
    grepl("TRANSLAT|RIBOSOM|EIF[0-9]|POLYSOME|RRNA|NONSENSE MEDIATED|SRP DEPENDENT", d) ~ "Translation / ribosome",
    grepl("PROTEASOME|UBIQUITIN|DEGRAD|AUTOPHAG|UNFOLDED|MISFOLD|CHAPERON|HSP|PROTEIN FOLDING|SUMO", d) ~ "Proteostasis",
    grepl("CELL CYCLE|CHECKPOINT|G2 M|MITOTIC|GTSE|CYCLIN|DNA DAMAGE|DNA REPAIR|CENTROMERE|KINETOCHORE", d) ~ "Cell cycle / DNA damage",
    grepl("CHROMATIN|HISTONE|TRANSCRIPTION|POLYMERASE|DNA BINDING|REGULATORY REGION|SPLICEOSOM|MRNA|RNA METABOL|RNA BIOSYNTH|RNA PROCESS|NUCLEIC", d) ~ "Transcription / RNA",
    grepl("MITOCHOND|OXIDATIVE PHOSPHOR|RESPIRAT|ELECTRON TRANS|NADH|OXIDOREDUC|ATP SYNTH|TCA|CITRIC ACID|GLYCOLYSIS|AEROBIC|PROTON TRANSMEMBRANE|ENERGY DERIVATION|METABOL|BIOSYNTH|CATABOL|LIPID|STEROL|CHOLESTEROL|FATTY ACID|GLUCOSE|GLYCAN|NUCLEOTIDE", d) ~ "Metabolism / energy",
    grepl("ION CHANNEL|CHANNEL|TRANSMEMBRANE TRANSPORT|TRANSPORTER|ACTION POTENTIAL|CALCIUM|POTASSIUM|SODIUM|CATION|ANION|SARCOLEMMA|MEMBRANE POTENTIAL", d) ~ "Ion transport / channel",
    grepl("VESICLE|ENDOSOM|LYSOSOM|VACUOLE|GOLGI|ENDOCYTOS|EXOCYTOS|SECRET|TRAFFICK|MELANOSOME|PIGMENT GRANULE|ENDOPLASMIC RETICULUM", d) ~ "Vesicle / trafficking",
    grepl("CYTOSKELET|ACTIN|MICROTUBULE|TUBULIN|FILOPODIUM|CILIUM|DYNEIN|KINESIN|MYOSIN|LAMELLIPOD|SPECTRIN|EXTRACELLULAR MATRIX|COLLAGEN|INTEGRIN|ADHESION|LAMININ|BASEMENT MEMBRANE", d) ~ "Cytoskeleton / ECM",
    grepl("IMMUNE|INFLAMM|CYTOKIN|INTERFERON|VIRAL|INFECTION|SARS|SCRAPIE|PRION|SCLEROSIS|NEURODEGENER|HUNTINGTON|PARKINSON|ALZHEIM|CANCER|CARCINOMA|TUMOR|DIABET|CARDIOMYOPATHY|SALMONELLA|LEUKEMIA|MELANOMA", d) ~ "Immune / disease",
    grepl("DEVELOPMENT|MORPHOGENESIS|DIFFERENTIATION|BEHAVIOR|GLIOGENESIS|NEUROGENESIS", d) ~ "Development",
    TRUE ~ "Other"
  )
}

prep <- function(d, cell_tag) {
  nes_col  <- intersect(c("log2_fold_enrichment", "mean_nes", "weighted_mean_NES"), colnames(d))[1]
  padj_col <- intersect(c("meta_fdr", "weighted_meta_padj"), colnames(d))[1]
  stopifnot(!is.na(nes_col), !is.na(padj_col))
  d$.nes  <- d[[nes_col]]
  d$.padj <- d[[padj_col]]
  d %>% mutate(
    cell = cell_tag,
    neg_log10_p = -log10(pmax(.padj, 1e-300)),
    direction = case_when(
      .padj < P_THRESH & .nes >  NES_THRESH ~ "Upregulated",
      .padj < P_THRESH & .nes < -NES_THRESH ~ "Downregulated",
      TRUE                                  ~ "Not significant"),
    theme_raw = classify_theme(paste(ID, Description)),
    theme_plot = ifelse(direction == "Not significant", "Not significant", theme_raw),
    Description = ifelse(nchar(Description) > 45,
                         paste0(substr(Description, 1, 42), "..."), Description)
  )
}

top_labels <- function(d) {
  d %>% filter(direction != "Not significant") %>%
    group_by(cell, direction) %>%
    arrange(desc(neg_log10_p), .by_group = TRUE) %>%
    slice_head(n = LABEL_TOP_N) %>% ungroup()
}

GREY_LV <- c("Other", "Not significant")

base_layers <- function(level_order) {
  list(
    geom_vline(xintercept = 0, colour = "black", alpha = 0.3),
    geom_vline(xintercept = c(-NES_THRESH, NES_THRESH), linetype = "dashed", colour = "gray70"),
    geom_hline(yintercept = -log10(P_THRESH), linetype = "dashed", colour = "gray70"),
    # grey background points (Other / n.s.): small + transparent, drawn first
    geom_point(data = function(x) dplyr::filter(x, theme_plot %in% GREY_LV),
               aes(colour = theme_plot), size = 0.8, alpha = 0.22),
    # themed points: full size + opacity, drawn on top
    geom_point(data = function(x) dplyr::filter(x, !theme_plot %in% GREY_LV),
               aes(colour = theme_plot), size = 2.2, alpha = 0.85),
    scale_colour_manual(values = THEME_COLOURS, limits = level_order, drop = FALSE,
                        name = "Functional theme",
                        guide = guide_legend(override.aes = list(size = 3.5, alpha = 1), nrow = 2)),
    coord_cartesian(clip = "off"),
    labs(x = "NES", y = "-log10(meta-FDR)")
  )
}

repel_layer <- function(lab) {
  geom_text_repel(
    data = lab, aes(label = Description), colour = "black",
    nudge_x = ifelse(lab$direction == "Upregulated", 0.6, -0.6),
    direction = "y", hjust = ifelse(lab$direction == "Upregulated", 0, 1),
    size = 2.9, fontface = "bold", max.overlaps = Inf, seed = 42,
    box.padding = 0.6, point.padding = 0.3, min.segment.length = 0,
    segment.size = 0.3, segment.colour = "gray55"
  )
}

main <- function() {
  out_dir <- file.path(OUT$fig5ce, "directional_volcanos")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cells <- list(
    list(file = "358_CA3_do_meta_analysis.csv", tag = "5c  Dorsal CA3 pyramidal neurons", slug = "5c_CA3_358"),
    list(file = "376_Astro_meta_analysis.csv",  tag = "5d  Hippocampal astrocytes",       slug = "5d_Astro"),
    list(file = "78_Sst_HPF_meta_analysis.csv", tag = "5e  Hippocampal SST interneurons", slug = "5e_Sst")
  )

  dat <- lapply(cells, function(c)
    prep(read.csv(file.path(OUT$fig5ce, c$file), stringsAsFactors = FALSE), c$tag))

  present <- unique(unlist(lapply(dat, function(d)
    d$theme_plot[d$direction != "Not significant"])))
  level_order <- c(intersect(names(THEME_COLOURS), present), "Not significant")
  message("[03e] shared legend themes: ", paste(level_order, collapse = ", "))

  # ---- individual panels (no legend), consistent colours ------------------
  for (i in seq_along(dat)) {
    d <- dat[[i]]; d$theme_plot <- factor(d$theme_plot, levels = level_order)
    nu <- sum(d$direction == "Upregulated"); nd <- sum(d$direction == "Downregulated")
    p <- ggplot(d, aes(.nes, neg_log10_p)) + base_layers(level_order) +
      labs(title = cells[[i]]$tag, subtitle = sprintf("%d up, %d down", nu, nd)) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(size = 13, face = "bold"),
            plot.subtitle = element_text(size = 9.5, colour = "gray35"),
            axis.title = element_text(size = 11, face = "bold"),
            legend.position = "none", panel.grid.minor = element_blank(),
            plot.margin = unit(c(0.4, 1.4, 0.4, 1.4), "cm"))
    ggsave(file.path(out_dir, sprintf("volcano_themed_%s.png", cells[[i]]$slug)), p, width = 7.5, height = 6, dpi = 300)
    ggsave(file.path(out_dir, sprintf("volcano_themed_%s.pdf", cells[[i]]$slug)), p, width = 7.5, height = 6)
  }

  # ---- combined: facet_wrap = ONE shared legend ---------------------------
  big <- bind_rows(dat)
  big$cell <- factor(big$cell, levels = sapply(cells, function(c) c$tag))
  big$theme_plot <- factor(big$theme_plot, levels = level_order)
  lab <- top_labels(big); lab$cell <- factor(lab$cell, levels = levels(big$cell))

  combined <- ggplot(big, aes(.nes, neg_log10_p)) + base_layers(level_order) +
    facet_wrap(~cell, nrow = 1, scales = "free") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(size = 13, face = "bold"),
          axis.title = element_text(size = 12, face = "bold"),
          legend.position = "bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(1.4, "cm"),
          plot.margin = unit(c(0.4, 1.2, 0.4, 1.2), "cm"))

  ggsave(file.path(out_dir, "volcano_themed_combined_5cde.png"), combined, width = 20, height = 8, dpi = 300)
  ggsave(file.path(out_dir, "volcano_themed_combined_5cde.pdf"), combined, width = 20, height = 8)
  message("[03e] themed volcanos written to ", out_dir)
}

if (sys.nframe() == 0L) main()
