# 01_harmonize_cell_types.R — core 201 → 23 harmonization via correlation-based hierarchical clustering.
# Condensed from Set 1/0Full Analysis/Step3_Harmonization/harmonize_cell_types_v2.R (339 lines).

suppressPackageStartupMessages({ library(tidyverse); library(pheatmap) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

HEIGHT_CUT <- 0.2   # paper §3.2: height cutoff on 1 - correlation distance

# ---- signature computation ------------------------------------------------

compute_signatures <- function(ref) {
  expr <- ref$expression; meta <- ref$metadata
  ct_col <- ref$cell_type_col %||% "cluster_label"
  cts <- unique(meta[[ct_col]])
  sig <- matrix(0, nrow(expr), length(cts),
                dimnames = list(rownames(expr), cts))
  for (ct in cts) {
    cells <- which(meta[[ct_col]] == ct)
    if (length(cells) > 0) sig[, ct] <- rowMeans(expr[, cells, drop = FALSE])
  }
  sig
}

# ---- name parsing ---------------------------------------------------------

parse_name <- function(name) {
  major <- case_when(
    grepl("CA1|CA2|CA3|DG|Mossy", name)                         ~ "Excitatory",
    grepl("Lamp5|Vip|Sst|Pvalb|Sncg|Serpinf1|Meis2|Ndnf", name)  ~ "Inhibitory",
    grepl("Astro", name)                                         ~ "Astrocyte",
    grepl("Oligo|MOL|NFOL|OPC|COL|MFOL", name)                   ~ "Oligodendrocyte",
    grepl("Micro|PVM|Macrophage", name)                          ~ "Microglia",
    grepl("Endo|VLMC|Peri", name)                                ~ "Vascular",
    grepl("SUB|ProS", name)                                      ~ "Subiculum",
    grepl("NP|CT", name)                                         ~ "Cortical",
    TRUE                                                         ~ "Other"
  )
  subtype <- case_when(
    grepl("CA1", name)   ~ "CA1",  grepl("CA2", name)   ~ "CA2",
    grepl("CA3", name)   ~ "CA3",  grepl("DG", name)    ~ "DG",
    grepl("Lamp5", name) ~ "Lamp5", grepl("Vip", name)  ~ "Vip",
    grepl("Sst", name)   ~ "Sst",  grepl("Pvalb", name) ~ "Pvalb",
    grepl("Sncg", name)  ~ "Sncg", grepl("SUB", name)   ~ "SUB",
    grepl("ProS", name)  ~ "ProS", TRUE                 ~ NA_character_
  )
  region <- case_when(
    grepl("-do|_do", name) ~ "dorsal",
    grepl("-ve|_ve", name) ~ "ventral",
    grepl("HPF", name)     ~ "HPF",
    TRUE                   ~ NA_character_
  )
  list(major_type = major, subtype = subtype, region = region)
}

unify_name <- function(d) {
  d %>%
    group_by(cluster_id) %>%
    mutate(unified_name = case_when(
      n_distinct(major_type) == 1 & major_type[1] == "Excitatory" &
        n_distinct(subtype) == 1 & !is.na(subtype[1]) & subtype[1] != "" ~
        paste0(subtype[1], ifelse(any(region != "" & !is.na(region)),
                                   paste0("_", first(region[region != "" & !is.na(region)])), "")),
      n_distinct(major_type) == 1 & major_type[1] == "Inhibitory" &
        n_distinct(subtype) == 1 & !is.na(subtype[1]) & subtype[1] != "" ~
        paste0(subtype[1], "_interneuron"),
      n_distinct(major_type) == 1 & major_type[1] %in% c("Astrocyte","Oligodendrocyte","Microglia") ~
        major_type[1],
      n_distinct(major_type) == 1 & major_type[1] == "Vascular"  ~ "Vascular",
      n_distinct(major_type) == 1 & major_type[1] == "Subiculum" ~ "Subiculum",
      TRUE                                                        ~ paste0("Cluster_", cluster_id)
    )) %>%
    ungroup() %>%
    mutate(
      hierarchical_level_1 = case_when(
        major_type %in% c("Excitatory","Inhibitory","Cortical","Subiculum") ~ "Neuron",
        major_type %in% c("Astrocyte","Oligodendrocyte","Microglia")        ~ "Glia",
        major_type == "Vascular"                                             ~ "Vascular",
        TRUE                                                                  ~ "Other"
      ),
      hierarchical_level_2 = major_type,
      hierarchical_level_3 = ifelse(subtype != "", subtype, major_type),
      hierarchical_level_4 = original_name
    )
}

# ---- main -----------------------------------------------------------------

main <- function() {
  out_dir <- OUT$fig2
  # Reference RDS checkpoints contain $expression + $metadata + $cell_type_col
  # produced by the deconvolution load_reference_data() function (see §01).
  ref_dir <- file.path(HARM_CHECKPOINTS, "reference_data")
  female_refs <- list(
    mouse10x_2020       = file.path(ref_dir, "ref_mouse10x_2020_female.rds"),
    mouse_smartseq_2019 = file.path(ref_dir, "ref_mouse_smartseq_2019_female.rds"),
    yao_hippo_10x       = file.path(ref_dir, "ref_yao_hippo_10x_female.rds")
  )
  ref_data <- lapply(female_refs, readRDS)

  # Signatures per reference
  message("[harm] computing signatures")
  sigs <- lapply(ref_data, compute_signatures)

  # Intersect genes
  common_genes <- Reduce(intersect, lapply(sigs, rownames))
  message(sprintf("[harm] common genes: %d", length(common_genes)))
  aligned <- lapply(sigs, function(m) m[common_genes, , drop = FALSE])
  saveRDS(aligned, file.path(out_dir, "aligned_signatures.rds"))

  # Pairwise correlation matrices (per-pair) + combined 201×201
  pair_cors <- list()
  for (i in seq_along(aligned)) for (j in seq_along(aligned))
    pair_cors[[paste(names(aligned)[i], names(aligned)[j], sep = "_vs_")]] <-
      cor(aligned[[i]], aligned[[j]], method = "pearson")
  saveRDS(pair_cors, file.path(out_dir, "correlation_matrices.rds"))

  all_sig <- do.call(cbind, aligned)
  dataset_labels <- unlist(lapply(names(aligned), function(n) rep(n, ncol(aligned[[n]]))))
  combined_cor <- cor(all_sig, method = "pearson")

  # Hierarchical clustering at h = 0.2 (1 - r)
  hc <- hclust(as.dist(1 - combined_cor), method = "average")
  clusters <- cutree(hc, h = HEIGHT_CUT)

  cluster_df <- data.frame(
    original_name = colnames(all_sig),
    dataset       = dataset_labels,
    cluster_id    = clusters,
    stringsAsFactors = FALSE
  ) %>%
    mutate(parsed = map(original_name, parse_name),
           major_type = map_chr(parsed, ~ .x$major_type),
           subtype    = map_chr(parsed, ~ .x$subtype %||% ""),
           region     = map_chr(parsed, ~ .x$region  %||% "")) %>%
    select(-parsed) %>%
    unify_name() %>%
    select(dataset, original_name, unified_name, cluster_id,
           hierarchical_level_1, hierarchical_level_2,
           hierarchical_level_3, hierarchical_level_4,
           major_type, subtype, region)

  # Only keep mapping_table.rds (consumed by downstream scripts);
  # cell_type_mapping_table.csv is superseded by refined_cell_type_mapping.csv
  # produced in 02_refine_mapping.R.
  saveRDS(cluster_df, file.path(out_dir, "mapping_table.rds"))
  saveRDS(list(cor_matrix = combined_cor,
               dataset_labels = dataset_labels,
               cell_type_names = colnames(all_sig)),
          file.path(out_dir, "combined_correlation_data.rds"))

  message(sprintf("[harm] %d original types → %d clusters (h = %.2f)",
                  nrow(cluster_df), length(unique(cluster_df$cluster_id)), HEIGHT_CUT))

  # Correlation heatmap (Fig 2b)
  pdf(file.path(out_dir, "correlation_heatmap.pdf"), width = 20, height = 20)
  pheatmap(combined_cor,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           clustering_method = "average",
           cutree_rows = length(unique(clusters)),
           cutree_cols = length(unique(clusters)),
           fontsize_row = 4, fontsize_col = 4,
           show_rownames = FALSE, show_colnames = FALSE,
           main = "Cell Type Correlation Across Datasets")
  dev.off()

  invisible(cluster_df)
}

if (sys.nframe() == 0L) main()
