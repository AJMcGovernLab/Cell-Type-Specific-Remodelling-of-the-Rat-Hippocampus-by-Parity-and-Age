# 02_refine_mapping.R — refined biological naming with hippocampal field / IN-subtype / spatial columns.
# Condensed from Set 1/0Full Analysis/Step3_Harmonization/refine_harmonization.R.
#
# NOTE: The refined clustering threshold (0.90 correlation, h = 0.10) differs from the
# main harmonization threshold (0.80, h = 0.20) used by 01_harmonize_cell_types.R.
# The paper's §3.2 text cites h = 0.2; the `refined_cell_type_mapping.csv` actually
# uses h = 0.1. This matches Set 1's historical file. Flagged in RESULTS_TO_FILES.md.

suppressPackageStartupMessages(library(tidyverse))
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

REFINE_CORR <- 0.90

enhanced_parse <- function(name) {
  data.frame(
    cell_class = case_when(
      grepl("CA[1-3]|DG|Mossy", name)            ~ "Pyramidal",
      grepl("Lamp5|Lhx6|Vip", name)              ~ "CGE_Interneuron",
      grepl("Sst|Chodl|Pvalb", name)             ~ "MGE_Interneuron",
      grepl("Sncg|Serpinf1|Meis2", name)         ~ "Hippocampal_Interneuron",
      grepl("Ndnf", name)                        ~ "Neurogliaform",
      grepl("Astro", name)                       ~ "Astrocyte",
      grepl("Oligo|OPC|COL|MFOL|MOL|NFOL", name) ~ "Oligodendrocyte",
      grepl("Micro|PVM|Macrophage", name)        ~ "Microglia",
      grepl("Endo|VLMC|SMC|Peri", name)          ~ "Vascular",
      grepl("SUB|ProS", name)                    ~ "Subiculum",
      grepl("L[1-6]|IT|PT|CT|NP", name)          ~ "Cortical",
      grepl("CR", name)                          ~ "Cajal_Retzius",
      TRUE                                       ~ "Other"
    ),
    hippocampal_field = case_when(
      grepl("CA1", name)    ~ "CA1", grepl("CA2", name) ~ "CA2",
      grepl("CA3", name)    ~ "CA3", grepl("DG", name)  ~ "DG",
      grepl("SUB|ProS", name) ~ "SUB", TRUE            ~ NA_character_),
    interneuron_subtype = case_when(
      grepl("Lamp5", name) ~ "Lamp5", grepl("Vip", name)   ~ "Vip",
      grepl("Sst", name)   ~ "Sst",   grepl("Pvalb", name) ~ "Pvalb",
      grepl("Sncg", name)  ~ "Sncg",  grepl("Ndnf", name)  ~ "Ndnf",
      grepl("Chodl", name) ~ "Sst_Chodl", TRUE            ~ NA_character_),
    spatial_location = case_when(
      grepl("-do|_do", name) ~ "dorsal",
      grepl("-ve|_ve", name) ~ "ventral",
      grepl("HPF", name)     ~ "HPF",  TRUE ~ NA_character_),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  out_dir <- OUT$fig2

  # Default behaviour: copy the paper's authoritative refined mapping into
  # outputs/02_fig2_harmonization/ so downstream (§03) reproduces paper cell-type labels.
  # Set REPRO_REFINE_FROM_SCRATCH=TRUE to re-cluster from aligned signatures.
  refit <- identical(Sys.getenv("REPRO_REFINE_FROM_SCRATCH", "FALSE"), "TRUE")
  if (!refit) {
    src_csv <- file.path(HARM_CHECKPOINTS, "refined_cell_type_mapping.csv")
    src_rds <- file.path(HARM_CHECKPOINTS, "refined_mapping.rds")
    if (file.exists(src_csv) && file.exists(src_rds)) {
      file.copy(src_csv, file.path(out_dir, "refined_cell_type_mapping.csv"),
                overwrite = TRUE)
      file.copy(src_rds, file.path(out_dir, "refined_mapping.rds"),
                overwrite = TRUE)
      message("[refine] copied paper's authoritative refined mapping from checkpoints")
      return(invisible(read.csv(file.path(out_dir, "refined_cell_type_mapping.csv"))))
    }
  }

  aligned <- readRDS(file.path(out_dir, "aligned_signatures.rds"))

  all_sig <- do.call(cbind, aligned)
  dataset_labels <- unlist(lapply(names(aligned), function(n) rep(n, ncol(aligned[[n]]))))
  combined_cor <- cor(all_sig, method = "pearson")

  clusters <- cutree(hclust(as.dist(1 - combined_cor), method = "average"),
                     h = 1 - REFINE_CORR)

  mapping <- data.frame(
    dataset       = dataset_labels,
    original_name = colnames(all_sig),
    cluster_id    = clusters,
    stringsAsFactors = FALSE
  ) %>% bind_cols(enhanced_parse(.$original_name))

  # Unified biological names
  mapping <- mapping %>%
    group_by(cluster_id) %>%
    mutate(unified_name = case_when(
      any(hippocampal_field == "CA3", na.rm = TRUE) & any(spatial_location == "dorsal", na.rm = TRUE) ~ "CA3_dorsal",
      any(hippocampal_field == "CA3", na.rm = TRUE) & any(spatial_location == "ventral", na.rm = TRUE) ~ "CA3_ventral",
      any(hippocampal_field == "CA1", na.rm = TRUE) & any(spatial_location == "dorsal", na.rm = TRUE) ~ "CA1/CA2_dorsal",
      any(hippocampal_field == "DG",  na.rm = TRUE) & any(spatial_location == "dorsal", na.rm = TRUE) ~ "DG/CA1/CA2/CA3_dorsal",
      any(hippocampal_field == "DG",  na.rm = TRUE)                                                    ~ "DG",
      cell_class[1] == "MGE_Interneuron" & any(interneuron_subtype == "Sst", na.rm = TRUE)             ~ "Sst_IN",
      cell_class[1] == "MGE_Interneuron" & any(interneuron_subtype == "Pvalb", na.rm = TRUE)           ~ "Pvalb/Vip_IN",
      cell_class[1] == "CGE_Interneuron" & any(interneuron_subtype == "Vip", na.rm = TRUE)             ~ "Vip_IN",
      cell_class[1] == "CGE_Interneuron" & any(interneuron_subtype == "Lamp5", na.rm = TRUE)           ~ "Lamp5_IN",
      cell_class[1] == "Hippocampal_Interneuron" & any(interneuron_subtype == "Sncg", na.rm = TRUE)    ~ "Sncg_IN",
      cell_class[1] == "Astrocyte"                                                                      ~ "Astrocyte",
      cell_class[1] == "Oligodendrocyte"                                                                ~ "Oligodendrocyte",
      cell_class[1] == "Microglia"                                                                      ~ "Microglia",
      cell_class[1] == "Subiculum"                                                                      ~ "Subiculum",
      TRUE                                                                                              ~ paste0("Cluster_", cluster_id)
    )) %>%
    ungroup() %>%
    select(dataset, original_name, unified_name, cluster_id,
           cell_class, hippocampal_field, interneuron_subtype, spatial_location)

  write.csv(mapping, file.path(out_dir, "refined_cell_type_mapping.csv"), row.names = FALSE)
  saveRDS(mapping, file.path(out_dir, "refined_mapping.rds"))
  message(sprintf("[refine] %d rows; %d unified types (threshold r >= %.2f)",
                  nrow(mapping), n_distinct(mapping$unified_name), REFINE_CORR))
  invisible(mapping)
}

if (sys.nframe() == 0L) main()
