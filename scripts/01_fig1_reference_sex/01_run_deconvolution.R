# 01_run_deconvolution.R — SCDC deconvolution for all 9 reference configurations.
# Condensed from Set 1/Final_Results_Summary/1_Reference_Sex_Influence/ScriptsPaper/run_deconvolution_final.R (717 → ~150 lines).

suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(Matrix)
  library(Biobase); library(rhdf5); library(SCDC)
})
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

MIN_GENES          <- 500
MIN_CELLS_PER_TYPE <- 20
MARKER_CUTOFF      <- 0.1

REFERENCES <- list(
  mouse10x_2020       = c("female", "male", "mixed"),
  mouse_smartseq_2019 = c("female", "male", "mixed"),
  yao_hippo_10x       = c("female", "male", "mixed")
)

# ---- I/O helpers ----------------------------------------------------------

load_bulk <- function() {
  counts <- fread(BULK_COUNTS, data.table = FALSE)
  genes  <- fread(GENE_META,   data.table = FALSE)
  list(counts = counts, genes = genes,
       samples = data.frame(sample_id = colnames(counts),
                            sample_clean = gsub("[^A-Za-z0-9_]", "_", colnames(counts)),
                            stringsAsFactors = FALSE))
}

load_reference <- function(dataset, sex) {
  dir <- REFS[[dataset]]
  h5  <- file.path(dir, sprintf("expression_matrix_%s.h5", sex))
  md  <- file.path(dir, sprintf("metadata_%s.csv", sex))
  if (!file.exists(h5) || !file.exists(md)) return(NULL)

  expr       <- rhdf5::h5read(h5, "/data/expression")
  gene_names <- rhdf5::h5read(h5, "/data/gene_names")
  cell_names <- rhdf5::h5read(h5, "/data/cell_names")
  if (nrow(expr) == length(cell_names) && ncol(expr) == length(gene_names)) expr <- t(expr)
  if (inherits(expr, "dgCMatrix")) expr <- as.matrix(expr)
  rownames(expr) <- gene_names; colnames(expr) <- cell_names

  meta <- fread(md, data.table = FALSE)
  ct_col <- intersect(c("cell_type","celltype","CellType","cluster_label",
                        "cell_type_designation_label","subclass_label"), colnames(meta))[1]

  list(expression = expr, metadata = meta, cell_type_col = ct_col,
       dataset = dataset, sex_type = sex)
}

# ---- Gene mapping (rat → mouse via external_gene_name) --------------------

map_genes <- function(bulk_genes, ref_genes, mapping) {
  rat2mouse <- setNames(mapping$external_gene_name, mapping$genes)
  bulk_mapped <- rat2mouse[bulk_genes]
  bulk_mapped[is.na(bulk_mapped)] <- bulk_genes[is.na(bulk_mapped)]
  common <- intersect(bulk_mapped, ref_genes)
  list(bulk_idx  = which(bulk_mapped %in% common),
       ref_idx   = which(ref_genes   %in% common),
       bulk_mapped = bulk_mapped,
       n_mapped = length(common))
}

# ---- Core: run SCDC for one (bulk, reference) pair ------------------------

run_scdc <- function(bulk, ref, mapping) {
  bulk_expr <- as.matrix(bulk$counts)
  ref_expr  <- ref$expression

  bulk_genes <- if ("external_gene_name" %in% colnames(bulk$genes)) bulk$genes$external_gene_name
                else if ("genes" %in% colnames(bulk$genes))         bulk$genes$genes
                else rownames(bulk_expr)

  gm <- map_genes(bulk_genes, rownames(ref_expr), mapping)
  if (gm$n_mapped < MIN_GENES) { message("Too few mapped genes; skipping"); return(NULL) }

  bulk_expr <- bulk_expr[gm$bulk_idx, ]
  ref_expr  <- ref_expr[gm$ref_idx,  ]
  rownames(bulk_expr) <- gm$bulk_mapped[gm$bulk_idx]
  common <- intersect(rownames(bulk_expr), rownames(ref_expr))
  bulk_expr <- bulk_expr[common, ]; ref_expr <- ref_expr[common, ]

  ct <- ref$metadata[[ref$cell_type_col]]
  keep_types <- names(which(table(ct) >= MIN_CELLS_PER_TYPE))
  if (length(keep_types) < 2) return(NULL)
  keep_cells <- ct %in% keep_types
  ref_expr <- ref_expr[, keep_cells]; ct <- ct[keep_cells]

  # SCDC expects subject pseudo-IDs
  n_subj <- min(5, floor(ncol(ref_expr) / 20))
  subj   <- paste0("subject", rep(seq_len(n_subj), length.out = ncol(ref_expr)))
  sc_eset <- ExpressionSet(
    assayData  = ref_expr,
    phenoData  = AnnotatedDataFrame(data.frame(cellType = ct, SubjectName = subj,
                                               row.names = colnames(ref_expr)))
  )
  bulk_eset <- ExpressionSet(assayData = bulk_expr)

  scdc <- SCDC::SCDC_prop(bulk_eset, sc_eset, ct.varname = "cellType",
                          sample.varname = "SubjectName",
                          ct.sub = keep_types, marker.cutoff = MARKER_CUTOFF)
  if (is.null(scdc$prop.est.mvw)) return(NULL)

  props <- as.data.frame(scdc$prop.est.mvw); props$sample <- rownames(props)
  props <- merge(props, bulk$samples, by.x = "sample", by.y = "sample_clean")

  list(proportions = props, cell_types = keep_types,
       n_genes_used = nrow(bulk_expr), n_cells_used = ncol(ref_expr),
       raw_result = scdc, reference = sprintf("%s_%s", ref$dataset, ref$sex_type),
       method = "SCDC")
}

# ---- Pipeline -------------------------------------------------------------

main <- function() {
  out_dir <- OUT$fig1; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  bulk    <- load_bulk()
  mapping <- bulk$genes[, c("genes", "external_gene_name")]
  mapping <- mapping[!is.na(mapping$genes) & !is.na(mapping$external_gene_name), ]
  mapping <- mapping[nzchar(mapping$genes) & nzchar(mapping$external_gene_name), ]

  all_results <- list(); summary_rows <- list()

  for (dataset in names(REFERENCES)) {
    for (sex in REFERENCES[[dataset]]) {
      ref_name <- sprintf("%s_%s", dataset, sex)
      message(sprintf("[scdc] %s", ref_name))

      ref <- load_reference(dataset, sex)
      if (is.null(ref)) next

      res <- run_scdc(bulk, ref, mapping)
      if (is.null(res)) next

      all_results[[ref_name]] <- res
      saveRDS(res, file.path(out_dir, sprintf("%s_scdc.rds", ref_name)))
      summary_rows[[ref_name]] <- data.frame(
        reference    = ref_name,
        n_samples    = nrow(res$proportions),
        n_cell_types = length(res$cell_types),
        n_genes      = res$n_genes_used,
        n_cells      = res$n_cells_used
      )
    }
  }

  saveRDS(all_results, file.path(out_dir, "all_results_scdc.rds"))
  write.csv(do.call(rbind, summary_rows),
            file.path(out_dir, "deconvolution_summary.csv"), row.names = FALSE)
  message(sprintf("[scdc] Wrote %d SCDC results to %s", length(all_results), out_dir))
  invisible(all_results)
}

if (sys.nframe() == 0L) main()
