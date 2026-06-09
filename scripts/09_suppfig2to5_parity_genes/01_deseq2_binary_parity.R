# 01_deseq2_binary_parity.R — DESeq2 binary parity (parous vs nulliparous) + age covariate.
# Condensed from Set 1/Final_Results_Summary/7_Parity_Gene_Expression/ScriptsPaper/
# phase2_corrected_anova.R (DESeq2 block).

suppressPackageStartupMessages({ library(DESeq2); library(data.table); library(dplyr) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

main <- function() {
  out_dir <- file.path(OUT$suppfig234, "DESeq2")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expr_df   <- fread(BULK_COUNTS, data.table = FALSE)
  gene_meta <- fread(GENE_META,   data.table = FALSE)
  md        <- read.csv(SAMPLE_META, stringsAsFactors = FALSE)

  stopifnot(nrow(expr_df) == nrow(gene_meta))
  stopifnot(all(colnames(expr_df) == md$sample))

  m <- as.matrix(expr_df)
  keep <- rowMeans(m) >= 1 & apply(m, 1, var) > 0 & rowMeans(m > 0) >= 0.1
  m    <- m[keep, , drop = FALSE]
  gene_meta <- gene_meta[keep, , drop = FALSE]
  rownames(m) <- gene_meta$genes

  # Counts CSV is normalized (library-scaled) floats; DESeq2 needs integer counts.
  counts <- round(m); storage.mode(counts) <- "integer"

  md$parity_binary <- factor(ifelse(md$parity == "Nulliparous",
                                    "Nulliparous", "Parous"),
                             levels = c("Nulliparous", "Parous"))
  md$age <- factor(md$age, levels = c("Young", "Old"))

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData   = md,
                                design    = ~ age + parity_binary)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, name = "parity_binary_Parous_vs_Nulliparous")

  out <- as.data.frame(res)
  out$gene_id            <- rownames(out)
  out$external_gene_name <- gene_meta$external_gene_name[match(out$gene_id,
                                                                gene_meta$genes)]
  out$significant        <- !is.na(out$padj)   & out$padj < 0.05
  out$significant_fdr20  <- !is.na(out$padj)   & out$padj < 0.20
  out$large_effect       <- !is.na(out$log2FoldChange) & abs(out$log2FoldChange) > 0.3
  out$weighted_score     <- 0.5 * abs(out$log2FoldChange) + 0.5 * -log10(out$padj)
  out$weighted_score[!is.finite(out$weighted_score)] <- 0

  col_order <- c("gene_id", "external_gene_name", "baseMean",
                 "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
                 "significant", "significant_fdr20", "large_effect",
                 "weighted_score")
  out <- out[order(out$pvalue), col_order]

  # Significant subset (FDR < 0.05) kept for backward compatibility.
  write.csv(out[out$significant, ],
            file.path(out_dir, "significant_genes_parity.csv"),
            row.names = FALSE)
  # Full results table used by Supp Fig 3 volcano plot.
  write.csv(out,
            file.path(out_dir, "all_genes_parity_deseq2.csv"),
            row.names = FALSE)

  message(sprintf("[deseq2] total=%d  p<.05=%d  FDR<.05=%d  |FC|>.3=%d",
                  nrow(out),
                  sum(out$pvalue < 0.05, na.rm = TRUE),
                  sum(out$significant),
                  sum(out$large_effect, na.rm = TRUE)))
  invisible(out)
}

if (sys.nframe() == 0L) main()
