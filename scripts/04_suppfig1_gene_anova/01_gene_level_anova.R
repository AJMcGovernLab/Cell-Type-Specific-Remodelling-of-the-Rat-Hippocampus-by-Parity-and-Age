# 01_gene_level_anova.R — gene-level three-way ANOVA (Age × Parity × Region) on 12,516 genes.
# Condensed from Set 1/Final_Results_Summary/7_Parity_Gene_Expression/ScriptsPaper/phase2_corrected_anova.R.
# Row indices 1..7 of the ANOVA summary table correspond to:
#   age, parity, region, age:parity, age:region, parity:region, age:parity:region.

suppressPackageStartupMessages({ library(dplyr); library(data.table) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

EFFECT_TERMS <- c("age", "parity", "region",
                  "age_parity", "age_region", "parity_region",
                  "age_parity_region")

# ---- single-gene ANOVA ----------------------------------------------------

anova_one_gene <- function(expr, md) {
  d <- data.frame(expression = as.numeric(expr),
                  age    = md$age, parity = md$parity, region = md$region)
  tryCatch({
    m  <- aov(expression ~ age * parity * region, data = d)
    a  <- summary(m)[[1]]                  # 8 rows: 7 terms + Residuals
    r2 <- summary(lm(expression ~ age * parity * region, data = d))$r.squared
    c(setNames(a[1:7, "Pr(>F)"],  paste0(EFFECT_TERMS, "_pval")),
      setNames(a[1:7, "F value"], paste0(EFFECT_TERMS, "_fstat")),
      model_rsquared = r2,
      residual_se    = sqrt(a[8, "Mean Sq"]))
  }, error = function(e) setNames(rep(NA_real_, 16),
     c(paste0(EFFECT_TERMS, "_pval"),
       paste0(EFFECT_TERMS, "_fstat"),
       "model_rsquared", "residual_se")))
}

# ---- main -----------------------------------------------------------------

main <- function() {
  out_dir <- OUT$suppfig1
  # Authoritative input: the paper's own pre-filtered 12,516 × 60 matrix
  # produced by phase1_data_validation_qc.R. No need to re-apply the filter.
  authoritative <- file.path(BULK_DIR, "filtered_expression_for_anova.csv")
  if (file.exists(authoritative)) {
    message(sprintf("[gene-anova] using paper's pre-filtered matrix: %s",
                    authoritative))
    expr_raw  <- fread(authoritative, data.table = FALSE)
    rownames(expr_raw) <- expr_raw[[1]]
    expr_df   <- expr_raw[, -1, drop = FALSE]
    gene_meta <- data.frame(genes = rownames(expr_df),
                            stringsAsFactors = FALSE)
    # Pull symbols from the shared gene-metadata table
    gm <- fread(GENE_META, data.table = FALSE)
    gene_meta$external_gene_name <- gm$external_gene_name[match(gene_meta$genes, gm$genes)]
  } else {
    # Fallback: re-filter normalized_filtered_counts.csv ourselves.
    expr_df   <- fread(BULK_COUNTS, data.table = FALSE)
    gene_meta <- fread(GENE_META,   data.table = FALSE)
    stopifnot(nrow(expr_df) == nrow(gene_meta))
    m <- as.matrix(expr_df)
    keep <- rowMeans(m) >= 1 & apply(m, 1, var) > 0 & rowMeans(m > 0) >= 0.1
    expr_df   <- expr_df[keep, , drop = FALSE]
    gene_meta <- gene_meta[keep, , drop = FALSE]
    rownames(expr_df) <- gene_meta$genes
  }

  md <- read.csv(SAMPLE_META, stringsAsFactors = TRUE)
  stopifnot(all(colnames(expr_df) == md$sample))

  message(sprintf("[gene-anova] running on %d genes × %d samples", nrow(expr_df), ncol(expr_df)))
  results_mat <- t(apply(as.matrix(expr_df), 1, anova_one_gene, md = md))
  results_df  <- as.data.frame(results_mat)
  results_df$gene_id    <- rownames(expr_df)
  results_df$gene_index <- seq_len(nrow(results_df))

  # FDR per effect term
  for (term in EFFECT_TERMS) {
    results_df[[paste0(term, "_fdr")]] <-
      p.adjust(results_df[[paste0(term, "_pval")]], method = "fdr")
    results_df[[paste0(term, "_significant")]] <-
      !is.na(results_df[[paste0(term, "_fdr")]]) &
      results_df[[paste0(term, "_fdr")]] < CONFIG$fdr_threshold
  }

  # Attach gene symbols
  results_df <- merge(results_df,
                      gene_meta[, c("genes", "external_gene_name")],
                      by.x = "gene_id", by.y = "genes", all.x = TRUE)

  # Consistent column order
  col_order <- c("gene_id", "external_gene_name", "gene_index",
                 unlist(lapply(EFFECT_TERMS, function(t)
                   paste0(t, c("_pval", "_fdr", "_significant", "_fstat")))),
                 "model_rsquared", "residual_se")
  results_df <- results_df[order(results_df$gene_index), col_order]

  write.csv(results_df,
            file.path(out_dir, "threeway_anova_results_CORRECTED.csv"),
            row.names = FALSE)

  # Per-effect significant-gene CSVs (ordered by p-value)
  for (term in EFFECT_TERMS) {
    sig <- results_df %>%
      filter(.data[[paste0(term, "_significant")]]) %>%
      arrange(.data[[paste0(term, "_pval")]])
    if (nrow(sig) > 0)
      write.csv(sig,
                file.path(out_dir,
                          sprintf("significant_genes_%s_CORRECTED.csv", term)),
                row.names = FALSE)
  }

  # Summary
  summary_df <- data.frame(
    effect            = c("Age","Parity","Region","Age:Parity","Age:Region","Parity:Region","Age:Parity:Region"),
    significant_genes = sapply(EFFECT_TERMS, function(t) sum(results_df[[paste0(t, "_significant")]], na.rm = TRUE)),
    percentage        = round(100 * sapply(EFFECT_TERMS, function(t)
                         sum(results_df[[paste0(t, "_significant")]], na.rm = TRUE)) / nrow(results_df), 2)
  )
  write.csv(summary_df,
            file.path(out_dir, "anova_significance_summary_CORRECTED.csv"),
            row.names = FALSE)

  message("[gene-anova] summary:"); print(summary_df)
  invisible(results_df)
}

if (sys.nframe() == 0L) main()
