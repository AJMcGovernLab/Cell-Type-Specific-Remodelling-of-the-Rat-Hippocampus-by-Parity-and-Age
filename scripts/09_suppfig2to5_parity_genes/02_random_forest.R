# 02_random_forest.R — 500-tree random forest classification of parous vs nulliparous
# with elbow-based feature selection (top 23 genes from Mean Decrease in Accuracy).
# Condensed from Set 1/Final_Results_Summary/7_Parity_Gene_Expression/ScriptsPaper/
# random_forest_parity_analysis.R + create_top23_gene_plot.R.

suppressPackageStartupMessages({
  library(randomForest); library(data.table); library(dplyr)
})
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

TOP_N <- 23   # paper's elbow point
RF_INPUTS_DIR <- file.path(CHECKPOINT_DIR, "rf_inputs")

# Refit mode (DEFAULT: FALSE). When FALSE, we load the saved .RData objects and
# extract importance — bit-for-bit reproduces the paper's top-23. When TRUE, we
# re-fit with set.seed(12345); this depends on R/randomForest version and will
# NOT reproduce bit-for-bit across patch versions. See TEST_RUN_RESULTS.md.
REFIT <- identical(Sys.getenv("REPRO_RF_REFIT", "FALSE"), "TRUE")

extract_importance <- function(rf, out_dir, filename) {
  imp <- importance(rf)
  df <- data.frame(
    gene                   = rownames(imp),
    mean_decrease_accuracy = imp[, "MeanDecreaseAccuracy"],
    mean_decrease_gini     = imp[, "MeanDecreaseGini"],
    stringsAsFactors       = FALSE
  )
  df <- df[order(-df$mean_decrease_accuracy), ]
  write.csv(df, file.path(out_dir, filename), row.names = FALSE)
  df
}

main <- function() {
  out_dir <- file.path(OUT$suppfig234, "Random_Forest")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  gene_meta <- fread(GENE_META, data.table = FALSE)

  if (REFIT) {
    message("[rf] REFIT mode (set REPRO_RF_REFIT=FALSE to load saved models)")
    ml_data <- read.csv(file.path(RF_INPUTS_DIR, "ml_expression_data.csv"),
                        check.names = FALSE, row.names = 1)
    ml_meta <- read.csv(file.path(RF_INPUTS_DIR, "ml_metadata.csv"),
                        stringsAsFactors = FALSE)
    ml_data <- as.matrix(ml_data[ml_meta$sample, ])
    binary_target <- factor(ml_meta$parity_binary,
                            levels = c("Nulliparous", "Parous"))
    multi_target  <- factor(ml_meta$parity_3group,
                            levels = c("Nulliparous", "Primiparous", "Biparous"))
    mtry <- floor(sqrt(ncol(ml_data)))

    set.seed(12345)
    cv_folds <- caret::createFolds(binary_target, k = 5,
                                   list = TRUE, returnTrain = FALSE)
    binary_rf <- randomForest(x = ml_data, y = binary_target,
                              ntree = 1000, mtry = mtry, importance = TRUE)
    multi_rf  <- randomForest(x = ml_data, y = multi_target,
                              ntree = 1000, mtry = mtry, importance = TRUE)
  } else {
    message("[rf] LOADING paper's saved models (bit-for-bit reproducible)")
    e <- new.env()
    load(file.path(RF_INPUTS_DIR, "binary_random_forest.RData"),     envir = e)
    load(file.path(RF_INPUTS_DIR, "multiclass_random_forest.RData"), envir = e)
    binary_rf <- get("binary_rf",     envir = e)
    multi_rf  <- get("multiclass_rf", envir = e)
  }

  # Paper's top-23 (Supp Fig 5) is derived from the MULTICLASS model only.
  # binary RF is run for completeness but not written out.
  multi_df <- extract_importance(multi_rf, out_dir,
                                 "multiclass_feature_importance.csv")
  message(sprintf("[rf] binary OOB err     = %.3f",
                  binary_rf$err.rate[binary_rf$ntree, "OOB"]))
  message(sprintf("[rf] multiclass OOB err = %.3f",
                  multi_rf$err.rate[multi_rf$ntree, "OOB"]))

  top <- head(multi_df, TOP_N)
  top$external_gene_name <- gene_meta$external_gene_name[match(top$gene,
                                                                gene_meta$genes)]
  top$display_name <- ifelse(is.na(top$external_gene_name) |
                               !nzchar(top$external_gene_name),
                             top$gene, top$external_gene_name)
  top <- top[, c("gene", "external_gene_name", "display_name",
                 "mean_decrease_accuracy", "mean_decrease_gini")]
  write.csv(top, file.path(out_dir, "top23_genes_with_names.csv"),
            row.names = FALSE)
  message(sprintf("[rf] top-%d MDA range %.3f - %.3f", TOP_N,
                  min(top$mean_decrease_accuracy),
                  max(top$mean_decrease_accuracy)))
  invisible(list(multiclass = multi_df, top23 = top))
}

if (sys.nframe() == 0L) main()
