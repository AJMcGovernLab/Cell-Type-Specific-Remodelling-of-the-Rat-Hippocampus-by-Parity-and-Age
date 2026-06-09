# config.R — central path configuration for the Parity paper reproduction.
# Every script under scripts/ sources this file; no hard-coded paths elsewhere.
#
# Layout (Repository/ is the project root):
#
#   Repository/
#   ├── data/                 authoritative inputs
#   ├── checkpoints/          pre-computed intermediates
#   ├── scripts/              this folder (config.R + 9 section folders + 00_utilities)
#   ├── outputs/              produced by running run_all.R or per-section run.R
#   └── manuscript/, docs/

`%||%` <- function(a, b) if (is.null(a)) b else a

.locate_repo_root <- function() {
  # 1. explicit env var
  env <- Sys.getenv("REPO_ROOT", unset = "")
  if (nzchar(env) && dir.exists(env) &&
      file.exists(file.path(env, "scripts", "config.R"))) {
    return(normalizePath(env, winslash = "/"))
  }
  # 2. sys.frame() — when sourced interactively from within R
  sf <- try(sys.frame(1)$ofile, silent = TRUE)
  if (!inherits(sf, "try-error") && !is.null(sf)) {
    d <- dirname(normalizePath(sf, winslash = "/", mustWork = FALSE))
    # config.R lives in Repository/scripts/, so repo root is its parent
    if (file.exists(file.path(d, "config.R"))) {
      return(normalizePath(file.path(d, ".."), winslash = "/"))
    }
  }
  # 3. --file= when run via Rscript; walk up until we find scripts/config.R
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(f)) {
    d <- dirname(normalizePath(f[1], winslash = "/", mustWork = FALSE))
    for (i in 0:5) {
      cand <- do.call(file.path, c(list(d), as.list(rep("..", i))))
      cand <- normalizePath(cand, winslash = "/", mustWork = FALSE)
      if (file.exists(file.path(cand, "scripts", "config.R"))) return(cand)
    }
  }
  # 4. hardcoded fallback (this machine)
  fb <- "f:/Parity/Final/Repository"
  if (file.exists(file.path(fb, "scripts", "config.R"))) return(fb)
  stop("Could not locate Repository root; set REPO_ROOT env var.")
}

REPO_ROOT      <- .locate_repo_root()
DATA_DIR       <- file.path(REPO_ROOT, "data")
BULK_DIR       <- file.path(DATA_DIR, "bulk_rnaseq")
REF_DIR        <- file.path(DATA_DIR, "references")
CHECKPOINT_DIR <- file.path(REPO_ROOT, "checkpoints")
OUTPUT_ROOT    <- file.path(REPO_ROOT, "outputs")

BULK_COUNTS    <- file.path(BULK_DIR, "normalized_filtered_counts.csv")
SAMPLE_META    <- file.path(BULK_DIR, "sample_metadata.csv")
GENE_META      <- file.path(DATA_DIR, "gene_metadata_with_rat_names.csv")

REFS <- list(
  mouse10x_2020       = file.path(REF_DIR, "mouse10x_2020"),
  mouse_smartseq_2019 = file.path(REF_DIR, "mouse_smartseq_2019"),
  yao_hippo_10x       = file.path(REF_DIR, "yao_hippo_10x")
)

SCDC_CHECKPOINTS <- file.path(CHECKPOINT_DIR, "scdc_deconvolution")
HARM_CHECKPOINTS <- file.path(CHECKPOINT_DIR, "harmonization")

OUT <- list(
  fig1       = file.path(OUTPUT_ROOT, "01_fig1_reference_sex"),
  fig2       = file.path(OUTPUT_ROOT, "02_fig2_harmonization"),
  fig3       = file.path(OUTPUT_ROOT, "03_fig3_harmonized_anova"),
  suppfig1   = file.path(OUTPUT_ROOT, "04_suppfig1_gene_anova"),
  fig4       = file.path(OUTPUT_ROOT, "05_fig4_age_region_enrichment"),
  fig5ab     = file.path(OUTPUT_ROOT, "06_fig5ab_parity_proportions"),
  suppfig234 = file.path(OUTPUT_ROOT, "09_suppfig2to5_parity_genes"),
  fig5ce     = file.path(OUTPUT_ROOT, "10_fig5_parity_enrichment"),
  fig5fi     = file.path(OUTPUT_ROOT, "11_fig5_pathway_overlaps")
)
for (d in OUT) dir.create(d, recursive = TRUE, showWarnings = FALSE)

CONFIG <- list(
  fdr_threshold = 0.05,
  age_levels    = c("7mo", "13mo"),
  region_levels = c("dorsal", "ventral"),
  parity_levels = c("nulliparous", "primiparous", "biparous")
)

message(sprintf("[config] REPO_ROOT=%s", REPO_ROOT))
message(sprintf("[config] OUTPUT_ROOT=%s", OUTPUT_ROOT))
