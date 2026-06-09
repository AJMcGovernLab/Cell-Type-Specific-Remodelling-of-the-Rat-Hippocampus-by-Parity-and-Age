# adapter.R — sourced at the top of scripts copied verbatim from the original tree.
# Redefines common absolute paths to point at the Reproduction/ layout so we don't
# have to touch every write.csv inside the legacy scripts. Legacy scripts typically
# assume a `setwd()` to their original location plus relative `"results/..."` paths;
# this adapter creates a staging dir and cd's into it.

source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

setup_legacy_dirs <- function(section_key) {
  stage <- OUT[[section_key]]
  dir.create(file.path(stage, "results"),  recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(stage, "figures"),  recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(stage, "summaries"),recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(stage, "data"),     recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(stage, "plots"),    recursive = TRUE, showWarnings = FALSE)
  setwd(stage)
  stage
}

# Path re-mapping for legacy hard-coded absolute paths.
# Used before sourcing a legacy script; R then resolves any read.csv calls
# against these mapped locations via an `on.exit` unwind.
resolve_legacy_path <- function(p) {
  repl <- c(
    "/home/ajukearth/Parity/Transfer/Final/"   = paste0(REPO_ROOT, "/"),
    "/home/ajukearth/Parity/Transfer/Parity/"  = paste0(REPO_ROOT, "/"),
    "/home/ajukearth/Parity/"                  = paste0(REPO_ROOT, "/"),
    "F:/7Parity Database Comparison/Final/"    = paste0(REPO_ROOT, "/")
  )
  for (k in names(repl)) p <- gsub(k, repl[[k]], p, fixed = TRUE)
  p
}
