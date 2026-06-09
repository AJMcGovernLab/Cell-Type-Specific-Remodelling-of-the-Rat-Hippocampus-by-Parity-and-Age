# helpers.R — small shared utilities used by the section scripts.

source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

# Save a data.frame as CSV into the given section output directory.
write_out <- function(df, section_key, filename) {
  dir <- OUT[[section_key]]
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, filename)
  write.csv(df, path, row.names = FALSE)
  message(sprintf("[write] %s", path))
  invisible(path)
}

# Save an R object as RDS into the given section output directory.
save_rds_out <- function(obj, section_key, filename) {
  dir <- OUT[[section_key]]
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, filename)
  saveRDS(obj, path)
  message(sprintf("[save] %s", path))
  invisible(path)
}

# Apply Benjamini-Hochberg FDR to a numeric vector, preserving NAs.
fdr_bh <- function(p) p.adjust(p, method = "BH")

# Consistent factor-level coercion for design variables.
coerce_design_factors <- function(df) {
  df$age    <- factor(df$age,    levels = CONFIG$age_levels)
  df$region <- factor(df$region, levels = CONFIG$region_levels)
  df$parity <- factor(df$parity, levels = CONFIG$parity_levels)
  df
}
