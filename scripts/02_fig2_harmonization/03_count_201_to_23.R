# 03_count_201_to_23.R — counts + labels written by get_201_to_23_counts.R and check_official_harmonization.R.

suppressPackageStartupMessages(library(tidyverse))
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

main <- function() {
  out_dir <- OUT$fig2
  mapping <- read.csv(file.path(out_dir, "refined_cell_type_mapping.csv"))

  # The "total_cells" 53,382 value comes from the reference datasets' row counts.
  # It's recomputed here from the SCDC checkpoint RDS to stay self-contained.
  female_refs <- c("mouse10x_2020_female_scdc.rds",
                   "mouse_smartseq_2019_female_scdc.rds",
                   "yao_hippo_10x_female_scdc.rds")
  total_cells <- sum(vapply(female_refs, function(f) {
    r <- readRDS(file.path(SCDC_CHECKPOINTS, f)); r$n_cells_used
  }, numeric(1)))

  summary_row <- data.frame(
    total_cells       = total_cells,
    original_types    = nrow(mapping),
    harmonized_types  = n_distinct(mapping$unified_name)
  )
  write.csv(summary_row,
            file.path(out_dir, "harmonization_201_to_23_counts.csv"),
            row.names = FALSE)

  official <- data.frame(
    harmonized_type = sort(unique(mapping$unified_name)),
    stringsAsFactors = FALSE
  ) %>% mutate(type_number = row_number())
  write.csv(official,
            file.path(out_dir, "official_23_harmonized_types.csv"),
            row.names = FALSE)

  message(sprintf("[counts] %d cells, %d → %d types",
                  total_cells, nrow(mapping), n_distinct(mapping$unified_name)))
}

if (sys.nframe() == 0L) main()
