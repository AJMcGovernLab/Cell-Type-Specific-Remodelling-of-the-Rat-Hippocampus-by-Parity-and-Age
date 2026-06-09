# 01_pathway_overlap_summary.R — 7-region Venn decomposition over the 4 per-cell-type
# parity enrichment meta-analysis CSVs.
#
# Writes `pathway_overlap_summary.csv` (Fig 5f-i source).
# Reads from OUT$fig5ce. If those files don't exist, falls back to the paper's
# authoritative copies in Set 1/Final_Results_Summary/8_Parity_Functional_Enrichment/.

suppressPackageStartupMessages({ library(dplyr) })
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

# Authoritative input: the v7 *significant* meta-analysis table, which already
# aggregates the four weighting methods × every database into one row per
# (pathway, cell-type group). The CA3 populations are combined upstream into
# "CA3_do_combined" so no dataset merging is needed here.
SIG_V7_REL <- file.path("functional_enrichment_v7",
                        "weighted_meta_analysis_significant_v7.csv")

main <- function() {
  out_dir <- OUT$fig5fi
  in_file <- file.path(CHECKPOINT_DIR, "enrichment_parity", SIG_V7_REL)
  message(sprintf("[overlaps] input: %s", in_file))
  d <- read.csv(in_file, stringsAsFactors = FALSE)

  ca3   <- unique(d$ID[d$meta_group == "CA3_do_combined"])
  astro <- unique(d$ID[d$meta_group == "376_Astro"])
  sst   <- unique(d$ID[d$meta_group == "78_Sst_HPF"])

  in_ca3   <- function(x) x %in% ca3
  in_astro <- function(x) x %in% astro
  in_sst   <- function(x) x %in% sst

  universe <- unique(c(ca3, astro, sst))
  regions <- tibble::tibble(pathway = universe) |>
    mutate(
      ca3 = in_ca3(pathway), astro = in_astro(pathway), sst = in_sst(pathway)
    )

  summary_df <- data.frame(
    Region = c("CA3 Dorsal only", "Astrocytes only", "Sst+ HPF only",
               "CA3 Dorsal & Astrocytes", "CA3 Dorsal & Sst+ HPF",
               "Astrocytes & Sst+ HPF", "All three cell types"),
    Count = c(
      sum(regions$ca3 & !regions$astro & !regions$sst),
      sum(!regions$ca3 & regions$astro & !regions$sst),
      sum(!regions$ca3 & !regions$astro & regions$sst),
      sum(regions$ca3 & regions$astro & !regions$sst),
      sum(regions$ca3 & !regions$astro & regions$sst),
      sum(!regions$ca3 & regions$astro & regions$sst),
      sum(regions$ca3 & regions$astro & regions$sst)
    )
  )
  # Paper normalizes %'s by nrow of the input significant table, not by the
  # 7-region sum (which ignores duplicates in intersection zones).
  summary_df$Percentage <- round(100 * summary_df$Count / nrow(d), 1)

  write.csv(summary_df,
            file.path(out_dir, "pathway_overlap_summary.csv"),
            row.names = FALSE)
  message(sprintf("[overlaps] total = %d pathways across 7 regions",
                  sum(summary_df$Count)))
  print(summary_df)

  # Also emit per-region pathway lists
  detail <- list(
    CA3_Dorsal_only            = regions$pathway[regions$ca3 & !regions$astro & !regions$sst],
    Astrocytes_only            = regions$pathway[!regions$ca3 & regions$astro & !regions$sst],
    Sst_HPF_only               = regions$pathway[!regions$ca3 & !regions$astro & regions$sst],
    CA3_Dorsal_and_Astrocytes  = regions$pathway[regions$ca3 & regions$astro & !regions$sst],
    CA3_Dorsal_and_Sst_HPF     = regions$pathway[regions$ca3 & !regions$astro & regions$sst],
    Astrocytes_and_Sst_HPF     = regions$pathway[!regions$ca3 & regions$astro & regions$sst],
    All_three_cell_types       = regions$pathway[regions$ca3 & regions$astro & regions$sst]
  )
  for (nm in names(detail)) {
    if (length(detail[[nm]]) > 0) {
      write.csv(data.frame(pathway_id = detail[[nm]]),
                file.path(out_dir, sprintf("%s_pathways.csv", nm)),
                row.names = FALSE)
    }
  }

  invisible(summary_df)
}

if (sys.nframe() == 0L) main()
