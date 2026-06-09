# 04_rf_pc_venn.R — Supp Fig 2: 4-way Venn of parity genes from
# RF top-23 ∪ PC5 ∪ PC6 ∪ PC8 driver gene sets.

suppressPackageStartupMessages({
  library(VennDiagram); library(grid); library(dplyr)
})
source(Sys.getenv("REPRO_CONFIG", "f:/Parity/Final/Repository/scripts/config.R"))

main <- function() {
  out_dir <- OUT$suppfig234
  rf_file <- file.path(out_dir, "Random_Forest", "top23_genes_with_names.csv")
  pc_dir  <- file.path(out_dir, "PC_Drivers")

  rf_df  <- read.csv(rf_file, stringsAsFactors = FALSE)
  rf  <- rf_df[[intersect(c("gene_id", "gene", "ensembl_id"), names(rf_df))[1]]]
  read_pc <- function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    d[[intersect(c("gene_id", "gene", "ensembl_id"), names(d))[1]]]
  }
  pc5 <- read_pc(file.path(pc_dir, "PC5_Parity_Binary_improved_drivers.csv"))
  pc6 <- read_pc(file.path(pc_dir, "PC6_Parity_Binary_improved_drivers.csv"))
  pc8 <- read_pc(file.path(pc_dir, "PC8_Parity_Binary_improved_drivers.csv"))

  sets <- list(`RF top 23` = rf, PC5 = pc5, PC6 = pc6, PC8 = pc8)
  overlap_sizes <- data.frame(
    set       = names(sets),
    n_genes   = lengths(sets)
  )
  write.csv(overlap_sizes,
            file.path(out_dir, "rf_vs_pc_parity_set_sizes.csv"),
            row.names = FALSE)

  # 3-way + 4-way intersection counts
  three_way <- length(Reduce(intersect, list(pc5, pc6, pc8)))
  four_way  <- length(Reduce(intersect, list(rf, pc5, pc6, pc8)))
  message(sprintf("[venn] |RF|=%d  |PC5|=%d  |PC6|=%d  |PC8|=%d",
                  length(rf), length(pc5), length(pc6), length(pc8)))
  message(sprintf("[venn] PC5∩PC6∩PC8=%d  all-four=%d",
                  three_way, four_way))

  out_png <- file.path(out_dir, "supp_fig2_rf_pc_venn.png")
  grid.newpage()
  venn.plot <- venn.diagram(sets,
                            filename = NULL,
                            fill = c("#EF476F", "#06D6A0", "#118AB2", "#FFD166"),
                            alpha = 0.6,
                            cex = 1.2, cat.cex = 1.2,
                            main = "Supp Fig 2 | RF vs PC5/PC6/PC8 parity drivers")
  png(out_png, width = 6, height = 6, units = "in", res = 600)
  grid.draw(venn.plot)
  dev.off()
  message(sprintf("[venn] wrote %s", out_png))
}

if (sys.nframe() == 0L) main()
