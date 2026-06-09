# install_packages.R — one-time dependency installation for the Parity reproduction pipeline.
# Target: R 4.3.0 (Methods §2.9).

cran_pkgs <- c(
  "tidyverse",          # 2.0.0
  "data.table",
  "ggplot2",            # 3.4.2
  "patchwork",
  "cowplot",
  "ggrepel",
  "RColorBrewer",
  "viridis",
  "scales",
  "reshape2",
  "gridExtra",
  "pheatmap",
  "factoextra",
  "FactoMineR",
  "dendextend",
  "randomForest",       # 4.7-1.1
  "Matrix",
  "VennDiagram",
  "eulerr",
  "ggvenn"
)

bioc_pkgs <- c(
  "DESeq2",             # 1.40.0
  "clusterProfiler",    # 4.8.0
  "ReactomePA",
  "org.Mm.eg.db",
  "org.Rn.eg.db",
  "enrichplot",
  "limma",
  "edgeR",
  "sva",
  "ComplexHeatmap",
  "rhdf5",
  "Biobase"
)

github_pkgs <- c("meichendong/SCDC")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
for (p in cran_pkgs) install_if_missing(p)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE)

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
for (p in github_pkgs) {
  short <- sub(".*/", "", p)
  if (!requireNamespace(short, quietly = TRUE)) devtools::install_github(p)
}

message("All package installations complete.")
