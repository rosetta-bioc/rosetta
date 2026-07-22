# install.R — Set up R/Bioconductor dependencies for Rosetta
# Run once: Rscript install.R

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(c(
    "DESeq2",
    "edgeR",
    "limma",
    "clusterProfiler",
    "org.Hs.eg.db",
    "ReactomePA",
    "apeglm"
), ask = FALSE, update = FALSE)

