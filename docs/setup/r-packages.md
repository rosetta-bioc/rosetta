# R Packages

Rosetta wraps R/Bioconductor packages. You only need to install the R packages for the wrappers you use.

---

## Package matrix

| rosetta function | R packages required | Install command |
|-----------------|--------------------|-----------------| 
| `deseq2()`, `run_deseq2()`, `get_results()` | DESeq2 | `BiocManager::install("DESeq2")` |
| `lfc_shrink(type="apeglm")` | DESeq2, apeglm | `BiocManager::install(c("DESeq2", "apeglm"))` |
| `lfc_shrink(type="ashr")` | DESeq2, ashr | `BiocManager::install(c("DESeq2", "ashr"))` |
| `edger()` | edgeR | `BiocManager::install("edgeR")` |
| `limma_voom()` | limma, edgeR | `BiocManager::install(c("limma", "edgeR"))` |
| `ORA.enrich_go()` | clusterProfiler, org.*.eg.db | `BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))` |
| `ORA.enrich_kegg()` | clusterProfiler | `BiocManager::install("clusterProfiler")` |
| `ORA.enrich_pathway()` | ReactomePA | `BiocManager::install("ReactomePA")` |
| `GSEA.gse_go()` | clusterProfiler, org.*.eg.db | `BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))` |
| `GSEA.gse_kegg()` | clusterProfiler | `BiocManager::install("clusterProfiler")` |
| `Seurat()` | Seurat, SeuratObject | `install.packages("Seurat")` |
| `Phyloseq()` | phyloseq | `BiocManager::install("phyloseq")` |

---

## Install everything

To install all supported R packages at once:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "DESeq2",
    "edgeR",
    "limma",
    "clusterProfiler",
    "org.Hs.eg.db",
    "apeglm",
    "ashr",
    "phyloseq",
    "ReactomePA"
))

install.packages("Seurat")
```

Or use the included install script:

```bash
Rscript install.R
```

---

## Minimal install (DE only)

For just differential expression (DESeq2 + edgeR + limma):

```r
BiocManager::install(c("DESeq2", "edgeR", "limma"))
```

This is sufficient for `deseq2()`, `edger()`, `limma_voom()`, and the Quick API equivalents.

---

## Organism annotation databases

For GO enrichment analysis, you need the OrgDb package for your organism:

| Organism | Package | Install |
|----------|---------|---------|
| Human | `org.Hs.eg.db` | `BiocManager::install("org.Hs.eg.db")` |
| Mouse | `org.Mm.eg.db` | `BiocManager::install("org.Mm.eg.db")` |
| Rat | `org.Rn.eg.db` | `BiocManager::install("org.Rn.eg.db")` |
| Drosophila | `org.Dm.eg.db` | `BiocManager::install("org.Dm.eg.db")` |
| Zebrafish | `org.Dr.eg.db` | `BiocManager::install("org.Dr.eg.db")` |
| Yeast | `org.Sc.sgd.db` | `BiocManager::install("org.Sc.sgd.db")` |
| Arabidopsis | `org.At.tair.db` | `BiocManager::install("org.At.tair.db")` |

---

## Checking what's installed

From Python:

```python
from rosetta._deps import is_installed

packages = ["DESeq2", "edgeR", "limma", "clusterProfiler", 
            "Seurat", "phyloseq", "apeglm", "ashr"]

for pkg in packages:
    status = "✓" if is_installed(pkg) else "✗"
    print(f"  {status} {pkg}")
```

From R:

```r
packages <- c("DESeq2", "edgeR", "limma", "clusterProfiler", 
              "Seurat", "phyloseq", "apeglm", "ashr")

for (pkg in packages) {
    installed <- requireNamespace(pkg, quietly = TRUE)
    cat(sprintf("  %s %s\n", ifelse(installed, "✓", "✗"), pkg))
}
```

---

## Version requirements

| R Package | Minimum version | Notes |
|-----------|----------------|-------|
| DESeq2 | 1.30+ | For `lfcShrink` with apeglm support |
| edgeR | 3.32+ | For `voomLmFit` (edgeR v4) |
| limma | 3.46+ | For `voomLmFit` integration |
| Seurat | 5.0+ | Recommended for latest pipeline |
| clusterProfiler | 4.0+ | For current API |
| phyloseq | 1.34+ | Stable API |

Rosetta targets current Bioconductor releases. If you're using an older R/Bioconductor version, some features may not be available.
