# Posit Cloud Setup

Rosetta works in Posit Cloud (formerly RStudio Cloud) with minimal setup.

## Quick Start

In a Posit Cloud Python project or Quarto notebook:

```bash
pip install rosetta-bioc
```

Then create `install.R` to ensure Bioconductor packages are available:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "edgeR", "limma", "clusterProfiler"), ask = FALSE)
```

Run once:
```bash
Rscript install.R
```

## Usage in Quarto

```{python}
import rosetta as rb

# DESeq2 differential expression
results = rb.run_deseq2(count_matrix, col_data, design="~ condition")
sig = rb.get_results(results, alpha=0.05)
```

## Notes

- Posit Cloud uses Ubuntu + R 4.x — `rpy2` works out of the box
- Bioconductor packages install to the project library (persists across sessions)
- For large datasets, use a Posit Cloud Plus plan (more RAM)
