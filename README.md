# 🪨 rosetta

**Python interface to R bioinformatics packages — pandas in, pandas out.**

`rosetta` handles the rpy2 boilerplate, type conversion, and R dependency management so you can call DESeq2, edgeR, limma, and other Bioconductor tools from Python without writing any R.

## What it does

- Converts pandas DataFrames to R data.frames/matrices and back automatically
- Manages R package detection and provides clear errors when packages are missing
- Uses namespaced R calls (`DESeq2::DESeq()`, not global environment lookups) for safety
- Validates inputs on the Python side before crossing the R boundary

## What it doesn't do

- Reimplement any statistics — all computation happens in the original R packages
- Eliminate rpy2 — it's the bridge layer, with all its known fragility
- Expose every parameter — v0.1 covers default pipelines; granular options (shrinkage, contrasts) are in progress

## Quick Start

```bash
pip install rosetta
```

```python
import rosetta as rb

# DESeq2 differential expression
results = rb.deseq2(
    counts=counts_df,        # pandas DataFrame (genes × samples)
    metadata=metadata_df,    # pandas DataFrame with conditions
    design="~ condition"
)

results.head()
#          baseMean  log2FoldChange  lfcSE    pvalue      padj
# GeneA     1523.4           2.31   0.41  3.2e-08   1.1e-06
# GeneB      892.1          -1.87   0.38  5.6e-07   8.4e-06
```

## Supported Packages

| R Package | Status | Python Function |
|-----------|--------|----------------|
| DESeq2 | ✅ Ready | `rb.deseq2()` |
| edgeR | ✅ Ready | `rb.edger()` |
| limma | ✅ Ready | `rb.limma_voom()` |
| Seurat | ✅ Ready | `rb.seurat()` |
| clusterProfiler | ✅ Ready | `rb.enrichment()` |
| phyloseq | ✅ Ready | `rb.phyloseq()` |

## Requirements

- Python 3.9+
- R 4.0+ (with Bioconductor)
- rpy2

## Glossary

New to bioinformatics or unsure about a term? See [GLOSSARY.md](GLOSSARY.md) — everything explained from first principles.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for details.

## License

MIT

## Acknowledgments

Built on [rpy2](https://rpy2.github.io/) and the R/Bioconductor ecosystem. All credit for the statistical methods goes to the original R package authors.
