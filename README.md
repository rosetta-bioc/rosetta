# 🪨 rosetta

**Python interface to R/Bioconductor — pandas in, pandas out, `.report()` when you're done.**

[![PyPI](https://img.shields.io/pypi/v/rosetta-bioc)](https://pypi.org/project/rosetta-bioc/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://img.shields.io/badge/tests-170%2B%20passing-brightgreen)]()

```bash
pip install rosetta-bioc
```

## 30-second demo

```python
import rosetta as rb

# DESeq2 differential expression — one call, pandas out
results = rb.deseq2(counts_df, metadata_df, design="~ condition")
results.report()
```
```
DESeq2 Results Summary
──────────────────────────────
Total genes tested:      12,000
Significant (padj<0.05): 843 (7.0%)
  ↑ Upregulated:         428
  ↓ Downregulated:       415
LFC range:               [-4.71, 3.50]
```

That's it. No R code. No rpy2 boilerplate. No type conversion. Just results.

## What it wraps

| R Package | Python | What it does |
|-----------|--------|--------------|
| DESeq2 | `rb.deseq2()` | Differential expression (negative binomial) |
| edgeR | `rb.edger()` | Quasi-likelihood differential expression |
| limma | `rb.limma_voom()` | Linear models + TREAT significance |
| clusterProfiler | `rb.enrich_go()` | GO/KEGG/Reactome pathway enrichment |
| phyloseq | `rb.phyloseq()` | Microbiome diversity analysis |
| Seurat | `rb.seurat()` | Single-cell RNA-seq |

All functions return a `RosettaDataFrame` (pandas DataFrame subclass) with a `.report()` method.

## Modular DESeq2 API

For more control, use the step-by-step interface:

```python
from rosetta.wrappers.deseq2 import run_deseq2, get_results, lfc_shrink

dds = run_deseq2(counts_df, metadata_df, design="~ condition")
res = get_results(dds, contrast=["condition", "treated", "control"], alpha=0.05)
shrunk = lfc_shrink(dds, coef="condition_treated_vs_control", type="apeglm")

res.report()
shrunk.report()
```

## Enrichment analysis

```python
import rosetta as rb

# Over-representation analysis
go_results = rb.enrich_go(gene_list, org_db="org.Hs.eg.db", ont="BP")
go_results.report()

# KEGG pathways
kegg = rb.enrich_kegg(gene_list, organism="hsa")
kegg.report()
```

## Setup

**Python side:**
```bash
pip install rosetta-bioc
```

**R side** (one-time):
```bash
Rscript install.R
```

Or manually:
```r
BiocManager::install(c("DESeq2", "edgeR", "limma", "clusterProfiler"))
```

**Posit Cloud:** See [docs/posit-cloud.md](docs/posit-cloud.md) for zero-config setup.

## Requirements

- Python 3.9+
- R 4.0+ with Bioconductor
- rpy2 ≥ 3.5

## Philosophy

1. **Rosetta calls R — it doesn't reimplement it.** All statistics run in the original, validated R packages.
2. **Pandas in, pandas out.** No R objects leak into your Python workflow.
3. **Fail early, fail clearly.** Input validation happens in Python before crossing the R boundary.
4. **`.report()` everything.** Results should be immediately interpretable without manual inspection.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Good first issues are labeled — start with [Issue #1: `report()` enhancements](https://github.com/rosetta-bioc/rosetta/issues/1).

## Acknowledgments

Built on [rpy2](https://rpy2.github.io/) and the extraordinary R/Bioconductor ecosystem. All credit for the statistical methods goes to the original R package authors.

GSoC 2026 · MIT License · [Nodes Bio](https://nodes.bio)
