# DESeq2 in Python — From 40 Lines of rpy2 to 3

*Run Bioconductor's gold-standard differential expression analysis without leaving your pandas workflow.*

---

## The Problem

DESeq2 is the most widely used tool for RNA-seq differential expression analysis. It's also written in R. If your analysis pipeline is in Python — and increasingly, most computational biology pipelines are — you have two options:

1. **Context-switch to R.** Export your data, open RStudio, run the analysis, import the results back. This breaks your notebook flow and makes reproducibility harder.

2. **Use rpy2 directly.** This keeps everything in Python, but the boilerplate is brutal:

```python
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import Formula
from rpy2.robjects.conversion import localconverter

pandas2ri.activate()
numpy2ri.activate()

deseq2 = importr("DESeq2")
base = importr("base")
stats = importr("stats")

# Convert pandas DataFrames to R objects
with localconverter(ro.default_converter + pandas2ri.converter):
    r_counts = ro.conversion.py2rpy(counts_df)
    r_metadata = ro.conversion.py2rpy(metadata_df)

# DESeq2 needs an integer matrix, not a data.frame
r_counts = base.as_matrix(r_counts)
storage_mode = ro.r("storage.mode")
storage_mode(r_counts, "integer")

# Build the DESeqDataSet
dds = deseq2.DESeqDataSetFromMatrix(
    countData=r_counts,
    colData=r_metadata,
    design=Formula("~ condition")
)

# Run the pipeline
dds = deseq2.DESeq(dds)
res = deseq2.results(dds)

# Convert back to pandas
with localconverter(ro.default_converter + pandas2ri.converter):
    results_df = ro.conversion.rpy2py(base.as_data_frame(res))

# Clean up index, column names, etc.
results_df.index.name = "gene"
```

That's 35-40 lines of code to do something conceptually simple: "run DESeq2 on my counts matrix." And it's fragile. Subtle bugs lurk in the type conversions — floats where DESeq2 expects integers, missing factor levels, rpy2 API changes between versions.

Every computational biologist who works in Python has written some version of this. Most have written it more than once.

## The Solution

[Rosetta](https://github.com/rosetta-bioc/rosetta) reduces the above to:

```python
import rosetta as rb

results = rb.deseq2(counts_df, metadata_df, design="~ condition")
```

Three lines. `counts_df` and `metadata_df` are pandas DataFrames. `results` is a pandas DataFrame. The full DESeq2 pipeline — size factor estimation, dispersion fitting, Wald test — runs in R exactly as it would natively. Rosetta handles all the type conversion, validation, and plumbing.

Call `.report()` for an instant summary:

```python
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

No manual inspection. No `results_df[results_df['padj'] < 0.05].shape[0]`. Just call `.report()`.

## Complete Working Example

Copy, paste, and run — this needs only Python, R, and DESeq2 installed:

```python
import pandas as pd
import numpy as np
import rosetta as rb

# Simulate RNA-seq counts: 1000 genes, 6 samples
np.random.seed(42)
counts = pd.DataFrame(
    np.random.negative_binomial(5, 0.1, size=(1000, 6)),
    index=[f"gene_{i}" for i in range(1000)],
    columns=["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"],
)

metadata = pd.DataFrame(
    {"condition": ["control"] * 3 + ["treated"] * 3},
    index=counts.columns,
)

# Run DESeq2 — one line
results = rb.deseq2(counts, metadata, design="~ condition")

# Explore results as a normal pandas DataFrame
sig = results.query("padj < 0.05").sort_values("log2FoldChange")
print(f"Significant genes: {len(sig)}")
print(sig.head(10))

# Or get the summary
results.report()
```

## Under the Hood

Rosetta is **not** a reimplementation of DESeq2. It's a translation layer. Here's what happens when you call `rb.deseq2()`:

1. **Validation** — Python-side checks: counts are non-negative integers, metadata index matches count columns, design formula references valid columns.
2. **Conversion** — pandas DataFrames are converted to R objects via rpy2 with correct types (integer matrix for counts, data.frame with factors for metadata).
3. **Execution** — The standard DESeq2 pipeline runs in R: `DESeqDataSetFromMatrix()` → `DESeq()` → `results()`.
4. **Return** — R results are converted back to a pandas DataFrame (a `RosettaDataFrame` subclass that adds `.report()`).

Want to see exactly what R code ran? Enable codegen:

```python
rb.codegen.enable()
results = rb.deseq2(counts, metadata, design="~ batch + condition")
```

```
  R> library(DESeq2)
  R> dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~ batch + condition)
  R> dds <- DESeq(dds)
  R> res <- results(dds, alpha=0.1)
```

You can paste this into R and get identical results. No black box.

## Three-Tier API

Rosetta offers three levels of abstraction depending on how much control you need:

### Tier 1 — Quick (one-liners for notebooks)

```python
results = rb.quick_deseq2(counts_df, metadata_df, design="~ condition")
clusters = rb.quick_seurat(matrix)
diversity = rb.quick_phyloseq(otu_table, sample_data)
```

### Tier 2 — Class-based (stateful, chainable workflows)

```python
seu = rb.Seurat(matrix).normalize().find_clusters().umap()
ps = rb.Phyloseq(otu_table, sample_data).rarefy(1000).ordinate("PCoA", "bray")
```

### Tier 3 — Functional (full control, explicit steps)

```python
from rosetta.wrappers.deseq2 import run_deseq2, get_results, lfc_shrink

dds = run_deseq2(counts, meta, design="~ batch + condition")
res = get_results(dds, contrast=["condition", "treated", "control"], alpha=0.05)
shrunk = lfc_shrink(dds, coef="condition_treated_vs_control", type="apeglm")
```

The functional tier exposes everything DESeq2 offers: contrasts, LFC thresholds (proper hypothesis testing, not post-hoc filtering), multiple shrinkage estimators, interaction terms, and multi-factor designs.

## Getting Started

**Install:**

```bash
pip install rosetta-bioc
```

**R dependencies** (one-time):

```r
BiocManager::install(c("DESeq2", "edgeR", "limma", "clusterProfiler"))
```

Or use the included setup script:

```bash
Rscript install.R
```

**Requirements:** Python 3.9+, R 4.0+, rpy2 ≥ 3.5.

**What's wrapped today:**

| R Package | What it does |
|-----------|--------------|
| DESeq2 | Differential expression (negative binomial) |
| edgeR | Quasi-likelihood DE |
| limma | Linear models + voom |
| clusterProfiler | GO/KEGG/Reactome enrichment |
| Seurat | Single-cell RNA-seq |
| phyloseq | Microbiome diversity |

---

Rosetta is MIT-licensed, developed as a GSoC 2026 project, and actively maintained. Contributions welcome.

- **GitHub:** https://github.com/rosetta-bioc/rosetta
- **PyPI:** https://pypi.org/project/rosetta-bioc/
- **Docs:** https://github.com/rosetta-bioc/rosetta#readme
