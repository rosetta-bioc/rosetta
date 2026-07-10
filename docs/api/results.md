# Results & Reporting

All rosetta functions return objects with a `.report()` method that prints a human-readable summary. There are two result types depending on the API tier.

---

## `RosettaDataFrame`

A pandas `DataFrame` subclass returned by Tier 1 Quick API (DESeq2, edgeR) and Tier 3 functional calls. It behaves exactly like a regular DataFrame — you can sort, filter, merge, save to CSV — but adds a `.report()` method.

### Usage

```python
from rosetta import deseq2

results = deseq2(counts, metadata, design="~ condition")

# It's a DataFrame — all pandas operations work
sig = results[results["padj"] < 0.05]
results.to_csv("deseq2_results.csv")

# Plus .report() for summaries
results.report()
```

### `.report(alpha=0.05)`

Prints and returns a formatted summary string. Automatically detects the result type by column names.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `alpha` | `float` | `0.05` | Significance threshold for counting hits. |

**Returns:** `str` — the report text.

### Report formats by result type

**DESeq2** (detected by `padj` + `log2FoldChange` columns):

```
DESeq2 Results Summary
──────────────────────────────
Total genes tested:      12,000
Significant (padj<0.05): 843 (7.0%)
  ↑ Upregulated:         428
  ↓ Downregulated:       415
LFC range:               [-4.71, 3.50]
```

**edgeR** (detected by `FDR` + `logFC` columns):

```
edgeR Results Summary
──────────────────────────────
Total genes tested:      12,000
Significant (FDR<0.05):  756 (6.3%)
  ↑ Upregulated:         390
  ↓ Downregulated:       366
logFC range:             [-3.82, 4.15]
```

**limma** (detected by `adj.P.Val` + `logFC` columns):

```
limma Results Summary
──────────────────────────────
Total genes tested:      12,000
Significant (adj.P<0.05): 912 (7.6%)
  ↑ Upregulated:         470
  ↓ Downregulated:       442
logFC range:             [-3.12, 3.89]
```

**Enrichment** (detected by `p.adjust` + `GeneRatio` columns):

```
Enrichment Results Summary
──────────────────────────────
Total terms tested:      1,245
Significant (p.adj<0.05): 87
Top enriched terms:
  • regulation of cell proliferation (p=1.23e-08)
  • immune response (p=4.56e-07)
  • signal transduction (p=8.90e-06)
```

**Unknown format:**

```
Rosetta Results: 500 rows × 6 columns
```

---

## `QuickResult`

Returned by `quick_seurat()` and `quick_phyloseq()`. Wraps a dictionary of results with dict-like access and a `.report()` method.

### Usage

```python
import rosetta as rb

results = rb.quick_seurat(sc_counts, resolution=0.8)

# Dict-like access
results["clusters"]          # pd.Series
results["umap"]              # pd.DataFrame
results["variable_features"] # list[str]

# Standard dict protocol
"clusters" in results  # True
results.keys()         # dict_keys(["clusters", "umap", "variable_features"])
len(results)           # 3

# Report
results.report()
```

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `.data` | `dict` | The underlying result dictionary. |
| `.method` | `str` | Quick API method name (`"seurat"`, `"phyloseq"`, etc.). |
| `.metadata` | `dict` | Additional metadata about the analysis parameters. |

### `.report()`

Dispatches to method-specific formatting:

**Seurat:**

```
Seurat Quick Analysis Summary
──────────────────────────────
Total cells:             5,000
Clusters found:          8
Cluster sizes:
  Cluster 0: 892 cells
  Cluster 1: 756 cells
  Cluster 2: 643 cells
  ...
UMAP dimensions:         2
Variable features:       2,000
```

**Phyloseq:**

```
Phyloseq Diversity Summary
──────────────────────────────
Samples:                 24
Diversity metrics:
  Shannon: mean=2.341, sd=0.456, range=[1.234, 3.012]
  Simpson: mean=0.891, sd=0.034, range=[0.812, 0.945]
```

### Dict-like methods

| Method | Description |
|--------|-------------|
| `result[key]` | Get value by key (raises `KeyError` if missing). |
| `result.get(key, default)` | Get value with optional default. |
| `key in result` | Check if key exists. |
| `result.keys()` | All available keys. |
| `result.values()` | All values. |
| `result.items()` | Key-value pairs. |
| `len(result)` | Number of entries. |

---

## Codegen integration

Both result types work with rosetta's `codegen` module. Enable it before running analysis to see the equivalent R code:

```python
import rosetta as rb
rb.codegen.enable()

results = rb.quick_deseq2(counts, metadata, design="~ condition")
# Prints: R> library(DESeq2)
#         R> dds <- DESeqDataSetFromMatrix(...)
#         R> dds <- DESeq(dds)
#         R> res <- results(dds, alpha=0.05)

# Get the R code as a string
r_code = rb.codegen.last()
```
