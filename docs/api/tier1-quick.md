# Tier 1: Quick API

The Quick API provides one-call functions for the most common analyses. Each returns either a `RosettaDataFrame` (for DESeq2/edgeR) or a `QuickResult` (for Seurat/phyloseq), both of which support `.report()`.

---

## `quick_deseq2()`

Fit a DESeq2 model and extract results in one call.

```python
import rosetta as rb

results = rb.quick_deseq2(counts, metadata, design="~ condition", alpha=0.05)
results.report()
```

### Signature

```python
rb.quick_deseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str = "~ condition",
    alpha: float = 0.05,
    **kwargs
) -> RosettaDataFrame
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `counts` | `pd.DataFrame` | *required* | Gene count matrix (genes × samples). Row index = gene IDs, columns = sample IDs. |
| `metadata` | `pd.DataFrame` | *required* | Sample metadata. Row index must match `counts` column names. |
| `design` | `str` | `"~ condition"` | R formula string for the experimental design. |
| `alpha` | `float` | `0.05` | FDR significance cutoff passed to `DESeq2::results()`. |
| `**kwargs` | | | Additional arguments forwarded to `get_results()` (e.g., `contrast`, `lfc_threshold`). |

### Returns

`RosettaDataFrame` with columns: `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`.

### Example with contrast

```python
results = rb.quick_deseq2(
    counts, metadata,
    design="~ condition",
    contrast=["condition", "treated", "control"],
    alpha=0.01,
)
```

---

## `quick_edger()`

Run the edgeR quasi-likelihood pipeline in one call.

```python
results = rb.quick_edger(counts, metadata, design="~ condition")
results.report()
```

### Signature

```python
rb.quick_edger(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str = "~ condition",
    **kwargs
) -> RosettaDataFrame
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `counts` | `pd.DataFrame` | *required* | Gene count matrix (genes × samples). |
| `metadata` | `pd.DataFrame` | *required* | Sample metadata with index matching count columns. |
| `design` | `str` | `"~ condition"` | R formula string for the experimental design. |
| `**kwargs` | | | Additional arguments forwarded to `edger()` (e.g., `contrast`, `lfc`). |

### Returns

`RosettaDataFrame` with columns: `logFC`, `logCPM`, `F`, `PValue`, `FDR`.

### Example with LFC threshold (TREAT)

```python
results = rb.quick_edger(counts, metadata, design="~ condition", lfc=1.0)
```

---

## `quick_seurat()`

Run the standard Seurat single-cell pipeline and return structured results.

```python
results = rb.quick_seurat(sc_counts, n_variable_features=2000, resolution=0.8)
results.report()
```

### Signature

```python
rb.quick_seurat(
    counts: pd.DataFrame,
    **kwargs
) -> QuickResult
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `counts` | `pd.DataFrame` | *required* | Cell × gene count matrix (or gene × cell — Seurat transposes internally). |
| `n_variable_features` | `int` | `2000` | Number of highly variable features to select. |
| `n_pcs` | `int` | `10` | Number of principal components for PCA/neighbors/UMAP. |
| `resolution` | `float` | `0.5` | Clustering resolution (higher = more clusters). |

All keyword arguments are passed to `Seurat.run_standard_pipeline()`.

### Returns

`QuickResult` with dict-like access:

| Key | Type | Description |
|-----|------|-------------|
| `"clusters"` | `pd.Series` | Cluster assignment per cell |
| `"umap"` | `pd.DataFrame` | UMAP coordinates (2 columns) |
| `"variable_features"` | `list[str]` | Selected highly variable gene names |

---

## `quick_phyloseq()`

Compute alpha diversity metrics from an OTU table.

```python
results = rb.quick_phyloseq(otu_table, sample_data=meta, measures=["Shannon", "Simpson"])
results.report()
```

### Signature

```python
rb.quick_phyloseq(
    otu_table: pd.DataFrame,
    sample_data: pd.DataFrame | None = None,
    measures: list[str] = ["Shannon"],
    **kwargs
) -> QuickResult
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `otu_table` | `pd.DataFrame` | *required* | OTU/ASV count table (taxa × samples). |
| `sample_data` | `pd.DataFrame` or `None` | `None` | Sample metadata. |
| `measures` | `list[str]` | `["Shannon"]` | Diversity metrics to compute (e.g., `"Shannon"`, `"Simpson"`, `"Chao1"`, `"Observed"`). |
| `**kwargs` | | | Additional arguments passed to `Phyloseq()` constructor (e.g., `tax_table`). |

### Returns

`QuickResult` with dict-like access:

| Key | Type | Description |
|-----|------|-------------|
| `"diversity"` | `pd.DataFrame` | Diversity metrics per sample |
| `"measures"` | `list[str]` | Which metrics were computed |
