# Tier 2: Class-based API

The class-based API provides stateful, chainable wrappers around R's Seurat and phyloseq packages. Build up analysis state with method chaining, then extract results as pandas DataFrames.

---

## `Seurat`

Wraps the R Seurat package for single-cell RNA-seq analysis.

### Constructor

```python
import rosetta as rb

seu = rb.Seurat(counts)
```

```python
rb.Seurat(counts: pd.DataFrame)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `counts` | `pd.DataFrame` | Gene × cell count matrix. Must be non-negative and non-empty. |

**Raises:** `RDataError` if counts are empty or contain negative values.

---

### `.run_standard_pipeline()`

Execute the standard Seurat workflow: Normalize → FindVariableFeatures → Scale → PCA → FindNeighbors → FindClusters → UMAP.

```python
seu = rb.Seurat(counts).run_standard_pipeline(
    n_variable_features=2000,
    n_pcs=10,
    resolution=0.5,
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_variable_features` | `int` | `2000` | Number of highly variable features to select. |
| `n_pcs` | `int` | `10` | Number of PCs for dimensionality reduction. |
| `resolution` | `float` | `0.5` | Leiden/Louvain clustering resolution. |

**Returns:** `self` (for chaining).

---

### `.run_sctransform()`

Normalize using SCTransform (variance-stabilizing transformation).

```python
seu = rb.Seurat(counts).run_sctransform()
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `**kwargs` | | | Additional arguments passed to `Seurat::SCTransform()`. |

**Returns:** `self` (for chaining).

---

### `.find_markers()`

Find differentially expressed genes between two groups of cells.

```python
markers = seu.find_markers(ident_1="0", ident_2="1")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ident_1` | `str` | *required* | Identity class for the first group. |
| `ident_2` | `str` or `None` | `None` | Identity class for the second group. If `None`, compares against all other cells. |
| `group_by` | `str` or `None` | `None` | Metadata column to use for cell grouping. |
| `**kwargs` | | | Additional arguments passed to `Seurat::FindMarkers()`. |

**Returns:** `pd.DataFrame` with marker gene statistics.

---

### `.get_results()`

Extract the final analysis results as a Python dictionary.

```python
results = seu.get_results()
```

**Returns:** `dict` with keys:

| Key | Type | Description |
|-----|------|-------------|
| `"clusters"` | `pd.Series` | Cluster assignment per cell |
| `"umap"` | `pd.DataFrame` | UMAP embedding coordinates |
| `"variable_features"` | `list[str]` | Selected highly variable genes |

---

### Full example

```python
import rosetta as rb

seu = (
    rb.Seurat(sc_counts)
    .run_standard_pipeline(n_pcs=20, resolution=0.8)
)

results = seu.get_results()
print(f"Found {results['clusters'].nunique()} clusters in {len(results['clusters'])} cells")

# Find markers between cluster 0 and cluster 1
markers = seu.find_markers("0", "1")
print(markers.head())
```

---

## `Phyloseq`

Wraps the R phyloseq package for microbiome analysis.

### Constructor

```python
import rosetta as rb

ps = rb.Phyloseq(otu_table, sample_data=sample_meta, tax_table=taxonomy)
```

```python
rb.Phyloseq(
    otu_table: pd.DataFrame,
    sample_data: pd.DataFrame | None = None,
    tax_table: pd.DataFrame | None = None,
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `otu_table` | `pd.DataFrame` | *required* | OTU/ASV count table (taxa × samples). Must be non-negative and non-empty. |
| `sample_data` | `pd.DataFrame` or `None` | `None` | Sample metadata with index matching OTU table columns. |
| `tax_table` | `pd.DataFrame` or `None` | `None` | Taxonomy table with index matching OTU table rows. |

**Raises:** `RDataError` if OTU table is empty or contains negative values.

---

### `.estimate_richness()`

Calculate alpha diversity metrics.

```python
diversity = ps.estimate_richness(measures=["Shannon", "Simpson", "Chao1"])
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `measures` | `list[str]` or `None` | `None` (all) | Diversity metrics to compute. Options: `"Observed"`, `"Chao1"`, `"ACE"`, `"Shannon"`, `"Simpson"`, `"InvSimpson"`, `"Fisher"`. |

**Returns:** `pd.DataFrame` with one row per sample, one column per metric.

---

### `.run_ordination()`

Perform ordination (dimensionality reduction) on the community data.

```python
coords = ps.run_ordination(method="PCoA", distance="bray")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `method` | `str` | `"PCoA"` | Ordination method: `"PCoA"`, `"NMDS"`, `"DCA"`, `"CCA"`, `"RDA"`. |
| `distance` | `str` | `"bray"` | Distance metric: `"bray"`, `"jaccard"`, `"unifrac"`, `"wunifrac"`. |
| `**kwargs` | | | Additional arguments passed to `phyloseq::ordinate()`. |

**Returns:** `pd.DataFrame` of ordination coordinates (samples × axes).

---

### Full example

```python
import rosetta as rb
import pandas as pd
import numpy as np

# Create OTU table: 50 taxa, 10 samples
otu = pd.DataFrame(
    np.random.poisson(10, size=(50, 10)),
    index=[f"OTU_{i}" for i in range(50)],
    columns=[f"sample_{i}" for i in range(10)],
)

sample_meta = pd.DataFrame(
    {"site": ["gut"] * 5 + ["skin"] * 5},
    index=otu.columns,
)

ps = rb.Phyloseq(otu, sample_data=sample_meta)

# Alpha diversity
diversity = ps.estimate_richness(measures=["Shannon", "Simpson"])
print(diversity)

# Beta diversity ordination
pcoa = ps.run_ordination(method="PCoA", distance="bray")
print(pcoa.head())
```
