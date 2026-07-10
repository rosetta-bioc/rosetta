# Single-Cell Analysis (Seurat) Guide

Complete guide to single-cell RNA-seq analysis using rosetta's Seurat wrapper.

---

## Overview

Rosetta wraps the R Seurat package for standard scRNA-seq workflows:

- Normalization (log-normalize or SCTransform)
- Feature selection
- PCA, UMAP
- Clustering
- Marker gene identification

All operations are chainable and return pandas objects.

---

## Quick analysis (Tier 1)

For a one-call standard pipeline:

```python
import rosetta as rb

results = rb.quick_seurat(sc_counts, n_variable_features=2000, resolution=0.5)
results.report()

# Access results
clusters = results["clusters"]    # pd.Series: cell → cluster
umap = results["umap"]            # pd.DataFrame: UMAP_1, UMAP_2
features = results["variable_features"]  # list of gene names
```

Output:

```
Seurat Quick Analysis Summary
──────────────────────────────
Total cells:             5,000
Clusters found:          8
Cluster sizes:
  Cluster 0: 892 cells
  Cluster 1: 756 cells
  ...
UMAP dimensions:         2
Variable features:       2,000
```

---

## Class-based workflow (Tier 2)

For more control, use the `Seurat` class:

### Standard pipeline

```python
import rosetta as rb

seu = (
    rb.Seurat(sc_counts)
    .run_standard_pipeline(
        n_variable_features=2000,
        n_pcs=15,
        resolution=0.8,
    )
)

results = seu.get_results()
```

The standard pipeline runs these steps in order:

1. `NormalizeData()` — log-normalization
2. `FindVariableFeatures()` — select highly variable genes
3. `ScaleData()` — z-score normalization
4. `RunPCA()` — principal component analysis
5. `FindNeighbors()` — build SNN graph
6. `FindClusters()` — Leiden/Louvain clustering
7. `RunUMAP()` — dimensionality reduction for visualization

### SCTransform normalization

For an alternative normalization that handles technical variation better:

```python
seu = rb.Seurat(sc_counts).run_sctransform()
```

SCTransform replaces the Normalize → FindVariableFeatures → Scale steps with a single variance-stabilizing transformation. After SCTransform, you would typically continue with PCA/UMAP/clustering through the standard pipeline.

---

## Finding marker genes

Identify differentially expressed genes between clusters or groups:

```python
# Markers for cluster 0 vs all other cells
markers_0 = seu.find_markers(ident_1="0")

# Markers between specific clusters
markers_0v1 = seu.find_markers(ident_1="0", ident_2="1")

# Use a metadata column for grouping
markers_by_condition = seu.find_markers(
    ident_1="treated",
    ident_2="control",
    group_by="condition",
)
```

Returns a DataFrame with Seurat's marker statistics (p_val, avg_log2FC, pct.1, pct.2, p_val_adj).

---

## Input data format

The `counts` parameter expects a pandas DataFrame:

- **Rows:** genes/features
- **Columns:** cells/barcodes
- **Values:** raw UMI counts (non-negative integers)

```python
import pandas as pd

# From a CSV
sc_counts = pd.read_csv("filtered_counts.csv", index_col=0)

# From AnnData (scanpy)
sc_counts = pd.DataFrame(
    adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
    index=adata.obs_names,
    columns=adata.var_names,
).T  # Transpose to genes × cells
```

---

## Parameters

### `run_standard_pipeline()`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_variable_features` | `int` | `2000` | Number of highly variable genes to select. Higher = more features, slower. |
| `n_pcs` | `int` | `10` | Number of PCs. Use elbow plot heuristic to choose. Typical: 10–50. |
| `resolution` | `float` | `0.5` | Clustering resolution. Higher = more clusters. Typical: 0.1–2.0. |

### Choosing resolution

| Resolution | Expected clusters (for ~5000 cells) |
|------------|--------------------------------------|
| 0.1–0.3 | 3–5 broad populations |
| 0.5–0.8 | 8–12 (typical) |
| 1.0–2.0 | 15–25+ fine-grained subtypes |

---

## Complete workflow

```python
import rosetta as rb
import pandas as pd

# Load 10X Genomics data (or similar)
sc_counts = pd.read_csv("pbmc_counts.csv", index_col=0)

# Run analysis
seu = (
    rb.Seurat(sc_counts)
    .run_standard_pipeline(n_variable_features=3000, n_pcs=20, resolution=0.8)
)

# Get results
results = seu.get_results()

print(f"Cells: {len(results['clusters'])}")
print(f"Clusters: {results['clusters'].nunique()}")
print(f"Variable features: {len(results['variable_features'])}")

# Find marker genes for each cluster vs rest
for cluster_id in sorted(results["clusters"].unique()):
    markers = seu.find_markers(ident_1=str(cluster_id))
    top_markers = markers.head(5).index.tolist()
    print(f"Cluster {cluster_id}: {top_markers}")

# Export UMAP coordinates for external plotting
results["umap"].to_csv("umap_coords.csv")
```

---

## Requirements

R packages needed:

```r
install.packages("Seurat")
# Or via BiocManager:
BiocManager::install("Seurat")
```

Seurat v5+ is recommended. SeuratObject is installed automatically as a dependency.
