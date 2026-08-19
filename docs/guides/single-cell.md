# Single-Cell Analysis (Seurat) Guide

Complete guide to single-cell RNA-seq analysis using rosetta's Seurat wrapper.

---

## Overview

Rosetta wraps R's Seurat package via the class-based `Seurat` interface for standard
scRNA-seq workflows:

- Normalization (log-normalize or SCTransform)
- Feature selection (highly variable features)
- Scaling & PCA
- Graph-based clustering
- UMAP dimensionality reduction
- Differential expression & marker gene identification

All operations are chainable and method outputs return standard pandas DataFrames or
updated wrapper instances.

---

## Class-based workflow

### 1. Standard pipeline (All-in-one)

For a standard quick run, use `run_standard_pipeline()`:

```python
import rosetta as rb
seu = rb.Seurat(sc_counts).run_standard_pipeline(
n_variable_features=2000,
n_pcs=15,
resolution=0.8,
)
results = seu.get_results()

# Access formatted outputs:
clusters = results["clusters"] # pd.Series: cell → seurat_clusters
umap = results["umap"] # pd.DataFrame: UMAP_1, UMAP_2
var_features = results["variable_features"] # list of gene names
```

The standard pipeline runs these steps under the hood:

1. `NormalizeData()` — Log-normalization (`scale.factor=10000`)
2. `FindVariableFeatures()` — Select top highly variable genes (`nfeatures`)
3. `ScaleData()` — Z-score normalization across selected features
4. `RunPCA()` — Principal component analysis (`npcs`)
5. `FindNeighbors()` — Construct SNN graph using specified PCs
6. `FindClusters()` — Graph-based clustering (`resolution`)
7. `RunUMAP()` — Non-linear dimensional reduction for visualization

---

### 2. Step-by-step modular pipeline

You can also run step-by-step pipeline methods explicitly:

```python
import rosetta as rb
seu = (
rb.Seurat(sc_counts)
.run_normalize(normalization_method="LogNormalize")
.run_find_variable_features(nfeatures=2000)
.run_scale_data()
.run_pca(npcs=20)
.run_find_neighbors(dims=range(1, 21), k_param=20)
.run_find_clusters(resolution=0.6)
.run_umap(dims=range(1, 21))
)
```

---

### 3. SCTransform normalization

For alternative variance-stabilizing transformation:

```python
seu = (
rb.Seurat(sc_counts)
.run_sctransform(variable_features_n=3000)
.run_pca(npcs=30)
.run_find_neighbors(dims=range(1, 31))
.run_find_clusters(resolution=0.8)

.run_umap(dims=range(1, 31))
)
```

---

## Finding marker genes

Identify differentially expressed genes using `find_markers()`:

```python

# Markers for cluster 0 vs all other cells
markers_0 = seu.find_markers(ident_1="0")

# Markers between specific clusters
markers_0v1 = seu.find_markers(ident_1="0", ident_2="1")

# Group by a metadata column before finding markers
markers_by_condition = seu.find_markers(
ident_1="treated",
ident_2="control",
group_by="condition",
)
```

Returns a pandas DataFrame containing Seurat's marker statistics (`p_val`,
`avg_log2FC`, `pct.1`, `pct.2`, `p_val_adj`).

---

## Input data requirements

The `counts` parameter expects a pandas DataFrame:

- **Rows:** genes / features
- **Columns:** cells / barcodes
- **Values:** non-negative raw count matrix (integers)

```python
import pandas as pd

# Load CSV matrix
sc_counts = pd.read_csv("filtered_counts.csv", index_col=0)

# Convert from AnnData (scanpy)
sc_counts = pd.DataFrame(
adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
index=adata.obs_names,
columns=adata.var_names,
).T # Transpose to genes × cells
```

---

## Pipeline parameters

### `run_standard_pipeline()`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_variable_features` | `int` | `2000` | Number of highly variable features (`nfeatures`). |
| `n_pcs` | `int` | `10` | Number of principal components used for neighbors and UMAP. |
| `resolution` | `float` | `0.5` | Clustering resolution parameter for `FindClusters`. |

---

## Complete workflow

```python
import rosetta as rb
import pandas as pd

# 1. Load data
sc_counts = pd.read_csv("pbmc_counts.csv", index_col=0)

# 2. Run analysis
seu = rb.Seurat(sc_counts).run_standard_pipeline(
n_variable_features=3000,
n_pcs=20,
resolution=0.8,
)

# 3. Extract formatted results
results = seu.get_results()
print(f"Total cells: {len(results['clusters'])}")
print(f"Clusters found: {results['clusters'].nunique()}")
print(f"Variable features: {len(results['variable_features'])}")

# 4. Find marker genes for each cluster
for cluster_id in sorted(results["clusters"].unique()):
markers = seu.find_markers(ident_1=str(cluster_id))
top_5 = markers.head(5).index.tolist()
print(f"Cluster {cluster_id} markers: {top_5}")

# 5. Export UMAP embeddings
results["umap"].to_csv("umap_coords.csv")
```

---

## Requirements

R packages needed:

```r
install.packages(c("Seurat", "SeuratObject"))
```