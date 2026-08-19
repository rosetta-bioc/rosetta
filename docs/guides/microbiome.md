# Microbiome Analysis (phyloseq) Guide

Complete guide to microbiome diversity analysis using rosetta's phyloseq wrapper.

---

## Overview

Rosetta wraps R's phyloseq package via the class-based `Phyloseq` interface for standard
microbiome analyses:
- **Alpha diversity** (within-sample richness and evenness)
- **Beta diversity &amp; Ordination** (PCoA, NMDS, etc.)
Input accepts OTU/ASV count tables as pandas DataFrames, returning standard pandas
DataFrames for easy integration with Python tools.

---

## Class-based workflow

### Creating a phyloseq object

```python
import rosetta as rb
import pandas as pd
import numpy as np

# OTU table: taxa (rows) × samples (columns)
otu_table = pd.DataFrame(
np.random.poisson(10, size=(100, 20)),
index=[f"OTU_{i}" for i in range(100)],
columns=[f"sample_{i}" for i in range(20)],
)

# Sample metadata (index matching otu_table columns)
sample_meta = pd.DataFrame({
"site": ["gut"] * 10 + ["skin"] * 10,
"subject": [f"subj_{i % 5}" for i in range(20)],
}, index=otu_table.columns)

# Taxonomy table (index matching otu_table rows, optional)
tax_table = pd.DataFrame({

"Kingdom": ["Bacteria"] * 100,
"Phylum": np.random.choice(["Firmicutes", "Bacteroidetes", "Proteobacteria"], 100),
"Genus": [f"Genus_{i}" for i in range(100)],
}, index=otu_table.index)

# Create phyloseq object
ps = rb.Phyloseq(otu_table, sample_data=sample_meta, tax_table=tax_table)
```

---

## Alpha diversity

Calculate alpha diversity metrics using `estimate_richness()`:

```python
diversity = ps.estimate_richness(measures=["Shannon", "Simpson", "Chao1", "Observed"])
print(diversity.head())
```

### Available metrics

| Metric | Description |
|--------|-------------|
| `"Observed"` | Number of observed taxa |
| `"Chao1"` | Estimated total richness (accounts for unseen taxa) |
| `"ACE"` | Abundance-based coverage estimator |
| `"Shannon"` | Richness + evenness (information entropy) |
| `"Simpson"` | Probability two individuals are different taxa |
| `"InvSimpson"` | Inverse Simpson |
| `"Fisher"` | Fisher&#39;s alpha |

### Group comparisons

```python
diversity = ps.estimate_richness(measures=["Shannon"])

# Merge with metadata for analysis
diversity["site"] = sample_meta["site"]
print(diversity.groupby("site")["Shannon"].describe())
```

---

## Beta diversity & ordination

Perform ordination analysis using `run_ordination()`, which returns sample coordinates as
a pandas DataFrame:

```python
# PCoA with Bray-Curtis distance
coords = ps.run_ordination(method="PCoA", distance="bray")
print(coords.head())
```

### Supported parameters (`run_ordination`)

Supported parameters filter down to phyloseq's `ordinate` arguments:

| Parameter | Type | Description |
|-----------|------|-------------|
| `method` | `str` | Ordination method (`"PCoA"`, `"NMDS"`, `"DCA"`, `"CCA"`, `"RDA"`) |
| `distance` | `str` | Distance metric (`"bray"`, `"jaccard"`, etc.) |
| `k` | `int` | Number of dimensions |
| `trymax` | `int` | Maximum iterations (e.g. for NMDS) |

---

## Input data requirements

- **OTU/ASV Table:** Must be a non-empty `pd.DataFrame` with non-negative integer
counts. Taxa must be on rows (`taxa_are_rows=True`) and samples on columns.
- **Sample Metadata:** `pd.DataFrame` where the index matches OTU table column names.
- **Taxonomy Table:** `pd.DataFrame` where the index matches OTU table row names.

---

## Complete workflow

```python
import rosetta as rb
import pandas as pd

# Load datasets
otu = pd.read_csv("otu_table.csv", index_col=0)
meta = pd.read_csv("sample_metadata.csv", index_col=0)

tax = pd.read_csv("taxonomy.csv", index_col=0)

# Initialize wrapper
ps = rb.Phyloseq(otu, sample_data=meta, tax_table=tax)

# 1. Alpha diversity
alpha = ps.estimate_richness(measures=["Shannon", "Simpson", "Chao1"])
alpha["treatment"] = meta["treatment"]
print(alpha.groupby("treatment")["Shannon"].agg(["mean", "std"]))

# 2. Ordination (PCoA)
pcoa = ps.run_ordination(method="PCoA", distance="bray")
print(pcoa.head())

# Export
alpha.to_csv("alpha_diversity.csv")
pcoa.to_csv("pcoa_coordinates.csv")
```

---

## Requirements

R package needed:

```r
BiocManager::install("phyloseq")
```