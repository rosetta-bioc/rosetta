# Microbiome Analysis (phyloseq) Guide

Complete guide to microbiome diversity analysis using rosetta's phyloseq wrapper.

---

## Overview

Rosetta wraps R's phyloseq package for standard microbiome analyses:

- Alpha diversity (within-sample richness)
- Beta diversity (between-sample dissimilarity)
- Ordination (PCoA, NMDS)

Input is OTU/ASV count tables as pandas DataFrames. Output is pandas DataFrames.

---

## Quick analysis (Tier 1)

For a one-call diversity analysis:

```python
import rosetta as rb

results = rb.quick_phyloseq(
    otu_table,
    sample_data=sample_meta,
    measures=["Shannon", "Simpson", "Chao1"],
)
results.report()
```

Output:

```
Phyloseq Diversity Summary
──────────────────────────────
Samples:                 24
Diversity metrics:
  Shannon: mean=2.341, sd=0.456, range=[1.234, 3.012]
  Simpson: mean=0.891, sd=0.034, range=[0.812, 0.945]
  Chao1: mean=156.2, sd=23.4, range=[112.0, 198.0]
```

Access the diversity DataFrame:

```python
diversity_df = results["diversity"]
print(diversity_df.head())
```

---

## Class-based workflow (Tier 2)

For interactive exploration with multiple analyses:

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

# Sample metadata
sample_meta = pd.DataFrame({
    "site": ["gut"] * 10 + ["skin"] * 10,
    "subject": [f"subj_{i % 5}" for i in range(20)],
}, index=otu_table.columns)

# Taxonomy table (optional)
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

Alpha diversity measures within-sample richness and evenness.

```python
diversity = ps.estimate_richness(measures=["Shannon", "Simpson", "Chao1", "Observed"])
print(diversity)
```

### Available metrics

| Metric | Description | Range |
|--------|-------------|-------|
| `"Observed"` | Number of observed taxa | 0 to ∞ |
| `"Chao1"` | Estimated total richness (accounts for unseen taxa) | ≥ Observed |
| `"ACE"` | Abundance-based coverage estimator | ≥ Observed |
| `"Shannon"` | Richness + evenness (information entropy) | 0 to ln(S) |
| `"Simpson"` | Probability two individuals are different taxa | 0 to 1 |
| `"InvSimpson"` | Inverse Simpson (effective number of taxa) | 1 to S |
| `"Fisher"` | Fisher's alpha (log-series model) | 0 to ∞ |

### Comparing groups

```python
diversity = ps.estimate_richness(measures=["Shannon"])

# Merge with metadata for group comparisons
diversity["site"] = sample_meta["site"]
print(diversity.groupby("site")["Shannon"].describe())
```

---

## Beta diversity & ordination

Beta diversity measures between-sample dissimilarity. Ordination reduces this to 2–3 dimensions for visualization.

```python
# PCoA with Bray-Curtis distance
coords = ps.run_ordination(method="PCoA", distance="bray")
print(coords.head())
```

### Ordination methods

| Method | Description | Best for |
|--------|-------------|----------|
| `"PCoA"` | Principal Coordinates Analysis | Default choice. Linear ordination of distances. |
| `"NMDS"` | Non-metric MDS | Non-linear relationships. Stress < 0.2 is acceptable. |
| `"DCA"` | Detrended Correspondence Analysis | Gradient data (ecology). |
| `"CCA"` | Canonical Correspondence Analysis | Constrained ordination with environmental variables. |
| `"RDA"` | Redundancy Analysis | Linear constrained ordination. |

### Distance metrics

| Metric | Description |
|--------|-------------|
| `"bray"` | Bray-Curtis dissimilarity (abundance-based, default) |
| `"jaccard"` | Jaccard index (presence/absence) |
| `"unifrac"` | Unweighted UniFrac (requires phylogenetic tree) |
| `"wunifrac"` | Weighted UniFrac (abundance + phylogeny) |

---

## Input data format

### OTU/ASV table

- **Rows:** taxa (OTUs, ASVs, or species)
- **Columns:** samples
- **Values:** non-negative integer counts

```python
# From QIIME2 output
otu_table = pd.read_csv("feature-table.csv", index_col=0)

# From DADA2 output (R → CSV)
otu_table = pd.read_csv("seqtab_nochim.csv", index_col=0).T
```

### Sample metadata

- **Index:** sample IDs matching OTU table columns
- **Columns:** experimental variables (site, treatment, timepoint, etc.)

### Taxonomy table (optional)

- **Index:** taxa IDs matching OTU table rows
- **Columns:** taxonomic ranks (Kingdom, Phylum, Class, Order, Family, Genus, Species)

---

## Complete workflow

```python
import rosetta as rb
import pandas as pd

# Load data
otu = pd.read_csv("otu_table.csv", index_col=0)
meta = pd.read_csv("sample_metadata.csv", index_col=0)
tax = pd.read_csv("taxonomy.csv", index_col=0)

# Create phyloseq object
ps = rb.Phyloseq(otu, sample_data=meta, tax_table=tax)

# Alpha diversity
alpha = ps.estimate_richness(measures=["Shannon", "Simpson", "Chao1"])
print("Alpha diversity per sample:")
print(alpha)

# Group comparison
alpha["treatment"] = meta["treatment"]
print("\nShannon by treatment group:")
print(alpha.groupby("treatment")["Shannon"].agg(["mean", "std"]))

# Beta diversity ordination
pcoa = ps.run_ordination(method="PCoA", distance="bray")
print("\nPCoA coordinates:")
print(pcoa.head())

# Export for visualization
alpha.to_csv("alpha_diversity.csv")
pcoa.to_csv("pcoa_coordinates.csv")
```

---

## Requirements

R packages needed:

```r
BiocManager::install("phyloseq")
```

For UniFrac distances, you also need a phylogenetic tree in your phyloseq object (not currently supported through the rosetta wrapper — use R directly for tree-based analyses).
