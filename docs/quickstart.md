# Quickstart

This guide walks through rosetta's three API tiers with working examples.

## Prerequisites

- Python 3.9+
- R 4.0+ with Bioconductor
- `pip install rosetta-bioc`
- R packages installed (at minimum: DESeq2). See [R Packages](setup/r-packages.md).

## Sample data

All examples below use this simulated RNA-seq dataset:

```python
import pandas as pd
import numpy as np

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
```

---

## Tier 1 — Quick API

One function call, one result. Designed for notebooks and exploratory analysis.

### DESeq2

```python
import rosetta as rb

results = rb.quick_deseq2(counts, metadata, design="~ condition")
results.report()
```

### edgeR

```python
results = rb.quick_edger(counts, metadata, design="~ condition")
results.report()
```

### Seurat (single-cell)

```python
results = rb.quick_seurat(sc_counts, n_variable_features=2000, resolution=0.5)
results.report()

# Access results dict-style
results["clusters"]
results["umap"]
```

### phyloseq (microbiome)

```python
results = rb.quick_phyloseq(otu_table, sample_data=sample_meta, measures=["Shannon", "Simpson"])
results.report()

# Access diversity DataFrame
results["diversity"]
```

---

## Tier 2 — Class-based API

Build up state, chain methods, extract results when ready. Best for interactive workflows where you want to inspect intermediate steps.

### Seurat

```python
import rosetta as rb

seu = (
    rb.Seurat(sc_counts)
    .run_standard_pipeline(n_variable_features=2000, n_pcs=10, resolution=0.5)
)

# Extract structured results
results = seu.get_results()
results["clusters"]  # pd.Series of cluster assignments
results["umap"]      # pd.DataFrame of UMAP coordinates

# Or use SCTransform normalization
seu = rb.Seurat(sc_counts).run_sctransform()

# Find marker genes
markers = seu.find_markers(ident_1="0", ident_2="1")
```

### Phyloseq

```python
import rosetta as rb

ps = rb.Phyloseq(otu_table, sample_data=sample_meta, tax_table=taxonomy)

# Alpha diversity
diversity = ps.estimate_richness(measures=["Shannon", "Simpson", "Chao1"])

# Beta diversity ordination
coords = ps.run_ordination(method="PCoA", distance="bray")
```

---

## Tier 3 — Functional API

Explicit steps, full parameter access. You control the pipeline — rosetta handles type conversion.

### DESeq2 (step-by-step)

```python
from rosetta.wrappers.deseq2 import run_deseq2, get_results, lfc_shrink

# Fit model
dds = run_deseq2(counts, metadata, design="~ condition")

# Extract results with contrast
res = get_results(dds, contrast=["condition", "treated", "control"], alpha=0.05)

# LFC shrinkage
shrunk = lfc_shrink(dds, coef="condition_treated_vs_control", type="apeglm")

res.report()
shrunk.report()
```

### edgeR

```python
from rosetta import edger

results = edger(counts, metadata, design="~ condition", lfc=1.0)
results.report()
```

### limma-voom

```python
from rosetta import limma_voom

results = limma_voom(counts, metadata, design="~ batch + condition")
results.report()
```

### Enrichment (ORA)

```python
from rosetta import enrich_go, enrich_kegg

go = enrich_go(gene_list, org_db="org.Hs.eg.db", ont="BP")
kegg = enrich_kegg(gene_list, organism="hsa")
go.report()
```

### Enrichment (GSEA)

```python
from rosetta import GSEA

# Prepare ranked gene list from DESeq2 results
ranked = GSEA.prepare_gene_list(deseq2_results, gene_col="index", fc_col="log2FoldChange")

# Run GSEA
gsea_go = GSEA.gse_go(ranked, organism="org.Hs.eg.db", ont="BP")
gsea_kegg = GSEA.gse_kegg(ranked, organism="hsa")
```

---

## See the R code

Turn on `codegen` to see exactly what R commands rosetta executes:

```python
import rosetta as rb
rb.codegen.enable()

from rosetta.wrappers.deseq2 import run_deseq2, get_results
dds = run_deseq2(counts, metadata, design="~ batch + condition")
res = get_results(dds, lfc_threshold=1.0)
```

Output:

```
  R> library(DESeq2)
  R> dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~ batch + condition)
  R> dds <- DESeq(dds)
  R> res <- results(dds, alpha=0.1, lfcThreshold=1.0)
```

Retrieve it programmatically:

```python
code = rb.codegen.last()  # Returns as a string
```

---

## What's next?

- **[DESeq2 Guide](guides/deseq2.md)** — contrasts, shrinkage, multifactor designs
- **[edgeR Guide](guides/edger.md)** — QL pipeline with TREAT thresholds
- **[Enrichment Guide](guides/enrichment.md)** — GO, KEGG, Reactome, custom gene sets
- **[API Reference](api/tier1-quick.md)** — full parameter documentation
