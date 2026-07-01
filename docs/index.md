# 🪨 rosetta

**Python interface to R/Bioconductor — pandas in, pandas out, `.report()` when you're done.**

---

## What is rosetta?

Rosetta gives you the full power of R/Bioconductor's statistical methods — DESeq2, edgeR, limma, clusterProfiler, Seurat, phyloseq — without writing R code, wrestling with rpy2 type conversions, or leaving your Python workflow.

- **One function, one result.** Call `deseq2()`, get a pandas DataFrame back.
- **Three API tiers.** One-liners for notebooks, class-based pipelines for interactive work, step-by-step functions for full control.
- **`.report()` everything.** Every result has a human-readable summary built in.
- **See the R code.** `codegen` mode prints the equivalent R commands so you can verify or reproduce.

Rosetta calls R — it doesn't reimplement it. All statistics run in the original, peer-reviewed R packages.

---

## Install

```bash
pip install rosetta-bioc
```

Requires R 4.0+ with Bioconductor. See [Installation](setup/install.md) for full setup instructions.

---

## 30-second demo

```python
import rosetta as rb

# DESeq2 differential expression — one call, pandas out
results = rb.deseq2(counts_df, metadata_df, design="~ condition")
results.report()
```

Output:

```
DESeq2 Results Summary
──────────────────────────────
Total genes tested:      12,000
Significant (padj<0.05): 843 (7.0%)
  ↑ Upregulated:         428
  ↓ Downregulated:       415
LFC range:               [-4.71, 3.50]
```

---

## Complete working example

```python
import pandas as pd
import numpy as np
from rosetta import deseq2

# Simulate RNA-seq counts: 1000 genes, 6 samples (3 control, 3 treated)
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

results = deseq2(counts=counts, metadata=metadata, design="~ condition")
print(results.sort_values("padj").head(10))
```

---

## Three-tier API at a glance

| Tier | Style | Use case |
|------|-------|----------|
| 1 — Quick | `quick_*()` | One-liners for notebooks |
| 2 — Class-based | `Seurat()`, `Phyloseq()` | Stateful, chainable workflows |
| 3 — Functional | `run_deseq2()` + `get_results()` | Full control, multi-step pipelines |

See the [Quickstart](quickstart.md) for examples of each tier.

---

## Supported R packages

| R Package | rosetta wrapper | What it does |
|-----------|-----------------|--------------|
| DESeq2 | `deseq2()`, `run_deseq2()`, `get_results()` | Differential expression (negative binomial) |
| edgeR | `edger()` | Quasi-likelihood differential expression |
| limma | `limma_voom()` | Linear models + voom precision weights |
| clusterProfiler | `ORA`, `GSEA` | GO/KEGG/Reactome pathway enrichment |
| Seurat | `Seurat()` | Single-cell RNA-seq |
| phyloseq | `Phyloseq()` | Microbiome diversity analysis |
