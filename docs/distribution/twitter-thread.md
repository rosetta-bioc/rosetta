# Twitter/X Thread — rosetta-bioc launch

---

## 1/6

Running DESeq2 from Python shouldn't require 40 lines of rpy2 boilerplate.

But if you've ever tried to call Bioconductor from a Jupyter notebook, you know the pain: type conversions, factor levels, integer matrices, cryptic R tracebacks…

We fixed it. 🧵

---

## 2/6

This is what DESeq2 in Python looks like today:

```python
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
pandas2ri.activate()
deseq2 = importr("DESeq2")
base = importr("base")
# ... 30 more lines of conversion,
# type coercion, formula objects,
# localconverter context managers ...
```

Every. Single. Time.

---

## 3/6

Here's what it looks like with rosetta:

```python
import rosetta as rb

results = rb.deseq2(counts_df, metadata_df, design="~ condition")
```

That's it. Pandas in. Pandas out. Full DESeq2 pipeline runs in R — no reimplementation.

---

## 4/6

Every result has `.report()`:

```
DESeq2 Results Summary
──────────────────────────────
Total genes tested:      12,000
Significant (padj<0.05): 843 (7.0%)
  ↑ Upregulated:         428
  ↓ Downregulated:       415
LFC range:               [-4.71, 3.50]
```

No more `df[df['padj'] < 0.05].shape[0]` to see what you got.

---

## 5/6

It's not just DESeq2. Rosetta wraps:

• edgeR — quasi-likelihood DE
• limma — voom + TREAT
• clusterProfiler — GO/KEGG/GSEA enrichment
• Seurat — single-cell (normalize → cluster → UMAP)
• phyloseq — microbiome diversity

Same pattern: pandas in, pandas out, `.report()` when done.

---

## 6/6

```bash
pip install rosetta-bioc
```

GitHub: https://github.com/rosetta-bioc/rosetta

MIT licensed. Python 3.9+ / R 4.0+. GSoC 2026 project.

`codegen` mode shows you the exact R code being run — no black box. Verify everything.

Feedback welcome. What R/Bioc packages do you wish had a Python interface?
