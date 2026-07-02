# Biostars Post — Tutorial

---

**Title:** Tutorial: Running DESeq2, edgeR, and limma from Python without rpy2 boilerplate

**Tags:** deseq2, rna-seq, python, edger, limma, differential-expression, bioconductor

---

**Type:** Tutorial

---

One of the most common questions on Biostars is "how do I run DESeq2 from Python?" The answer usually involves 30-40 lines of rpy2 conversion code — managing integer matrices, factor levels, localconverters, and R formula objects.

I built [rosetta-bioc](https://github.com/rosetta-bioc/rosetta) to eliminate that boilerplate. It calls the real R packages (no reimplementation), handles all the type conversion automatically, and returns pandas DataFrames.

## Quick example — DESeq2

```python
import rosetta as rb

results = rb.deseq2(counts_df, metadata_df, design="~ condition")
results.report()
```

```
DESeq2 Results Summary
──────────────────────────────
Total genes tested:      12,000
Significant (padj<0.05): 843 (7.0%)
  ↑ Upregulated:         428
  ↓ Downregulated:       415
LFC range:               [-4.71, 3.50]
```

## Full control when you need it

```python
from rosetta.wrappers.deseq2 import run_deseq2, get_results, lfc_shrink

dds = run_deseq2(counts, meta, design="~ batch + condition")
res = get_results(dds, contrast=["condition", "treated", "control"], lfc_threshold=1.0)
shrunk = lfc_shrink(dds, coef="condition_treated_vs_control", type="apeglm")
```

## It's not just DESeq2

| R Package | Python API | What it does |
|-----------|-----------|--------------|
| DESeq2 | `rb.quick_deseq2()` | Differential expression |
| edgeR | `rb.quick_edger()` | Quasi-likelihood DE |
| limma | `rb.limma_voom()` | Linear models + TREAT |
| clusterProfiler | `rb.enrich_go()`, `rb.ORA`, `rb.GSEA` | Pathway enrichment |
| Seurat | `rb.quick_seurat()` | Single-cell RNA-seq |
| phyloseq | `rb.quick_phyloseq()` | Microbiome diversity |
| VariantAnnotation | `rb.VCF()` | VCF reading + annotation |

## Normalization transforms

```python
from rosetta.wrappers.normalize import vst, rlog, tmm_normalize

# DESeq2 variance-stabilizing transformation
transformed = vst(counts_df, metadata_df)

# edgeR TMM normalization → log-CPM
normalized = tmm_normalize(counts_df)
```

## Show me the R code

Don't trust a black box? Turn on codegen:

```python
rb.codegen.enable()
dds = rb.wrappers.deseq2.run_deseq2(counts, meta, design="~ batch + condition")
```
```
  R> library(DESeq2)
  R> dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~ batch + condition)
  R> dds <- DESeq(dds)
```

## Plotting

```python
from rosetta.plots import volcano, ma_plot, pca

volcano(results, alpha=0.05, lfc_cutoff=1.0)
pca(counts_df, metadata_df, color_by="condition")
```

## Install

```bash
pip install rosetta-bioc
```

Requires Python 3.9+ and R 4.0+ with the relevant Bioconductor packages installed on the R side.

- **GitHub:** https://github.com/rosetta-bioc/rosetta
- **PyPI:** https://pypi.org/project/rosetta-bioc/
- **Docs:** https://rosetta-bioc.github.io/rosetta/

MIT licensed, GSoC 2026 project. Feedback welcome — especially on which R/Bioconductor packages to wrap next.

---

# Biostars Answer Template (for existing questions about DESeq2/Python/rpy2)

---

**For questions like:**
- "How to run DESeq2 in Python"
- "rpy2 DESeq2 conversion error"
- "Call Bioconductor from Python"
- "DESeq2 alternative in Python"

---

**Answer:**

If you're looking for a simpler interface than raw rpy2, [rosetta-bioc](https://pypi.org/project/rosetta-bioc/) wraps the conversion layer:

```python
import rosetta as rb

results = rb.deseq2(counts_df, metadata_df, design="~ condition")
# returns pandas DataFrame with baseMean, log2FoldChange, pvalue, padj
```

It still runs the real R/DESeq2 under the hood (via rpy2) — no reimplementation. Handles the integer matrix conversion, formula objects, and R↔pandas bridging automatically.

For multi-factor designs, contrasts, and LFC shrinkage:

```python
from rosetta.wrappers.deseq2 import run_deseq2, get_results, lfc_shrink

dds = run_deseq2(counts, meta, design="~ batch + condition")
res = get_results(dds, contrast=["condition", "treated", "control"])
shrunk = lfc_shrink(dds, coef="condition_treated_vs_control", type="apeglm")
```

`pip install rosetta-bioc` — requires R 4.0+ with DESeq2 installed.

GitHub: https://github.com/rosetta-bioc/rosetta
