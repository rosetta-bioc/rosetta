# rosetta-bioc Social Posts

---

## Reddit — r/bioinformatics

**Title:**
> I got tired of copy-pasting between Python and R, so I wrapped DESeq2/edgeR/limma in a pandas API — and added a codegen mode that shows you every R line it runs

**Body:**

Hey r/bioinformatics — built something I've wanted for a while and figured others might too.

**`rosetta-bioc`** — run DESeq2, edgeR, limma-voom, clusterProfiler, phyloseq, and Seurat from Python. Pandas in, pandas out.

```python
import rosetta as rb

results = rb.deseq2(counts_df, metadata_df, design="~ batch + condition")
results.report()
```
```
DESeq2 Results Summary
──────────────────────────────
Total genes tested:      12,000
Significant (padj<0.05): 843 (7.0%)
  ↑ Upregulated:         428
  ↓ Downregulated:       415
```

**The part I actually think is interesting:**

```python
rb.codegen.enable()
results = rb.deseq2(counts_df, metadata_df, design="~ batch + condition")
```
```
R> library(DESeq2)
R> dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~ batch + condition)
R> dds <- DESeq(dds)
R> res <- results(dds, alpha=0.05)
```

`rb.codegen.last()` returns the R code as a string so you can paste it directly into R and reproduce independently. It's not reimplementing anything — it's literally calling R via rpy2. The stats are 100% the original packages.

It also supports the things that actually come up in real analyses: multi-factor designs, LFC thresholding, apeglm/ashr shrinkage, contrasts, voom with `duplicateCorrelation`, etc.

```bash
pip install rosetta-bioc
Rscript install.R  # installs the R packages
```

MIT license. Would love feedback, especially from anyone who has strong opinions about what a Python↔Bioconductor bridge should or shouldn't do.

---

## Reddit — r/bioinformatics**tools**

**Title:**
> rosetta-bioc — Python wrapper for DESeq2, edgeR, limma, clusterProfiler, phyloseq, Seurat. Pandas in, pandas out. Codegen shows the R code it runs.

**Body:**

**rosetta-bioc** wraps R/Bioconductor packages so you can call them from Python without writing any R.

| R Package | Python call | What it does |
|---|---|---|
| DESeq2 | `rb.deseq2()` | Differential expression |
| edgeR | `rb.edger()` | Quasi-likelihood DE |
| limma | `rb.limma_voom()` | Linear models + TREAT |
| clusterProfiler | `rb.enrich_go()` | GO/KEGG/Reactome enrichment |
| phyloseq | `rb.phyloseq()` | Microbiome diversity |
| Seurat | `rb.seurat()` | Single-cell RNA-seq |

**Codegen mode** — see exactly what R is running:

```python
rb.codegen.enable()
results = rb.deseq2(counts_df, meta_df, design="~ batch + condition")
# prints every R call live; rb.codegen.last() returns it as a string
```

**`.report()`** — instant human-readable summary on any result object.

```bash
pip install rosetta-bioc
Rscript install.R
```

- Python 3.9+ · R 4.0+ · rpy2 ≥ 3.5
- MIT license · 170+ tests passing
- GitHub: https://github.com/rosetta-bioc/rosetta

---

## Twitter / X — Thread

**Tweet 1:**
> 🧬 Just shipped `rosetta-bioc`: run DESeq2, edgeR, limma, and clusterProfiler from Python with zero R boilerplate.
>
> Pandas in. Pandas out. One line.
>
> But the real killer feature? 🧵

**Tweet 2:**
> You can see every line of R code it runs.
>
> ```python
> rb.codegen.enable()
> results = rb.deseq2(counts_df, meta_df, design="~ batch + condition")
> ```
> ```
> R> library(DESeq2)
> R> dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~ batch + condition)
> R> dds <- DESeq(dds)
> R> res <- results(dds, alpha=0.05)
> ```
>
> `rb.codegen.last()` returns it as a string. Paste straight into R to reproduce independently.
>
> It's not a black box. It's R — with a Python face.

**Tweet 3:**
> `.report()` gives you an instant human-readable summary:
>
> ```
> DESeq2 Results Summary
> ──────────────────────
> Total genes tested:       12,000
> Significant (padj<0.05):  843 (7.0%)
>   ↑ Upregulated:          428
>   ↓ Downregulated:        415
> ```
>
> No manual inspection. No magic numbers. Just results.

**Tweet 4:**
    > ```bash
    > pip install rosetta-bioc
    > Rscript install.R  # one-time
    > ```
    >
    > MIT license. 170+ tests passing.
    > GitHub → https://github.com/rosetta-bioc/rosetta
    >
    > Built by @NodesBio + @GSoC 2026 🙏

---

## LinkedIn

We built a Python wrapper for R/Bioconductor that does something unusual: it shows its work.

**rosetta-bioc** lets you run DESeq2, edgeR, limma, and clusterProfiler directly from Python — pandas DataFrames in, pandas DataFrames out. No R code, no rpy2 boilerplate, no type juggling.

But the feature that makes bioinformaticians actually trust it:

**`rb.codegen.enable()`** — every R function call gets printed in real time. Turn it on and you see exactly what's running under the hood. Copy it into R and reproduce your results independently.

Because in science, a black box isn't a feature. It's a bug.

```bash
pip install rosetta-bioc
```

MIT · 170+ tests · Google Summer of Code 2026

→ https://github.com/rosetta-bioc/rosetta
