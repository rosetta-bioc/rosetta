# Reddit Post — r/bioinformatics

---

**Title:** I got tired of writing 40 lines of rpy2 boilerplate for DESeq2, so I made a Python wrapper — looking for feedback

---

**Body:**

Every time I needed to run DESeq2 from a Python notebook, I'd end up writing something like this:

```python
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import Formula
from rpy2.robjects.conversion import localconverter

pandas2ri.activate()
numpy2ri.activate()

deseq2 = importr("DESeq2")
base = importr("base")

with localconverter(ro.default_converter + pandas2ri.converter):
    r_counts = ro.conversion.py2rpy(counts_df)
    r_metadata = ro.conversion.py2rpy(metadata_df)

r_counts = base.as_matrix(r_counts)
storage_mode = ro.r("storage.mode")
storage_mode(r_counts, "integer")

dds = deseq2.DESeqDataSetFromMatrix(
    countData=r_counts,
    colData=r_metadata,
    design=Formula("~ condition")
)
dds = deseq2.DESeq(dds)
res = deseq2.results(dds)

with localconverter(ro.default_converter + pandas2ri.converter):
    results_df = ro.conversion.rpy2py(base.as_data_frame(res))
```

That's ~30-40 lines before I even get a DataFrame I can work with. And it breaks in subtle ways — float vs integer matrices, factor levels, rpy2 version changes.

I built a library called **rosetta** that reduces this to:

```python
import rosetta as rb

results = rb.deseq2(counts_df, metadata_df, design="~ condition")
```

That's it. Pandas in, pandas out. The full DESeq2 pipeline runs in R — rosetta doesn't reimplement any statistics, it just handles the conversion layer via rpy2 internally.

**What it wraps:**

- **DESeq2** — differential expression (quick API + full step-by-step with contrasts, shrinkage, lfcThreshold)
- **edgeR** — quasi-likelihood DE
- **limma** — voom + TREAT
- **clusterProfiler** — GO/KEGG/Reactome enrichment (ORA and GSEA)
- **Seurat** — single-cell RNA-seq (normalize, cluster, UMAP)
- **phyloseq** — microbiome diversity

Every result has a `.report()` method that prints a human-readable summary (n significant, up/down split, LFC range, etc.).

There's also a `codegen` mode that prints the equivalent R code so you can verify what's actually running, or paste it into R to reproduce independently.

**Status:** v0.2.2, MIT license, GSoC 2026 project. Works on Python 3.9+ / R 4.0+.

- PyPI: `pip install rosetta-bioc`
- GitHub: https://github.com/rosetta-bioc/rosetta

I'd appreciate feedback on the API design, missing features, or anything that feels wrong. If you work with any of these R packages from Python, I'd love to know what pain points you hit that this could address (or doesn't yet).

Not affiliated with any of the original R package authors — just a wrapper that tries to make the Python ↔ R bridge less painful.
