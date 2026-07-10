# Enrichment Analysis Guide

Over-Representation Analysis (ORA) and Gene Set Enrichment Analysis (GSEA) using clusterProfiler through rosetta.

---

## Overview

Rosetta provides two enrichment approaches:

| Approach | Class | When to use |
|----------|-------|-------------|
| **ORA** | `rosetta.ORA` | You have a gene list (e.g., significant DEGs) |
| **GSEA** | `rosetta.GSEA` | You have ranked results (all genes with fold changes) |

Both wrap R's clusterProfiler package and return pandas DataFrames with `.report()` support.

---

## Over-Representation Analysis (ORA)

ORA tests whether a predefined gene set is over-represented in your list of significant genes compared to background expectation.

### Gene Ontology (GO)

```python
import rosetta as rb

# From your DE results, extract significant gene IDs
sig_genes = results[results["padj"] < 0.05].index.tolist()

# GO Biological Process enrichment
go_bp = rb.enrich_go(sig_genes, organism="org.Hs.eg.db", ont="BP")
go_bp.report()
```

Output:

```
Enrichment Results Summary
──────────────────────────────
Total terms tested:      1,245
Significant (p.adj<0.05): 87
Top enriched terms:
  • regulation of cell proliferation (p=1.23e-08)
  • immune response (p=4.56e-07)
  • signal transduction (p=8.90e-06)
```

### GO ontologies

| Value | Ontology | Description |
|-------|----------|-------------|
| `"BP"` | Biological Process | What the gene product does biologically |
| `"MF"` | Molecular Function | Biochemical activity of the gene product |
| `"CC"` | Cellular Component | Where the gene product is located |
| `"ALL"` | All three | Tests all ontologies together |

### KEGG pathways

```python
kegg = rb.enrich_kegg(sig_genes, organism="hsa")
kegg.report()
```

Common organism codes:

| Code | Organism |
|------|----------|
| `"hsa"` | Human |
| `"mmu"` | Mouse |
| `"rno"` | Rat |
| `"dme"` | Drosophila |
| `"sce"` | Yeast |

### Reactome pathways

```python
reactome = rb.ORA.enrich_pathway(sig_genes)
```

!!! note
    Requires the `ReactomePA` R package: `BiocManager::install("ReactomePA")`

### Custom gene sets

Use your own term-to-gene mapping:

```python
import pandas as pd

# Two-column DataFrame: term → gene
custom_sets = pd.DataFrame({
    "term": ["pathway_A", "pathway_A", "pathway_B", "pathway_B", "pathway_B"],
    "gene": ["gene1", "gene2", "gene3", "gene4", "gene5"],
})

results = rb.ORA.enrich_custom(sig_genes, term2gene=custom_sets)
```

### Adjusting parameters

```python
go = rb.enrich_go(
    sig_genes,
    organism="org.Hs.eg.db",
    ont="BP",
    pvalue_cutoff=0.01,    # Stricter significance
    min_gs_size=15,        # Skip very small gene sets
    max_gs_size=300,       # Skip very large (less specific) gene sets
)
```

---

## Gene Set Enrichment Analysis (GSEA)

GSEA uses all genes ranked by fold change — no arbitrary significance cutoff needed. It detects coordinated changes in gene sets even when individual genes don't reach significance.

### Step 1: Prepare ranked gene list

```python
from rosetta import GSEA
from rosetta.wrappers.deseq2 import run_deseq2, get_results

# Run DE analysis
dds = run_deseq2(counts, metadata, design="~ condition")
res = get_results(dds)

# Create ranked list: gene names → log2FoldChange, sorted descending
ranked = GSEA.prepare_gene_list(res, gene_col="index", fc_col="log2FoldChange")
print(ranked.head())
# gene_42     3.21
# gene_157    2.98
# gene_891    2.76
# ...
```

!!! tip
    `gene_col="index"` uses the DataFrame's row index as gene identifiers. If your gene IDs are in a column, pass the column name instead.

### Step 2: Run GSEA

```python
# GO GSEA
gsea_go = GSEA.gse_go(ranked, organism="org.Hs.eg.db", ont="BP")

# KEGG GSEA
gsea_kegg = GSEA.gse_kegg(ranked, organism="hsa")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gene_list` | `pd.Series` | *required* | Named numeric series (gene → fold change), sorted descending. |
| `organism` | `str` | `"org.Hs.eg.db"` / `"hsa"` | Organism annotation. |
| `ont` | `str` | `"BP"` | GO ontology (for `gse_go`). |
| `pvalue_cutoff` | `float` | `0.05` | Adjusted p-value threshold for reporting. |
| `eps` | `float` | `1e-10` | Boundary for calculating p-values (numerical stability). |

---

## ORA vs GSEA — when to use which

| Feature | ORA | GSEA |
|---------|-----|------|
| Input | Gene list (significant DEGs) | All genes with fold changes |
| Requires cutoff | Yes (you pick padj threshold) | No |
| Detects subtle shifts | No | Yes |
| Interpretation | "Are my DEGs enriched for X?" | "Is gene set X coordinately up/down?" |
| Recommended for | Small, focused gene lists | Exploratory, unbiased analysis |

**Best practice:** Run both. ORA confirms enrichment in your significant hits; GSEA catches pathway-level effects that no single gene drives.

---

## Complete workflow

```python
import rosetta as rb
from rosetta import GSEA
from rosetta.wrappers.deseq2 import run_deseq2, get_results

# 1. Differential expression
dds = run_deseq2(counts, metadata, design="~ condition")
de_results = get_results(dds, alpha=0.05)

# 2. ORA on significant genes
sig_genes = de_results[de_results["padj"] < 0.05].index.tolist()
go_ora = rb.enrich_go(sig_genes, organism="org.Hs.eg.db", ont="BP")
go_ora.report()

# 3. GSEA on all genes
ranked = GSEA.prepare_gene_list(de_results, gene_col="index")
gsea_results = GSEA.gse_go(ranked, organism="org.Hs.eg.db", ont="BP")

# 4. Export
go_ora.to_csv("go_ora_results.csv")
gsea_results.to_csv("gsea_go_results.csv")
```

---

## Required R packages

| Function | R package |
|----------|-----------|
| `enrich_go()`, `gse_go()` | `clusterProfiler`, `org.Hs.eg.db` (or appropriate OrgDb) |
| `enrich_kegg()`, `gse_kegg()` | `clusterProfiler` |
| `enrich_pathway()` | `ReactomePA` |
| `enrich_custom()` | `clusterProfiler` |

Install all:

```r
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ReactomePA"))
```
