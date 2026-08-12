# Tier 3: Functional API & R Escape Hatch

The functional API and escape hatch give you step-by-step control over each analysis pipeline and direct access to underlying R objects. You manage intermediate objects, call individual functions, and have full access to every parameter and script the underlying R packages expose.
---

## DESeq2

### `deseq2()`

Legacy/convenience function that runs the full DESeq2 pipeline in one call.

```python
from rosetta import deseq2

results = deseq2(counts, metadata, design="~ condition")
```

```python
deseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str = "~ condition",
    **kwargs
) -> RosettaDataFrame
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `counts` | `pd.DataFrame` | *required* | Raw count matrix (genes × samples). |
| `metadata` | `pd.DataFrame` | *required* | Sample metadata, index matches count columns. |
| `design` | `str` | `"~ condition"` | R formula string. |
| `**kwargs` | | | Passed to `DESeq2::results()`. |

**Returns:** `RosettaDataFrame` with `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`.

---

### `run_deseq2()`

Fit the DESeq2 model. Returns the fitted R object for downstream use with `get_results()` and `lfc_shrink()`.

```python
from rosetta.wrappers.deseq2 import run_deseq2

dds = run_deseq2(counts, metadata, design="~ batch + condition")
```

```python
run_deseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str,
) -> DESeqDataSet  # R object
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `counts` | `pd.DataFrame` | Raw count matrix (genes × samples). Non-negative integers. |
| `metadata` | `pd.DataFrame` | Sample metadata. Index must contain all count column names. |
| `design` | `str` | R formula (e.g., `"~ condition"`, `"~ batch + condition"`). |

**Returns:** Fitted DESeqDataSet R object.

**Raises:**

- `RDataError` — negative values, mismatched samples, or model fitting failure.
- `RFormulaError` — invalid R formula syntax.

---

### `get_results()`

Extract results from a fitted DESeqDataSet.

```python
from rosetta.wrappers.deseq2 import get_results

res = get_results(dds, contrast=["condition", "treated", "control"], alpha=0.05)
```

```python
get_results(
    dds,
    contrast: list | None = None,
    lfc_threshold: float = 0.0,
    alpha: float = 0.1,
) -> RosettaDataFrame
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `dds` | R object | *required* | Fitted DESeqDataSet from `run_deseq2()`. |
| `contrast` | `list` or `None` | `None` | Three-element list: `[factor, numerator, denominator]`. |
| `lfc_threshold` | `float` | `0.0` | Log2 fold change threshold for hypothesis testing (not post-hoc filtering). |
| `alpha` | `float` | `0.1` | FDR cutoff for independent filtering. |

**Returns:** `RosettaDataFrame` with `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`.

---

### `lfc_shrink()`

Apply log2 fold change shrinkage for visualization and ranking.

```python
from rosetta.wrappers.deseq2 import lfc_shrink

shrunk = lfc_shrink(dds, coef="condition_treated_vs_control", type="apeglm")
```

```python
lfc_shrink(
    dds,
    coef: str,
    type: str = "apeglm",
    **kwargs
) -> pd.DataFrame
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `dds` | R object | *required* | Fitted DESeqDataSet. |
| `coef` | `str` | *required* | Coefficient name (use `get_results_names()` to discover valid names). |
| `type` | `str` | `"apeglm"` | Shrinkage method: `"apeglm"`, `"ashr"`, or `"normal"`. |

**Returns:** `pd.DataFrame` with shrunken log2FoldChange values.

!!! note
    `apeglm` requires the R package `apeglm` installed. `ashr` requires `ashr`.

---

### `get_results_names()`

Discover available coefficient names for a given design.

```python
from rosetta.wrappers.deseq2 import get_results_names

names = get_results_names(counts, metadata, design="~ batch + condition")
# e.g., ['Intercept', 'batch_B_vs_A', 'condition_treated_vs_control']
```

---

## edgeR

### `edger()`

Run the edgeR quasi-likelihood (QL) pipeline.

```python
from rosetta import edger

results = edger(counts, metadata, design="~ condition")
```

```python
edger(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str = "~ condition",
    contrast: list | None = None,
    lfc: float = 0,
    **kwargs
) -> RosettaDataFrame
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `counts` | `pd.DataFrame` | *required* | Raw count matrix (genes × samples). |
| `metadata` | `pd.DataFrame` | *required* | Sample metadata. |
| `design` | `str` | `"~ condition"` | R formula. |
| `contrast` | `list` or `None` | `None` | Numeric contrast vector (e.g., `[0, 1]`) or contrast matrix. |
| `lfc` | `float` | `0` | Log-FC threshold for `glmTreat()`. If > 0, uses TREAT instead of QL F-test. |
| `**kwargs` | | | Passed to `edgeR::glmQLFit()`. |

**Returns:** `RosettaDataFrame` with `logFC`, `logCPM`, `F`, `PValue`, `FDR`.

**Raises:** `RDataError`, `RFormulaError`.

---

## limma

### `limma_voom()`

Run limma-voom differential expression analysis using edgeR v4's `voomLmFit`.

```python
from rosetta import limma_voom

results = limma_voom(counts, metadata, design="~ condition")
```

```python
limma_voom(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str = "~ condition",
    contrast: list | None = None,
    decide_tests: bool = False,
    **kwargs
) -> pd.DataFrame
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `counts` | `pd.DataFrame` | *required* | Raw count matrix (genes × samples). |
| `metadata` | `pd.DataFrame` | *required* | Sample metadata. |
| `design` | `str` | `"~ condition"` | R formula. |
| `contrast` | `list` or `None` | `None` | Contrast specification for `limma::contrasts.fit()`. |
| `decide_tests` | `bool` | `False` | Whether to run `decideTests()` on the fit. |
| `**kwargs` | | | Passed to `edgeR::voomLmFit()`. |

**Returns:** `pd.DataFrame` with `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`.

---

## Enrichment — ORA

Over-Representation Analysis using clusterProfiler.

### `ORA.enrich_go()`

```python
from rosetta import ORA

go = ORA.enrich_go(gene_list, organism="org.Hs.eg.db", ont="BP")
```

```python
ORA.enrich_go(
    gene_list: list[str],
    organism: str = "org.Hs.eg.db",
    ont: str = "BP",
    pvalue_cutoff: float = 0.05,
    min_gs_size: int = 10,
    max_gs_size: int = 500,
    **kwargs
) -> pd.DataFrame
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gene_list` | `list[str]` | *required* | List of gene identifiers (Entrez IDs). |
| `organism` | `str` | `"org.Hs.eg.db"` | OrgDb annotation package name. |
| `ont` | `str` | `"BP"` | GO ontology: `"BP"` (Biological Process), `"MF"` (Molecular Function), `"CC"` (Cellular Component), or `"ALL"`. |
| `pvalue_cutoff` | `float` | `0.05` | Adjusted p-value threshold. |
| `min_gs_size` | `int` | `10` | Minimum gene set size. |
| `max_gs_size` | `int` | `500` | Maximum gene set size. |

### `ORA.enrich_kegg()`

```python
kegg = ORA.enrich_kegg(gene_list, organism="hsa")
```

```python
ORA.enrich_kegg(
    gene_list: list[str],
    organism: str = "hsa",
    pvalue_cutoff: float = 0.05,
    min_gs_size: int = 10,
    max_gs_size: int = 500,
    **kwargs
) -> pd.DataFrame
```

### `ORA.enrich_pathway()`

Reactome pathway enrichment (requires `ReactomePA` R package).

```python
reactome = ORA.enrich_pathway(gene_list)
```

### `ORA.enrich_custom()`

Custom gene set enrichment with user-provided term-to-gene mapping.

```python
custom = ORA.enrich_custom(gene_list, term2gene=my_mapping_df)
```

```python
ORA.enrich_custom(
    gene_list: list[str],
    term2gene: pd.DataFrame,
    min_gs_size: int = 10,
    max_gs_size: int = 500,
    **kwargs
) -> pd.DataFrame
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `term2gene` | `pd.DataFrame` | Two-column DataFrame: first column = term ID, second column = gene ID. |

---

## Enrichment — GSEA

Gene Set Enrichment Analysis using ranked gene lists.

### `GSEA.prepare_gene_list()`

Convert differential expression results into a ranked gene list for GSEA.

```python
from rosetta import GSEA

ranked = GSEA.prepare_gene_list(deseq2_results, gene_col="index", fc_col="log2FoldChange")
```

```python
GSEA.prepare_gene_list(
    df: pd.DataFrame,
    gene_col: str,
    fc_col: str = "log2FoldChange",
) -> pd.Series
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `df` | `pd.DataFrame` | *required* | DE results (e.g., from `get_results()`). |
| `gene_col` | `str` | *required* | Column containing gene IDs. Use `"index"` to use the DataFrame index. |
| `fc_col` | `str` | `"log2FoldChange"` | Column to rank by (sorted descending). |

**Returns:** `pd.Series` — named numeric vector sorted by fold change (descending).

### `GSEA.gse_go()`

```python
gsea_results = GSEA.gse_go(ranked, organism="org.Hs.eg.db", ont="BP")
```

```python
GSEA.gse_go(
    gene_list: pd.Series,
    organism: str = "org.Hs.eg.db",
    ont: str = "BP",
    pvalue_cutoff: float = 0.05,
    eps: float = 1e-10,
    **kwargs
) -> pd.DataFrame
```

### `GSEA.gse_kegg()`

```python
gsea_kegg = GSEA.gse_kegg(ranked, organism="hsa")
```

```python
GSEA.gse_kegg(
    gene_list: pd.Series,
    organism: str = "hsa",
    pvalue_cutoff: float = 0.05,
    eps: float = 1e-10,
    **kwargs
) -> pd.DataFrame
```

---

## Convenience aliases

These top-level aliases are available for quick access:

```python
import rosetta as rb

rb.enrich_go(...)       # → ORA.enrich_go()
rb.enrich_kegg(...)     # → ORA.enrich_kegg()
rb.enrich_pathway(...)  # → ORA.enrich_pathway()
rb.enrich_custom(...)   # → ORA.enrich_custom()
```
