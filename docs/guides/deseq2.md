# DESeq2 Guide

Complete guide to differential expression analysis with DESeq2 through rosetta. Covers basic usage, contrasts, LFC shrinkage, and multifactor designs.

---

## Basic analysis

The simplest case: two conditions, one comparison.

```python
import pandas as pd
import numpy as np
from rosetta.wrappers.deseq2 import run_deseq2, get_results

# Your count matrix: genes as rows, samples as columns
# counts = pd.read_csv("counts.csv", index_col=0)

# Sample metadata
metadata = pd.DataFrame(
    {"condition": ["control", "control", "control", "treated", "treated", "treated"]},
    index=["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"],
)

# Fit the model
dds = run_deseq2(counts, metadata, design="~ condition")

# Extract results (default: last factor level vs first)
results = get_results(dds, alpha=0.05)
results.report()

# Filter significant genes
sig = results[results["padj"] < 0.05].sort_values("padj")
print(f"Found {len(sig)} significant genes")
```

---

## Contrasts

Specify exactly which comparison you want using the `contrast` parameter.

### Two-group contrast

```python
# Explicit: treated vs control
res = get_results(dds, contrast=["condition", "treated", "control"])
```

The contrast format is `[factor_name, numerator, denominator]`. Positive `log2FoldChange` means higher in the numerator.

### Multi-level factor

```python
metadata = pd.DataFrame({
    "genotype": ["WT", "WT", "KO", "KO", "rescue", "rescue"],
}, index=["s1", "s2", "s3", "s4", "s5", "s6"])

dds = run_deseq2(counts, metadata, design="~ genotype")

# KO vs WT
res_ko = get_results(dds, contrast=["genotype", "KO", "WT"])

# rescue vs KO
res_rescue = get_results(dds, contrast=["genotype", "rescue", "KO"])
```

### Discovering available coefficients

```python
from rosetta.wrappers.deseq2 import get_results_names

names = get_results_names(counts, metadata, design="~ genotype")
print(names)
# ['Intercept', 'genotype_KO_vs_WT', 'genotype_rescue_vs_WT']
```

---

## LFC threshold (hypothesis testing)

Use `lfc_threshold` for proper statistical testing of whether the fold change exceeds a threshold. This is **not** the same as post-hoc filtering — it adjusts the null hypothesis from H₀: LFC = 0 to H₀: |LFC| ≤ threshold.

```python
# Test for genes with |LFC| > 1 (2-fold change)
res = get_results(dds, lfc_threshold=1.0, alpha=0.05)
res.report()
```

!!! warning
    Post-hoc filtering (`results[abs(results["log2FoldChange"]) > 1]`) does not control FDR correctly. Use `lfc_threshold` for statistically rigorous thresholding.

---

## LFC shrinkage

Shrinkage reduces noise in fold change estimates, especially for low-count genes. Use it for **visualization and ranking** (volcano plots, gene lists), not for significance testing.

```python
from rosetta.wrappers.deseq2 import run_deseq2, get_results, lfc_shrink

dds = run_deseq2(counts, metadata, design="~ condition")

# Get unshrunk results (for p-values)
res = get_results(dds, alpha=0.05)

# Get shrunk LFCs (for ranking/plotting)
shrunk = lfc_shrink(dds, coef="condition_treated_vs_control", type="apeglm")
```

### Shrinkage methods

| Method | Description | R package required |
|--------|-------------|-------------------|
| `"apeglm"` | Recommended. Adaptive t prior. Fast, accurate. | `apeglm` |
| `"ashr"` | Adaptive shrinkage. Handles multimodal distributions. | `ashr` |
| `"normal"` | Original DESeq2 method. No extra packages needed. | (built-in) |

### Finding the correct `coef` name

The `coef` parameter must match one of the model's result names:

```python
from rosetta.wrappers.deseq2 import get_results_names

names = get_results_names(counts, metadata, design="~ condition")
print(names)
# ['Intercept', 'condition_treated_vs_control']

shrunk = lfc_shrink(dds, coef="condition_treated_vs_control", type="apeglm")
```

---

## Multifactor designs

### Batch correction

Control for known batch effects while testing for the condition of interest:

```python
metadata = pd.DataFrame({
    "batch": ["A", "A", "B", "B", "A", "B"],
    "condition": ["ctrl", "ctrl", "ctrl", "treat", "treat", "treat"],
}, index=counts.columns)

dds = run_deseq2(counts, metadata, design="~ batch + condition")

# The condition effect is adjusted for batch
res = get_results(dds, contrast=["condition", "treat", "ctrl"])
res.report()
```

!!! tip
    DESeq2 uses the **last variable** in the formula as the primary comparison. `~ batch + condition` means "test condition, controlling for batch."

### Interaction terms

Test whether the effect of treatment differs between genotypes:

```python
metadata = pd.DataFrame({
    "genotype": ["WT", "WT", "KO", "KO", "WT", "KO"],
    "treatment": ["ctrl", "drug", "ctrl", "drug", "drug", "ctrl"],
}, index=counts.columns)

dds = run_deseq2(counts, metadata, design="~ genotype + treatment + genotype:treatment")
```

---

## Codegen — see the R code

Enable codegen to see exactly what R commands rosetta executes:

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

Retrieve as a string for reproducibility:

```python
r_code = rb.codegen.last()
print(r_code)
```

---

## Complete workflow

```python
import pandas as pd
import numpy as np
from rosetta.wrappers.deseq2 import run_deseq2, get_results, lfc_shrink

# 1. Load data
counts = pd.read_csv("counts.csv", index_col=0)
metadata = pd.read_csv("metadata.csv", index_col=0)

# 2. Fit model
dds = run_deseq2(counts, metadata, design="~ batch + condition")

# 3. Extract results
res = get_results(dds, contrast=["condition", "treated", "control"], alpha=0.05)
res.report()

# 4. LFC shrinkage for visualization
shrunk = lfc_shrink(dds, coef="condition_treated_vs_control", type="apeglm")

# 5. Export
res.to_csv("deseq2_results.csv")
sig_genes = res[res["padj"] < 0.05].index.tolist()
```

---

## Error handling

Rosetta validates inputs before crossing the R boundary:

```python
from rosetta._errors import RDataError, RFormulaError

try:
    dds = run_deseq2(counts, metadata, design="~ condition")
except RDataError as e:
    print(f"Data problem: {e}")
    # e.g., "Count matrix contains negative values"
    # e.g., "Count matrix columns must match metadata row names"
except RFormulaError as e:
    print(f"Formula problem: {e}")
    # e.g., "Invalid design formula '~ nonexistent_column'"
```
