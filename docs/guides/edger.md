# edgeR Guide

Complete guide to quasi-likelihood differential expression analysis with edgeR through rosetta. Covers basic usage, TREAT, custom contrasts, and multifactor designs.

---

## Overview

Rosetta wraps edgeR's recommended **quasi-likelihood (QL) pipeline** via the class-based `EdgeR` interface:

1. Create DGEList -> normalize (TMM) -> estimate dispersion -> fit QL model -> run test

This is the pipeline recommended by Gordon Smyth (edgeR author) for most RNA-seq experiments. It provides more reliable error rate control than the older exact test or LRT approaches.

---

## Basic analysis

The simplest case: two conditions, one comparison using the `EdgeR` class:

```python
import pandas as pd
from rosetta import EdgeR

# Your count matrix: genes as rows, samples as columns
# counts = pd.read_csv("counts.csv", index_col=0)
# metadata = pd.read_csv("metadata.csv", index_col=0)

# 1. Initialize and fit the model
model = EdgeR(counts, metadata, design="~ condition")

# 2. Run test and extract results
res_obj = model.run_test(lfc=0)
results = model.get_results(res_obj)
results.report()

# Filter significant genes
sig = results[results["FDR"] < 0.05].sort_values("FDR")
print(f"Found {len(sig)} significant genes")
```

---

## With LFC threshold (TREAT)

Use the `lfc` parameter in `run_test()` to test whether fold changes exceed a minimum threshold. This uses edgeR's `glmTreat()` function — a proper statistical test, not post-hoc filtering.

```python
# Test for genes with |logFC| > 1 (2-fold change)
res_obj = model.run_test(lfc=1.0)
results = model.get_results(res_obj)
results.report()
```

When `lfc > 0`, rosetta calls `glmTreat()` instead of `glmQLFTest()`. The null hypothesis changes from H₀: logFC = 0 to H₀: |logFC| ≤ threshold.

---

## Custom contrasts

### Numeric contrast vector

Specify which coefficient(s) to test using a numeric vector corresponding to the design matrix columns via the `contrast` parameter:

```python
# Test the second coefficient (e.g., condition effect in "~ condition")
res_obj = model.run_test(contrast=[0, 1])
results = model.get_results(res_obj)
```

### Multi-level factor

```python
metadata = pd.DataFrame({
    "group": ["A", "A", "B", "B", "C", "C"],
}, index=counts.columns)

model = EdgeR(counts, metadata, design="~ group")

# B vs A (second coefficient in a 3-level factor)
res_ba = model.get_results(model.run_test(contrast=[0, 1, 0]))

# C vs A
res_ca = model.get_results(model.run_test(contrast=[0, 0, 1]))

# B vs C (difference between non-reference levels)
res_bc = model.get_results(model.run_test(contrast=[0, 1, -1]))
```

---

## Multifactor designs

### Batch correction

```python
metadata = pd.DataFrame({
    "batch": ["A", "A", "B", "B", "A", "B"],
    "condition": ["ctrl", "ctrl", "ctrl", "treat", "treat", "treat"],
}, index=counts.columns)

model = EdgeR(counts, metadata, design="~ batch + condition")

# Test condition while adjusting for batch (assuming condition is the last coefficient)
res_obj = model.run_test(contrast=[0, 0, 1])
results = model.get_results(res_obj)
```

### Interaction model

```python
metadata = pd.DataFrame({
    "genotype": ["WT", "WT", "KO", "KO", "WT", "KO"],
    "treatment": ["ctrl", "drug", "ctrl", "drug", "drug", "ctrl"],
}, index=counts.columns)

model = EdgeR(counts, metadata, design="~ genotype * treatment")
res_obj = model.run_test()
results = model.get_results(res_obj)
```
---

## Combining TREAT with contrasts

You can use both `lfc` and `contrast` together inside `run_test()`:

```python
model = EdgeR(counts, metadata, design="~ batch + condition")

# Test if the condition effect has |logFC| > 0.5, adjusting for batch
res_obj = model.run_test(
    contrast=[0, 0, 1],  # condition coefficient
    lfc=0.5,
)
results = model.get_results(res_obj)
```

---

## Complete workflow

```python
import pandas as pd
from rosetta import EdgeR

# Load data
counts = pd.read_csv("counts.csv", index_col=0)
metadata = pd.read_csv("metadata.csv", index_col=0)

# Run edgeR QL pipeline using EdgeR class
model = EdgeR(counts, metadata, design="~ batch + condition")
res_obj = model.run_test()
results = model.get_results(res_obj)
results.report()

# Filter and export
sig = results[results["FDR"] < 0.05].sort_values("FDR")
print(f"\n{len(sig)} significant genes at FDR < 0.05")
print(sig.head(10))

sig.to_csv("edger_significant.csv")
```

---

## Error handling

```python
from rosetta import EdgeR
from rosetta._errors import RDataError, RFormulaError

try:
    model = EdgeR(counts, metadata, design="~ condition")
except RDataError as e:
    print(f"Data problem: {e}")
    # "Count matrix contains negative values"
    # "Count matrix columns must match metadata row names"
except RFormulaError as e:
    print(f"Formula problem: {e}")
```

---

## edgeR vs DESeq2 — when to use which

| Consideration | DESeq2 | edgeR QL |
|---------------|--------|----------|
| Small sample sizes (n < 5/group) | ✅ More conservative | Acceptable |
| Large sample sizes (n > 10/group) | Good | ✅ Faster, less conservative |
| LFC shrinkage | ✅ Built-in (apeglm) | Manual (via limma-voom) |
| TREAT thresholds | Via `lfcThreshold` | ✅ Via `glmTreat` |
| Speed | Slower | ✅ Faster |
| Batch effects | Both handle equally well | Both handle equally well |

In practice, both give similar results for well-powered experiments. Use whichever your field prefers for reviewers.