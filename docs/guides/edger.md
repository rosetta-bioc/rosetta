# edgeR Guide

Complete guide to quasi-likelihood differential expression analysis with edgeR through rosetta.

---

## Overview

Rosetta wraps edgeR's recommended **quasi-likelihood (QL) pipeline**:

1. Create DGEList → normalize (TMM) → fit QL model → test

This is the pipeline recommended by Gordon Smyth (edgeR author) for most RNA-seq experiments. It provides more reliable error rate control than the older exact test or LRT approaches.

!!! note
    Rosetta skips `estimateDisp()` in the QL pipeline. Per edgeR v4 documentation, this step is only needed for diagnostic plots — `glmQLFit()` estimates dispersions internally.

---

## Basic analysis

```python
from rosetta import edger

results = edger(counts, metadata, design="~ condition")
results.report()
```

Output:

```
edgeR Results Summary
──────────────────────────────
Total genes tested:      12,000
Significant (FDR<0.05):  756 (6.3%)
  ↑ Upregulated:         390
  ↓ Downregulated:       366
logFC range:             [-3.82, 4.15]
```

### What happens internally

1. `DGEList(counts)` — create the edgeR data object
2. `calcNormFactors()` — TMM normalization
3. `model.matrix(design, data=metadata)` — build design matrix
4. `glmQLFit(dge, design)` — fit quasi-likelihood model
5. `glmQLFTest(fit)` — quasi-likelihood F-test
6. `topTags(res, n=Inf)` — extract all results

---

## With LFC threshold (TREAT)

Use the `lfc` parameter to test whether fold changes exceed a minimum threshold. This uses edgeR's `glmTreat()` function — a proper statistical test, not post-hoc filtering.

```python
# Test for genes with |logFC| > 1 (2-fold change)
results = edger(counts, metadata, design="~ condition", lfc=1.0)
results.report()
```

When `lfc > 0`, rosetta calls `glmTreat()` instead of `glmQLFTest()`. The null hypothesis changes from H₀: logFC = 0 to H₀: |logFC| ≤ threshold.

---

## Custom contrasts

### Numeric contrast vector

Specify which coefficient(s) to test using a numeric vector corresponding to the design matrix columns:

```python
# Test the second coefficient (e.g., condition effect in "~ condition")
results = edger(counts, metadata, design="~ condition", contrast=[0, 1])
```

### Multi-level factor

```python
metadata = pd.DataFrame({
    "group": ["A", "A", "B", "B", "C", "C"],
}, index=counts.columns)

# B vs A (second coefficient in a 3-level factor)
results = edger(counts, metadata, design="~ group", contrast=[0, 1, 0])

# C vs A
results = edger(counts, metadata, design="~ group", contrast=[0, 0, 1])

# B vs C (difference between non-reference levels)
results = edger(counts, metadata, design="~ group", contrast=[0, 1, -1])
```

---

## Multifactor designs

### Batch correction

```python
metadata = pd.DataFrame({
    "batch": ["A", "A", "B", "B", "A", "B"],
    "condition": ["ctrl", "ctrl", "ctrl", "treat", "treat", "treat"],
}, index=counts.columns)

# Test condition while adjusting for batch
results = edger(counts, metadata, design="~ batch + condition")
```

### Interaction model

```python
metadata = pd.DataFrame({
    "genotype": ["WT", "WT", "KO", "KO", "WT", "KO"],
    "treatment": ["ctrl", "drug", "ctrl", "drug", "drug", "ctrl"],
}, index=counts.columns)

results = edger(counts, metadata, design="~ genotype * treatment")
```

---

## Combining TREAT with contrasts

You can use both `lfc` and `contrast` together:

```python
# Test if the condition effect has |logFC| > 0.5, adjusting for batch
results = edger(
    counts, metadata,
    design="~ batch + condition",
    contrast=[0, 0, 1],  # condition coefficient
    lfc=0.5,
)
```

---

## Complete workflow

```python
import pandas as pd
from rosetta import edger

# Load data
counts = pd.read_csv("counts.csv", index_col=0)
metadata = pd.read_csv("metadata.csv", index_col=0)

# Run edgeR QL pipeline
results = edger(counts, metadata, design="~ batch + condition")
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
from rosetta._errors import RDataError, RFormulaError

try:
    results = edger(counts, metadata, design="~ condition")
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
