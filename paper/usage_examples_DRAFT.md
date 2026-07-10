# Usage Examples

Rosetta exposes the same underlying R/Bioconductor computation at three levels
of abstraction. A user chooses a tier based on how much control they need over
intermediate objects. The examples below all perform DESeq2 differential
expression on the same count matrix (`counts_df`, genes × samples) and sample
metadata (`metadata_df`), so the ergonomic tradeoff is directly visible.

**Tier 1 — Quick API.** A single call with sensible defaults, intended for
exploratory analysis and notebooks. The result is a `RosettaDataFrame` whose
`.report()` method prints a human-readable summary.

```python
import rosetta as rb

results = rb.quick_deseq2(counts_df, metadata_df, design="~ condition")
results.report()
```

```text
DESeq2 Results Summary
──────────────────────────────
Total genes tested:      200
Significant (padj<0.05): 30 (15.0%)
  ↑ Upregulated:         30
  ↓ Downregulated:       0
LFC range:               [1.57, 3.07]
```

**Tier 2 — Class-based.** A stateful object that keeps the underlying R object
alive across calls, for multi-step single-cell or microbiome pipelines. Here a
`Seurat` object runs the standard clustering pipeline and returns results as
pandas objects.

```python
from rosetta import Seurat

sc = Seurat(counts_df)
sc.run_standard_pipeline(n_variable_features=2000, n_pcs=10, resolution=0.5)
markers = sc.find_markers(ident_1="0", ident_2="1")
results = sc.get_results()
```

**Tier 3 — Functional.** Explicit step-by-step functions that map closely to the
native R signatures, giving full control over each intermediate object — useful
when a user needs a custom contrast or shrinkage estimator.

```python
from rosetta import run_deseq2, get_results, lfc_shrink

dds = run_deseq2(counts_df, metadata_df, design="~ condition")
res = get_results(dds, contrast=["condition", "treated", "control"], alpha=0.05)
shrunk = lfc_shrink(dds, coef="condition_treated_vs_control", type="apeglm")
```

Tiers 1 and 3 wrap the identical DESeq2 code path; Tier 1 simply fixes the
common defaults (`run_deseq2` followed by `get_results`) into one call. Tier 2
is distinct in that it retains R session state, which is required for the
iterative single-cell workflow.
