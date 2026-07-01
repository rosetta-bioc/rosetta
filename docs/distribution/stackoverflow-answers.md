# Stack Overflow Answer Drafts — rosetta-bioc

---

## Question 1: How to run DESeq2 from Python using rpy2?

**URL:** https://stackoverflow.com/questions/XXXXX/how-to-run-deseq2-from-python-using-rpy2

### Answer

The existing rpy2 answers here are solid — rpy2 is the right foundation for calling R from Python. If you find yourself writing the same conversion boilerplate repeatedly (pandas → R matrix, R data.frame → pandas, managing `importr` and `ro.conversion`), there's a library that wraps all of that into a one-liner:

```python
import rosetta as rb

results = rb.deseq2(counts_df, metadata_df, design="~ condition")
print(results.head())  # it's a pandas DataFrame
```

`rosetta-bioc` uses rpy2 under the hood but handles all the type conversion automatically — you pass in pandas DataFrames and get pandas DataFrames back. The full DESeq2 pipeline (size factor estimation, dispersion, Wald test) runs in R exactly as it would natively.

If you need more control (contrasts, LFC shrinkage, lfcThreshold):

```python
from rosetta.wrappers.deseq2 import run_deseq2, get_results, lfc_shrink

dds = run_deseq2(counts_df, meta_df, design="~ batch + condition")
res = get_results(dds, contrast=["condition", "treated", "control"])
shrunk = lfc_shrink(dds, coef="condition_treated_vs_control", type="apeglm")
```

Install with `pip install rosetta-bioc` (requires R 4.0+ and DESeq2 installed on the R side).

GitHub: https://github.com/rosetta-bioc/rosetta

---

## Question 2: Converting DESeq2 results to pandas DataFrame in Python

**URL:** https://stackoverflow.com/questions/XXXXX/converting-deseq2-results-to-pandas-dataframe-python

### Answer

The manual conversion approach (extracting slots, rebuilding column by column) works but is fragile — it breaks when DESeq2 changes its results structure, and you lose metadata like the `alpha` and `lfcThreshold` used.

If your goal is "run DESeq2, get a pandas DataFrame," `rosetta-bioc` does exactly that:

```python
import rosetta as rb

results = rb.deseq2(counts_df, metadata_df, design="~ condition")
# results is already a pandas DataFrame with columns:
# baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
```

The returned DataFrame is a proper pandas object — you can `.sort_values()`, `.query()`, merge it, save to CSV, whatever you'd normally do. It also has a `.report()` method for a quick summary:

```python
results.report()
# DESeq2 Results Summary
# Total genes tested: 12,000
# Significant (padj<0.05): 843 (7.0%)
```

Under the hood it's still rpy2 calling DESeq2 in R — no reimplementation of the statistics. It just automates the pandas↔R conversion that everyone ends up writing by hand.

`pip install rosetta-bioc` — requires R with DESeq2 installed.

---

## Question 3: rpy2 DESeq2 — error with data type conversion

**URL:** https://stackoverflow.com/questions/XXXXX/rpy2-deseq2-error-data-type-conversion

### Answer

This type conversion error is a common pain point with rpy2 + DESeq2. The issue is usually that DESeq2 expects an integer matrix for counts, but the pandas → R conversion produces floats, or the metadata factors aren't being set correctly.

The fix in raw rpy2 is to explicitly cast:

```python
from rpy2.robjects import numpy2ri
# ensure integer matrix, set factor levels, etc.
```

If you'd rather skip the conversion debugging entirely, `rosetta-bioc` handles all the type coercion internally (integer enforcement, factor setting, column alignment):

```python
import rosetta as rb

results = rb.deseq2(counts_df, metadata_df, design="~ condition")
# Just works — returns a pandas DataFrame
```

It validates inputs on the Python side before passing to R, so you get clear Python errors ("counts must be non-negative integers", "metadata index must match counts columns") instead of cryptic R tracebacks.

The library uses rpy2 internally — it's not a replacement for rpy2, just a convenience layer on top for Bioconductor workflows. Supports multi-factor designs, contrasts, LFC shrinkage, and all the standard DESeq2 options.

`pip install rosetta-bioc` | GitHub: https://github.com/rosetta-bioc/rosetta
