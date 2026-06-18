# Feature: `report()` method on results

**Good first issue** — ideal for new contributors.

## Summary

Add a `.report()` method that produces a human-readable summary of differential expression or enrichment results. Currently, results are returned as plain `pd.DataFrame` objects. Users must manually inspect columns to understand significance.

## Desired behavior

```python
import rosetta as rb

results = rb.get_results(dds, alpha=0.05)
results.report()
```

Output:
```
DESeq2 Results Summary
──────────────────────
Total genes tested:     12,000
Significant (padj<0.05):   843 (7.0%)
  ↑ Upregulated (LFC>0):   412
  ↓ Downregulated (LFC<0): 431
LFC range: [-4.2, 5.8]
```

## Implementation approach

1. Create `rosetta/results.py` with a `RosettaDataFrame` subclass of `pd.DataFrame`
2. Add a `.report()` method that detects result type (DESeq2/edgeR/limma/enrichment) by column names
3. Update wrappers to return `RosettaDataFrame` instead of plain `pd.DataFrame`

Key columns to detect:
- DESeq2/edgeR/limma: `padj`, `log2FoldChange`, `pvalue`
- Enrichment: `p.adjust`, `GeneRatio`, `Description`

## Acceptance criteria

- [ ] `results.report()` prints a summary for DESeq2 results
- [ ] `results.report()` prints a summary for enrichment results
- [ ] Inherits all normal DataFrame behavior (slicing, to_csv, etc.)
- [ ] Unit test in `tests/test_report.py`
