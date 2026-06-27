# Validation Tests

These tests verify that rosetta output matches direct R output on published
datasets. They run real R (not mocked) and require Bioconductor packages to
be installed.

**Scope:** GSoC Phase 3, Week 9 (July 21 – August 25, 2026)

## What Goes Here

- `validate_deseq2.py` — rosetta DESeq2 results vs direct R DESeq2 on airway/pasilla
- `validate_edger.py` — rosetta edgeR results vs direct R edgeR
- `validate_limma.py` — rosetta limma-voom results vs direct R limma

## Acceptance Criteria

For each wrapper, rosetta output must match direct R output within:
- log2FoldChange: tolerance ± 1e-6
- padj: tolerance ± 1e-6
- Gene rankings: Spearman ρ > 0.999

## Running

```bash
pytest tests/validation/ -v
```

These are intentionally separate from the main test suite — they're slow
(real R calls, real datasets) and require Bioconductor. The main `pytest`
run uses mocked R calls for speed.
