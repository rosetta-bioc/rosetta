"""Property-based tests for DESeq2 wrapper invariants using Hypothesis."""

import pytest
from hypothesis import given, settings
from hypothesis.strategies import composite
from rosetta._deps import is_installed
from rosetta.wrappers.deseq2 import DESeq2
from tests.property.conftest import valid_count_matrix, valid_metadata

if not is_installed("DESeq2"):
    pytest.skip("DESeq2 package not available", allow_module_level=True)

@composite
def count_and_metadata_strategy(draw):
    counts = draw(valid_count_matrix(min_genes=5, max_genes=20, min_samples=4, max_samples=10))
    metadata = draw(valid_metadata(list(counts.columns)))
    return counts, metadata

@settings(max_examples=15, deadline=None)
@given(pair=count_and_metadata_strategy())
def test_deseq2_invariants(pair):
    counts, metadata = pair
    
    # 1. Initialize and fit the model
    model = DESeq2(counts, metadata, design="~ condition")
    model.run_deseq(verbose=False)
    results = model.get_results(alpha=0.05)
    
    # 2. Invariant 1: Result row count must match input gene count
    assert len(results) == len(counts), "Result row count must match input gene count"
    
    # 3. Invariant 2: padj values must be between 0 and 1 (ignoring NA/NaN)
    if "padj" in results.columns:
        valid_padj = results["padj"].dropna()
        assert ((valid_padj >= 0.0) & (valid_padj <= 1.0)).all(), "padj values must be between 0 and 1"
        
    # 4. Invariant 3: log2FoldChange must be finite values or NA (no unexpected NaN or Inf)
    if "log2FoldChange" in results.columns:
        valid_lfc = results["log2FoldChange"].dropna()
        assert valid_lfc.notnull().all(), "log2FoldChange cannot contain NaN"
        import numpy as np
        assert np.isfinite(valid_lfc).all(), "log2FoldChange must be finite"