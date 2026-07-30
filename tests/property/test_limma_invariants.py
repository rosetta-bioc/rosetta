"""Property-based tests for Limma-voom wrapper invariants using Hypothesis."""

import pytest
from hypothesis import given, settings
from hypothesis.strategies import composite
from tests.property.conftest import valid_count_matrix, valid_metadata
from rosetta._deps import is_installed
from rosetta.wrappers.limma import Limma

if not is_installed("limma"):
    pytest.skip("limma package not available", allow_module_level=True)

@composite
def limma_input_strategy(draw):
    counts = draw(valid_count_matrix(min_genes=5, max_genes=20, min_samples=4, max_samples=10))
    metadata = draw(valid_metadata(list(counts.columns)))
    return counts, metadata

@settings(max_examples=15, deadline=None)
@given(pair=limma_input_strategy())
def test_limma_invariants(pair):
    counts, metadata = pair
    
    # 1. Initialize Limma model (automatically runs voomLmFit under the hood)
    model = Limma(counts, metadata, design="~ condition")
    
    # 2. Run empirical Bayes moderation
    model.run_ebayes()
    
    # 3. Get results using topTable wrapper
    results = model.get_results(alpha=0.05)
    
    # 4. Invariant 1: Result row count must match input gene count
    assert len(results) == len(counts), "Result row count must match input gene count"
    
    # 5. Invariant 2: P-values must be between 0 and 1
    pvalue_col = "adj.P.Val" if "adj.P.Val" in results.columns else ("P.Value" if "P.Value" in results.columns else None)
    if pvalue_col:
        valid_pval = results[pvalue_col].dropna()
        assert ((valid_pval >= 0.0) & (valid_pval <= 1.0)).all(), f"{pvalue_col} values must be between 0 and 1"
        
    # 6. Invariant 3: logFC must be finite values
    lfc_col = "logFC" if "logFC" in results.columns else None
    if lfc_col:
        valid_lfc = results[lfc_col].dropna()
        assert valid_lfc.notnull().all(), f"{lfc_col} cannot contain NaN"
        import numpy as np
        assert np.isfinite(valid_lfc).all(), f"{lfc_col} must be finite"