"""Property-based tests for EdgeR wrapper invariants using Hypothesis."""

import pytest
from hypothesis import given, settings
from hypothesis.strategies import composite
from tests.property.conftest import valid_count_matrix, valid_metadata
from rosetta._deps import is_installed
from rosetta.wrappers.edger import EdgeR

if not is_installed("edgeR"):
    pytest.skip("edgeR package not available", allow_module_level=True)

@composite
def edger_input_strategy(draw):
    counts = draw(valid_count_matrix(min_genes=5, max_genes=20, min_samples=4, max_samples=10))
    metadata = draw(valid_metadata(list(counts.columns)))
    return counts, metadata

@settings(max_examples=15, deadline=None)
@given(pair=edger_input_strategy())
def test_edger_invariants(pair):
    counts, metadata = pair
    
    # 1. Initialize EdgeR model
    model = EdgeR(counts, metadata, design="~ condition")
    
    # 2. Run the statistical test and capture the returned test result object
    res_obj = model.run_test()
    
    # 3. Get results using the test result object
    results = model.get_results(res_obj)
    
    # 4. Invariant 1: Result row count must match input gene count
    assert len(results) == len(counts), "Result row count must match input gene count"
    
    # 5. Invariant 2: P-values must be between 0 and 1
    pvalue_col = "PValue" if "PValue" in results.columns else ("p.value" if "p.value" in results.columns else None)
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