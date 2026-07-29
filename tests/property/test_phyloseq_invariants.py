"""Property-based tests for Phyloseq wrapper invariants using Hypothesis."""

import pytest
import numpy as np
from hypothesis import given, settings
from hypothesis.strategies import composite
import hypothesis.strategies as st
from tests.property.conftest import valid_count_matrix, valid_metadata
from rosetta._deps import is_installed
from rosetta.wrappers.phyloseq import Phyloseq

if not is_installed("phyloseq"):
    pytest.skip("phyloseq package not available", allow_module_level=True)

@composite
def phyloseq_input_strategy(draw):
    """Generate valid OTU abundance matrix and corresponding metadata using shared strategies."""
    # Use valid_count_matrix as our OTU abundance table (Taxa x Samples)
    otu_table = draw(valid_count_matrix(min_genes=10, max_genes=30, min_samples=4, max_samples=8))
    
    # Generate matching metadata using valid_metadata
    metadata = draw(valid_metadata(list(otu_table.columns)))
    
    return otu_table, metadata

@settings(max_examples=15, deadline=None)
@given(pair=phyloseq_input_strategy())
def test_phyloseq_invariants(pair):
    otu_table, metadata = pair
    
    # 1. Initialize Phyloseq model
    model = Phyloseq(otu_table=otu_table, sample_data=metadata)
    
    # 2. Test Alpha Diversity (estimate_richness)
    richness_res = model.estimate_richness()
    
    # Invariant 1: Result sample count must match input metadata size
    assert len(richness_res) == len(metadata), "Result sample count must match input metadata size"
    
    # Invariant 2: Sample names (row names) must correspond to sample IDs
    assert set(richness_res.index).issubset(set(metadata.index)), "Result index must correspond to sample IDs"
    
    # Invariant 3: Diversity metrics must be finite numeric values
    for col in richness_res.columns:
        valid_vals = richness_res[col].dropna()
        if len(valid_vals) > 0:
            assert np.isfinite(valid_vals).all(), f"Diversity metric {col} must contain only finite values"

    # 3. Test Ordination Analysis (run_ordination)
    ord_res = model.run_ordination(method="PCoA", distance="bray")
    
    # Invariant 4: Ordination coordinate rows must match sample count
    assert len(ord_res) == len(metadata), "Ordination sample coordinate count must match input metadata size"
    
    # Invariant 5: Ordination coordinates must be finite
    assert np.isfinite(ord_res.select_dtypes(include=[np.number])).all().all(), "Ordination coordinates must be finite numbers"