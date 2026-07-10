"""Tests for rosetta.wrappers.edger."""

import pandas as pd
import pytest
from rosetta._errors import RDataError, RFormulaError
from rosetta.wrappers.edger import EdgeR

# --- Fixtures ---

@pytest.fixture
def edger_model(sample_counts, sample_metadata):
    """Fixture to provide an initialized EdgeR model."""
    return EdgeR(sample_counts, sample_metadata, design="~ condition")

# --- Input validation ---

def test_negative_counts_raises(sample_counts, sample_metadata):
    bad_counts = sample_counts.copy()
    bad_counts.iloc[0, 0] = -1
    with pytest.raises(RDataError, match="negative"):
        EdgeR(bad_counts, sample_metadata, design="~ condition")

def test_mismatched_samples_raises(sample_counts, sample_metadata):
    bad_meta = sample_metadata.rename(index={"S1": "X1"})
    with pytest.raises(RDataError, match="columns must match"):
        EdgeR(sample_counts, bad_meta, design="~ condition")

def test_bad_formula_raises(sample_counts, sample_metadata):
    with pytest.raises(RFormulaError):
        EdgeR(sample_counts, sample_metadata, design="not a formula ~~~")

# --- Full pipeline tests ---

def test_edger_pipeline(edger_model):
    """Test standard QL-test pipeline."""
    res_obj = edger_model.run_test(lfc=0)
    result = edger_model.get_results(res_obj)
    
    assert isinstance(result, pd.DataFrame)
    assert "logFC" in result.columns
    assert "FDR" in result.columns

def test_edger_with_contrast(edger_model):
    """Test contrast application."""
    res_obj = edger_model.run_test(contrast=[0, 1])
    result = edger_model.get_results(res_obj)
    
    assert isinstance(result, pd.DataFrame)
    assert "logFC" in result.columns

def test_edger_with_treat(edger_model):
    """Test glmTreat pathway."""
    res_obj = edger_model.run_test(lfc=1.0)
    result = edger_model.get_results(res_obj)
    
    assert isinstance(result, pd.DataFrame)
    assert "logFC" in result.columns