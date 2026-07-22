"""Tests for rosetta.wrappers.limma."""

import pandas as pd
import pytest
from rosetta._deps import is_installed
if not is_installed("limma"):
    pytest.skip("limma R package not installed", allow_module_level=True)
from rosetta._errors import RDataError, RFormulaError
from rosetta.wrappers.limma import Limma

# --- Fixtures ---

@pytest.fixture
def limma_model(sample_counts, sample_metadata):
    """Fixture to provide a fitted Limma model for testing."""
    model = Limma(sample_counts, sample_metadata, design="~ condition")
    model.run_ebayes(verbose=False)
    return model

# --- Input validation (no R/fitting required) ---

def test_negative_counts_raises(sample_counts, sample_metadata):
    bad_counts = sample_counts.copy()
    bad_counts.iloc[0, 0] = -1
    with pytest.raises(RDataError, match="negative"):
        Limma(bad_counts, sample_metadata, design="~ condition")

def test_mismatched_samples_raises(sample_counts, sample_metadata):
    bad_meta = sample_metadata.rename(index={"S1": "X1"})
    with pytest.raises(RDataError, match="columns must match"):
        Limma(sample_counts, bad_meta, design="~ condition")

def test_bad_formula_raises(sample_counts, sample_metadata):
    with pytest.raises(RFormulaError):
        Limma(sample_counts, sample_metadata, design="not a formula ~~~")

# --- Full pipeline tests ---

def test_limma_pipeline(limma_model):
    """Test standard model fitting and results extraction."""
    result = limma_model.get_results(coef=1)
    
    assert isinstance(result, pd.DataFrame)
    assert "logFC" in result.columns
    assert "adj.P.Val" in result.columns
    assert len(result) > 0

def test_apply_contrasts(sample_counts, sample_metadata):
    """Test contrast matrix application."""
    model = Limma(sample_counts, sample_metadata, design="~ condition")
    model.apply_contrasts("conditiontreated") 
    model.run_ebayes()
    
    result = model.get_results()
    assert isinstance(result, pd.DataFrame)
    assert "logFC" in result.columns
