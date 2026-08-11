"""Tests for rosetta.wrappers.deseq2."""

"""Tests for rosetta.wrappers.deseq2."""

import pandas as pd
import pytest
from rosetta._deps import is_installed
if not is_installed("DESeq2"):
    pytest.skip("DESeq2 R package not installed", allow_module_level=True)
from rosetta._errors import RDataError, RFormulaError
from rosetta.wrappers.deseq2 import DESeq2

# --- Fixtures ---

@pytest.fixture
def deseq_model(sample_counts, sample_metadata):
    """Fixture to provide a fitted DESeq2 model for testing."""
    return DESeq2(sample_counts, sample_metadata, design="~ condition")

# --- Input validation (no R/fitting required) ---

def test_negative_counts_raises(sample_counts, sample_metadata):
    bad_counts = sample_counts.copy()
    bad_counts.iloc[0, 0] = -1
    with pytest.raises(RDataError, match="negative"):
        DESeq2(bad_counts, sample_metadata, design="~ condition")

def test_mismatched_samples_raises(sample_counts, sample_metadata):
    bad_meta = sample_metadata.rename(index={"S1": "X1"})
    with pytest.raises(RDataError, match="columns must match"):
        DESeq2(sample_counts, bad_meta, design="~ condition")

def test_bad_formula_raises(sample_counts, sample_metadata):
    with pytest.raises(RFormulaError):
        DESeq2(sample_counts, sample_metadata, design="not a formula ~~~")

def test_invalid_shrink_method_raises(deseq_model):
    # Testing logic before model fitting
    with pytest.raises(ValueError, match="Invalid shrinkage type"):
        deseq_model.lfc_shrink(coef="condition_treated_vs_control", type="bad_method")

# --- Full pipeline tests ---

def test_deseq2_pipeline(deseq_model):
    """Test standard model fitting and results extraction."""
    deseq_model.run_deseq()
    result = deseq_model.get_results(alpha=0.1)
    
    assert isinstance(result, pd.DataFrame)
    assert "log2FoldChange" in result.columns
    assert "padj" in result.columns

def test_lfc_shrink_success(deseq_model):
    """Verify the lfc_shrink integration using the class object."""
    deseq_model.run_deseq(verbose=False)

    # Use r_obj (Tier 3) to discover coefficients
    coefs = list(deseq_model.deseq_pkg.resultsNames(deseq_model.r_obj))
    target_coef = next((c for c in coefs if "condition" in c), coefs[-1])

    # apeglm requires the optional Bioconductor package; fall back to "normal"
    # which ships with DESeq2 itself and is always available.
    shrink_type = "apeglm" if is_installed("apeglm") else "normal"
    result = deseq_model.lfc_shrink(coef=target_coef, type=shrink_type)

    assert isinstance(result, pd.DataFrame)
    assert "log2FoldChange" in result.columns