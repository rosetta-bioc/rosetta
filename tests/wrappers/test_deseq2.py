"""Tests for rosetta.wrappers.deseq2."""

import pandas as pd
import pytest

from rosetta._errors import RDataError, RFormulaError


def _deseq2_available():
    try:
        from rosetta._deps import is_installed
        return is_installed("DESeq2")
    except Exception:
        return False


def test_negative_counts_raises(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import deseq2
    bad_counts = sample_counts.copy()
    bad_counts.iloc[0, 0] = -1
    with pytest.raises(RDataError, match="negative"):
        deseq2(bad_counts, sample_metadata)


def test_mismatched_samples_raises(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import deseq2
    bad_meta = sample_metadata.rename(index={"S1": "X1"})
    with pytest.raises(RDataError, match="columns must match"):
        deseq2(sample_counts, bad_meta)


def test_bad_formula_raises(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import deseq2
    with pytest.raises(RFormulaError):
        deseq2(sample_counts, sample_metadata, design="not a formula ~~~")


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_deseq2_returns_dataframe(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import deseq2
    result = deseq2(sample_counts, sample_metadata, design="~ condition")
    assert isinstance(result, pd.DataFrame)
    assert "log2FoldChange" in result.columns
    assert "padj" in result.columns
    assert len(result) == len(sample_counts)


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_get_results_names_success(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import get_results_names
    
    names = get_results_names(sample_counts, sample_metadata, design="~ condition")
    
    assert isinstance(names, list)
    assert "Intercept" in names
    assert any("condition" in name for name in names)


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_get_results_names_bad_formula(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import get_results_names
    # Since 'invalid' is not in metadata, this should raise RDataError.
    with pytest.raises(RDataError):
        get_results_names(sample_counts, sample_metadata, design="~ invalid")


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_preview_design(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import preview_design
    # Verify that dds can be initialized without crashing
    dds = preview_design(sample_counts, sample_metadata)
    assert dds is not None


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_run_deseq2(sample_counts, sample_metadata):
    """Verify that run_deseq2 performs model fitting successfully."""
    from rosetta.wrappers.deseq2 import run_deseq2
    # 添加 design 參數
    dds = run_deseq2(sample_counts, sample_metadata, design="~ condition")
    assert dds is not None


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_lfc_shrink_success(sample_counts, sample_metadata):
    """Verify the lfc_shrink pipeline integration."""
    from rosetta.wrappers.deseq2 import run_deseq2, get_results_names, lfc_shrink
    
    # 1. Fit the model
    dds = run_deseq2(sample_counts, sample_metadata, design="~ condition")
    
    # 2. Get names and verify we have coefficients
    coefs = get_results_names(sample_counts, sample_metadata, design="~ condition")
    
    # Debug: Print coefficients if the test fails
    # Let's find the specific coefficient for 'condition'
    # Usually it's 'condition_treated_vs_control'
    target_coef = next((c for c in coefs if "condition" in c and "control" in c), None)
    
    if target_coef is None:
        pytest.fail(f"Could not find a valid condition coefficient in: {coefs}")
    
    # 3. Perform lfcShrink
    result = lfc_shrink(dds, coef=target_coef, type="apeglm")
    
    assert isinstance(result, pd.DataFrame)
    assert "log2FoldChange" in result.columns
    assert len(result) == len(sample_counts)

    print(f"DEBUG: Using coef={target_coef}")


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_lfc_shrink_invalid_type(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import run_deseq2, lfc_shrink
    
    dds = run_deseq2(sample_counts, sample_metadata, design="~ condition")
    
    # Verify that invalid type raises ValueError before hitting R
    with pytest.raises(ValueError, match="Invalid shrinkage type"):
        lfc_shrink(dds, coef="condition_treated_vs_control", type="invalid_method")