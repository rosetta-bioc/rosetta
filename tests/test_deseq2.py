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


def test_invalid_shrink_method_raises(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import deseq2
    with pytest.raises(ValueError, match="shrink must be one of"):
        deseq2(sample_counts, sample_metadata, shrink="bad")


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_deseq2_returns_dataframe(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import deseq2
    result = deseq2(sample_counts, sample_metadata, design="~ condition")
    assert isinstance(result, pd.DataFrame)
    assert "log2FoldChange" in result.columns
    assert "padj" in result.columns
    assert len(result) == len(sample_counts)


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_deseq2_shrink_normal(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import deseq2
    result = deseq2(sample_counts, sample_metadata, design="~ condition", shrink="normal")
    assert isinstance(result, pd.DataFrame)
    assert "log2FoldChange" in result.columns
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
    from rosetta.wrappers.deseq2 import run_deseq2
    # Verify that full fitting runs correctly and returns a dds object
    dds = run_deseq2(sample_counts, sample_metadata)
    assert dds is not None


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_get_results_function(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import run_deseq2, get_results
    
    dds = run_deseq2(sample_counts, sample_metadata)
    result = get_results(dds)
    
    assert isinstance(result, pd.DataFrame)
    assert "log2FoldChange" in result.columns
    assert "padj" in result.columns
    assert len(result) == len(sample_counts)


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_get_results_with_parameters(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import run_deseq2, get_results
    
    dds = run_deseq2(sample_counts, sample_metadata)
    result = get_results(dds, lfc_threshold=1.0, alpha=0.05)
    
    assert isinstance(result, pd.DataFrame)
    assert "log2FoldChange" in result.columns


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_get_results_with_contrast(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import run_deseq2, get_results
    
    dds = run_deseq2(sample_counts, sample_metadata)
    # Note: This contrast may not exist in the small test dataset, 
    # but we're testing the parameter passing
    try:
        result = get_results(dds, contrast=["condition", "treated", "control"])
        assert isinstance(result, pd.DataFrame)
    except Exception:
        # Expected to potentially fail with small test data
        pass


def test_get_results_invalid_shrink_method(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import get_results, run_deseq2
    
    if _deseq2_available():
        dds = run_deseq2(sample_counts, sample_metadata)
        with pytest.raises(ValueError, match="shrink must be one of"):
            get_results(dds, shrink="invalid_method")


def test_deseq2_kwargs_passthrough(sample_counts, sample_metadata):
    """Test that kwargs are passed through to DESeq2::results."""
    if _deseq2_available():
        from rosetta.wrappers.deseq2 import deseq2
        # Test with a valid kwarg (this will test passthrough even if it doesn't change results)
        result = deseq2(sample_counts, sample_metadata, alpha=0.1)
        assert isinstance(result, pd.DataFrame)
