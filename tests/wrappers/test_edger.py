"""Tests for rosetta.wrappers.edger."""

import pandas as pd
import pytest

from rosetta._errors import RDataError, RFormulaError


def _edger_available():
    try:
        from rosetta._deps import is_installed
        return is_installed("edgeR")
    except Exception:
        return False


def test_negative_counts_raises(sample_counts, sample_metadata):
    from rosetta.wrappers.edger import edger
    bad_counts = sample_counts.copy()
    bad_counts.iloc[0, 0] = -1
    with pytest.raises(RDataError, match="negative"):
        edger(bad_counts, sample_metadata)


def test_mismatched_samples_raises(sample_counts, sample_metadata):
    from rosetta.wrappers.edger import edger
    bad_meta = sample_metadata.rename(index={"S1": "X1"})
    with pytest.raises(RDataError, match="columns must match"):
        edger(sample_counts, bad_meta)


def test_bad_formula_raises(sample_counts, sample_metadata):
    from rosetta.wrappers.edger import edger
    with pytest.raises(RFormulaError):
        edger(sample_counts, sample_metadata, design="not a formula ~~~")


@pytest.mark.skipif(not _edger_available(), reason="edgeR not installed in R")
def test_edger_returns_dataframe(sample_counts, sample_metadata):
    from rosetta.wrappers.edger import edger
    result = edger(sample_counts, sample_metadata, design="~ condition")
    assert isinstance(result, pd.DataFrame)
    assert "logFC" in result.columns
    assert "FDR" in result.columns
    assert len(result) == len(sample_counts)


@pytest.mark.skipif(not _edger_available(), reason="edgeR not installed in R")
def test_edger_with_contrast(sample_counts, sample_metadata):
    from rosetta.wrappers.edger import edger
    # Define a simple contrast vector to compare conditions.
    # This assumes the design matrix columns can accept this contrast vector.
    contrast = [0, 1] 
    result = edger(sample_counts, sample_metadata, design="~ condition", contrast=contrast)
    assert isinstance(result, pd.DataFrame)
    assert "logFC" in result.columns


@pytest.mark.skipif(not _edger_available(), reason="edgeR not installed in R")
def test_edger_with_treat(sample_counts, sample_metadata):
    from rosetta.wrappers.edger import edger
    # Test if lfc > 0 triggers the glmTreat pathway correctly
    result = edger(sample_counts, sample_metadata, design="~ condition", lfc=1.0)
    assert isinstance(result, pd.DataFrame)
    assert "logFC" in result.columns
