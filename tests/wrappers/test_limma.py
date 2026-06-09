"""Tests for rosetta.wrappers.limma."""

import pandas as pd
import pytest

from rosetta._errors import RDataError, RFormulaError


def _limma_available():
    try:
        from rosetta._deps import is_installed
        return is_installed("limma") and is_installed("edgeR")
    except Exception:
        return False


def test_negative_counts_raises(sample_counts, sample_metadata):
    from rosetta.wrappers.limma import limma_voom
    bad_counts = sample_counts.copy()
    bad_counts.iloc[0, 0] = -1
    with pytest.raises(RDataError, match="negative"):
        limma_voom(bad_counts, sample_metadata)


def test_mismatched_samples_raises(sample_counts, sample_metadata):
    from rosetta.wrappers.limma import limma_voom
    bad_meta = sample_metadata.rename(index={"S1": "X1"})
    with pytest.raises(RDataError, match="columns must match"):
        limma_voom(sample_counts, bad_meta)


def test_bad_formula_raises(sample_counts, sample_metadata):
    from rosetta.wrappers.limma import limma_voom
    with pytest.raises(RFormulaError):
        limma_voom(sample_counts, sample_metadata, design="not a formula ~~~")


@pytest.mark.skipif(not _limma_available(), reason="limma/edgeR not installed in R")
def test_limma_voom_returns_dataframe(sample_counts, sample_metadata):
    from rosetta.wrappers.limma import limma_voom
    result = limma_voom(sample_counts, sample_metadata, design="~ condition")
    assert isinstance(result, pd.DataFrame)
    assert "logFC" in result.columns
    assert "adj.P.Val" in result.columns
    assert len(result) == len(sample_counts)
