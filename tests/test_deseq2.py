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
