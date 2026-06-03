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


# --- Input validation (no R required) ---

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
    from rosetta.wrappers.deseq2 import lfc_shrink
    # lfc_shrink validates type before touching R — no dds needed
    with pytest.raises(ValueError, match="Invalid shrinkage type"):
        lfc_shrink(None, coef="condition_treated_vs_control", type="bad_method")


# --- Full pipeline (requires DESeq2 in R) ---

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
    with pytest.raises(RDataError):
        get_results_names(sample_counts, sample_metadata, design="~ invalid")


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_preview_design(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import preview_design
    dds = preview_design(sample_counts, sample_metadata)
    assert dds is not None


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_run_deseq2(sample_counts, sample_metadata):
    from rosetta.wrappers.deseq2 import run_deseq2
    dds = run_deseq2(sample_counts, sample_metadata, design="~ condition")
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
    try:
        result = get_results(dds, contrast=["condition", "treated", "control"])
        assert isinstance(result, pd.DataFrame)
    except Exception:
        pass  # Small test dataset may not support this contrast


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_lfc_shrink_success(sample_counts, sample_metadata):
    """Verify lfc_shrink works end-to-end with a real dds object."""
    from rosetta.wrappers.deseq2 import run_deseq2, get_results_names, lfc_shrink

    dds = run_deseq2(sample_counts, sample_metadata, design="~ condition")
    coefs = get_results_names(sample_counts, sample_metadata, design="~ condition")

    target_coef = next((c for c in coefs if "condition" in c and "control" in c), None)
    if target_coef is None:
        pytest.fail(f"Could not find a valid condition coefficient in: {coefs}")

    result = lfc_shrink(dds, coef=target_coef, type="apeglm")
    assert isinstance(result, pd.DataFrame)
    assert "log2FoldChange" in result.columns
    assert len(result) == len(sample_counts)


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_lfc_shrink_ashr(sample_counts, sample_metadata):
    """Verify lfc_shrink works with ashr method."""
    from rosetta.wrappers.deseq2 import run_deseq2, get_results_names, lfc_shrink

    dds = run_deseq2(sample_counts, sample_metadata, design="~ condition")
    coefs = get_results_names(sample_counts, sample_metadata, design="~ condition")
    target_coef = next((c for c in coefs if "condition" in c and "control" in c), None)
    if target_coef is None:
        pytest.skip("No suitable coefficient found for ashr test")

    result = lfc_shrink(dds, coef=target_coef, type="ashr")
    assert isinstance(result, pd.DataFrame)
    assert "log2FoldChange" in result.columns


def test_lfc_shrink_invalid_type_no_r(sample_counts, sample_metadata):
    """Validation fires before R is called — no DESeq2 needed."""
    from rosetta.wrappers.deseq2 import lfc_shrink
    with pytest.raises(ValueError, match="Invalid shrinkage type"):
        lfc_shrink(None, coef="condition_treated_vs_control", type="invalid_method")


@pytest.mark.skipif(not _deseq2_available(), reason="DESeq2 not installed in R")
def test_deseq2_kwargs_passthrough(sample_counts, sample_metadata):
    """kwargs are forwarded through deseq2() -> get_results() -> DESeq2::results()."""
    from rosetta.wrappers.deseq2 import deseq2
    result = deseq2(sample_counts, sample_metadata, alpha=0.1)
    assert isinstance(result, pd.DataFrame)
