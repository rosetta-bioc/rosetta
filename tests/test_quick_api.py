"""Tests for Tier 1 Quick API in rosetta."""

import pytest
import pandas as pd
import numpy as np
from rosetta import quick_seurat, quick_phyloseq, quick_deseq2, quick_edger
from rosetta.results import RosettaDataFrame
from rosetta._deps import is_installed
from rosetta.wrappers.seurat import _seurat_available
from rosetta.wrappers.phyloseq import _phyloseq_available

@pytest.fixture
def sc_counts():
    """Generate a larger simulated gene expression matrix to satisfy Seurat requirements."""
    np.random.seed(42)
    # Increase cells from 20 to 50 to satisfy the default n_neighbors=30 requirement
    data = np.random.poisson(lam=2, size=(100, 50))
    return pd.DataFrame(data, index=[f"G{i}" for i in range(100)], columns=[f"C{i}" for i in range(50)])

@pytest.mark.skipif(not _seurat_available(), reason="Seurat not installed")
def test_quick_seurat_api(sc_counts):
    """Test the Quick Seurat API functionality."""
    result = quick_seurat(sc_counts, n_pcs=2)
    assert isinstance(result, dict)
    assert "clusters" in result

@pytest.mark.skipif(not _phyloseq_available(), reason="phyloseq not installed")
def test_quick_phyloseq_api(sc_counts):
    """Test the Quick Phyloseq API functionality."""
    # Verify the quick_phyloseq entry point calculates alpha diversity correctly
    result = quick_phyloseq(sc_counts, measures=["Shannon"])
    assert isinstance(result, pd.DataFrame)
    assert "Shannon" in result.columns


# --- Tests for quick_deseq2 ---

@pytest.mark.skipif(not is_installed("DESeq2"), reason="DESeq2 not installed")
def test_quick_deseq2_basic(sample_counts, sample_metadata):
    """Test quick_deseq2 returns a RosettaDataFrame with expected columns."""
    result = quick_deseq2(sample_counts, sample_metadata, design="~ condition")
    assert isinstance(result, RosettaDataFrame)
    for col in ["baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]:
        assert col in result.columns, f"Missing column: {col}"


@pytest.mark.skipif(not is_installed("DESeq2"), reason="DESeq2 not installed")
def test_quick_deseq2_alpha(sample_counts, sample_metadata):
    """Test that the alpha parameter is passed through."""
    result = quick_deseq2(sample_counts, sample_metadata, design="~ condition", alpha=0.1)
    assert isinstance(result, RosettaDataFrame)


@pytest.mark.skipif(not is_installed("DESeq2"), reason="DESeq2 not installed")
def test_quick_deseq2_report(sample_counts, sample_metadata):
    """Test that .report() works on quick_deseq2 output."""
    result = quick_deseq2(sample_counts, sample_metadata, design="~ condition")
    # .report() should not raise
    result.report()


# --- Tests for quick_edger ---

@pytest.mark.skipif(not is_installed("edgeR"), reason="edgeR not installed")
def test_quick_edger_basic(sample_counts, sample_metadata):
    """Test quick_edger returns a RosettaDataFrame with expected columns."""
    result = quick_edger(sample_counts, sample_metadata, design="~ condition")
    assert isinstance(result, RosettaDataFrame)
    for col in ["logFC", "logCPM", "F", "PValue", "FDR"]:
        assert col in result.columns, f"Missing column: {col}"


@pytest.mark.skipif(not is_installed("edgeR"), reason="edgeR not installed")
def test_quick_edger_report(sample_counts, sample_metadata):
    """Test that .report() works on quick_edger output."""
    result = quick_edger(sample_counts, sample_metadata, design="~ condition")
    # .report() should not raise
    result.report()

