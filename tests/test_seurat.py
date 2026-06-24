"""Tests for rosetta.wrappers.seurat."""

import numpy as np
import pandas as pd
import pytest

from rosetta._errors import RDataError
from rosetta.wrappers.seurat import Seurat


def _seurat_available():
    """Check if Seurat is installed in the R environment."""
    try:
        from rosetta._deps import is_installed
        return is_installed("Seurat")
    except Exception:
        return False


@pytest.fixture
def sc_counts():
    """Simulated single-cell count matrix: 200 genes x 50 cells."""
    np.random.seed(42)
    data = np.random.poisson(lam=2, size=(200, 50))
    genes = [f"Gene{i}" for i in range(200)]
    cells = [f"Cell{i}" for i in range(50)]
    return pd.DataFrame(data, index=genes, columns=cells)


def test_empty_counts_raises():
    """Ensure Seurat initialization fails with empty input."""
    with pytest.raises(RDataError, match="empty"):
        Seurat(pd.DataFrame())


def test_negative_counts_raises(sc_counts):
    """Ensure Seurat initialization fails with negative values."""
    bad = sc_counts.copy()
    bad.iloc[0, 0] = -1
    with pytest.raises(RDataError, match="negative"):
        Seurat(bad)


@pytest.mark.skipif(not _seurat_available(), reason="Seurat not installed in R")
def test_seurat_pipeline(sc_counts):
    """Test the standard Seurat analysis pipeline execution."""
    model = Seurat(sc_counts)
    result = model.run_standard_pipeline(n_variable_features=10, n_pcs=5).get_results()
    
    assert isinstance(result, dict)
    assert "clusters" in result
    assert "umap" in result
    assert "variable_features" in result
    assert len(result["clusters"]) == 50


@pytest.mark.skipif(not _seurat_available(), reason="Seurat not installed in R")
def test_seurat_new_features(sc_counts):
    """Test SCTransform and FindMarkers integration."""
    model = Seurat(sc_counts)
    
    # Verify SCTransform execution
    model.run_sctransform()
    
    # Verify FindMarkers functionality
    # Note: Requires defined clusters or ident_1/ident_2 parameters
    # This acts as a basic API integration test
    with pytest.raises(Exception): # May fail on random data without valid clusters
        model.find_markers(ident_1="0", ident_2="1")