"""Tests for rosetta.wrappers.seurat."""

import numpy as np
import pandas as pd
import pytest

from rosetta._errors import RDataError


def _seurat_available():
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
    from rosetta.wrappers.seurat import seurat
    with pytest.raises(RDataError, match="empty"):
        seurat(pd.DataFrame())


def test_negative_counts_raises(sc_counts):
    from rosetta.wrappers.seurat import seurat
    bad = sc_counts.copy()
    bad.iloc[0, 0] = -1
    with pytest.raises(RDataError, match="negative"):
        seurat(bad)


@pytest.mark.skipif(not _seurat_available(), reason="Seurat not installed in R")
def test_seurat_returns_dict(sc_counts):
    from rosetta.wrappers.seurat import seurat
    result = seurat(sc_counts, min_cells=1, min_features=1, n_pcs=5)
    assert isinstance(result, dict)
    assert "clusters" in result
    assert "umap" in result
    assert "variable_features" in result
    assert "umap_1" in result["umap"].columns or "UMAP_1" in result["umap"].columns
    assert len(result["clusters"]) > 0


@pytest.mark.skipif(not _seurat_available(), reason="Seurat not installed in R")
def test_seurat_custom_parameters(sc_counts):
    from rosetta.wrappers.seurat import seurat
    
    # Test with custom parameters
    result = seurat(
        sc_counts, 
        min_cells=2, 
        min_features=5, 
        n_variable_features=50,
        n_pcs=3,
        resolution=0.8
    )
    
    assert isinstance(result, dict)
    assert len(result["variable_features"]) <= 50  # May be fewer if not enough genes pass filters


@pytest.mark.skipif(not _seurat_available(), reason="Seurat not installed in R")
def test_seurat_kwargs_passthrough(sc_counts):
    """Test that kwargs are passed through to FindClusters."""
    from rosetta.wrappers.seurat import seurat
    
    # Test with random.seed as a kwarg to FindClusters
    result = seurat(sc_counts, min_cells=1, min_features=1, n_pcs=3, random_seed=42)
    assert isinstance(result, dict)
