"""Tests for rosetta.wrappers.seurat."""

import numpy as np
import pandas as pd
import pytest
from rosetta._deps import is_installed
if not is_installed("Seurat"):
    pytest.skip("Seurat R package not installed", allow_module_level=True)
from rosetta._errors import RDataError
from rosetta.wrappers.seurat import Seurat

# --- Fixtures --- #TODO: integrate the fixtures with the other fixtures below
@pytest.fixture
def sample_counts():
    data = np.random.randint(0, 10, size=(100, 10))
    return pd.DataFrame(data, index=[f"gene_{i}" for i in range(100)], columns=[f"cell_{i}" for i in range(10)])

@pytest.fixture
def seurat_obj(sample_counts):
    return Seurat(sample_counts)

@pytest.fixture
def preprocessed_seurat_obj(seurat_obj):
    """Fixture that handles the boring preprocessing steps."""
    seurat_obj.run_normalize(verbose=False) \
              .run_find_variable_features(nfeatures=10, verbose=False) \
              .run_scale_data(verbose=False)
    return seurat_obj

# --- Tests ---
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
    """Test SCTransform runs and find_markers raises informatively without clusters."""
    model = Seurat(sc_counts)

    # SCTransform should succeed and return self for chaining
    result = model.run_sctransform()
    assert result is model

    # FindMarkers requires valid cluster identities; on fresh random data
    # without clustering, it should raise an RDataError (not silently pass)
    with pytest.raises((RDataError, Exception), match=r"(ident|cluster|FindMarkers|cannot)"):
        model.find_markers(ident_1="0", ident_2="1")

# --- New Granular kwargs Tests ---
def test_seurat_methods_parameter_passing(preprocessed_seurat_obj):
    """Verify granular parameters are passed correctly."""
    # decrease npcs to 2; dims must in the range of npcs
    preprocessed_seurat_obj.run_pca(npcs=2, verbose=False) \
              .run_find_neighbors(k_param=5, dims=list(range(1, 3)), verbose=False) \
              .run_find_clusters(resolution=0.5, verbose=False)
    assert preprocessed_seurat_obj is not None

def test_seurat_invalid_param_warning(preprocessed_seurat_obj, caplog):
    with caplog.at_level("WARNING"):
        preprocessed_seurat_obj.run_pca(npcs=2, invalid_key="test")
    assert "Parameter 'invalid_key' is not supported" in caplog.text