# tests/stats/test_design.py
import numpy as np
import pytest
from unittest.mock import MagicMock, patch
from rosetta._deps import is_installed
from rosetta.stats.design import build_contrast_matrix

@pytest.mark.skipif(not is_installed("limma"), reason="limma R package not installed")
def test_build_contrast_matrix_with_fixtures(sample_metadata):
    # 1. Prepare data
    # Manually define column names consistent with design matrix logic.
    # Based on the 'condition' column in sample_metadata, R's model.matrix 
    # typically generates these specific column names.
    colnames = ["(Intercept)", "conditiontreated"] 
    contrast_str = "conditiontreated - (Intercept)" # A standard contrast format
    
    # 2. Execute the function
    contrast_mat = build_contrast_matrix(colnames, contrast_str)
    
    # 3. Validation
    assert contrast_mat is not None
    
    # Convert R object to numpy array for numerical validation
    mat_np = np.array(contrast_mat)
    
    # Expected: (Intercept)=-1, conditiontreated=1
    # Verify dimensions are correct
    assert mat_np.shape == (2, 1)
    assert mat_np[0, 0] == -1
    assert mat_np[1, 0] == 1


# --- mock-based tests (no live R session required) ---

def test_build_contrast_matrix_empty_string_returns_none():
    """Empty contrast_str → None, no R call made."""
    assert build_contrast_matrix(["A", "B"], "") is None
    assert build_contrast_matrix(["A", "B"], None) is None


def test_build_contrast_matrix_importr_none_raises():
    """When importr is None (subprocess backend), RuntimeError is raised."""
    import rosetta.stats.design as design_mod
    orig = design_mod.importr
    design_mod.importr = None
    try:
        with pytest.raises(RuntimeError, match="rpy2 is required"):
            build_contrast_matrix(["A", "B"], "A - B")
    finally:
        design_mod.importr = orig


def test_build_contrast_matrix_mock_success():
    """Happy path: mocked limma.makeContrasts returns a sentinel matrix."""
    import rosetta.stats.design as design_mod

    sentinel = object()
    mock_limma = MagicMock()
    mock_limma.makeContrasts.return_value = sentinel

    mock_ro = MagicMock()
    mock_ro.StrVector.return_value = ["A", "B"]

    with patch.object(design_mod, "importr", return_value=mock_limma), \
         patch.object(design_mod, "ro", mock_ro):
        result = build_contrast_matrix(["A", "B"], "A - B")

    assert result is sentinel
    mock_limma.makeContrasts.assert_called_once()


def test_build_contrast_matrix_r_error_raises_security_error():
    """R-side exception is re-raised as RosettaSecurityError."""
    import rosetta.stats.design as design_mod
    from rosetta._errors import RosettaSecurityError

    mock_limma = MagicMock()
    mock_limma.makeContrasts.side_effect = RuntimeError("R exploded")

    mock_ro = MagicMock()
    mock_ro.StrVector.return_value = ["A", "B"]

    with patch.object(design_mod, "importr", return_value=mock_limma), \
         patch.object(design_mod, "ro", mock_ro):
        with pytest.raises(RosettaSecurityError, match="Failed to build contrast matrix"):
            build_contrast_matrix(["A", "B"], "A - B")