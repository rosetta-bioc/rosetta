# tests/stats/test_design.py
import numpy as np
import pytest
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