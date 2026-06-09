"""Design matrix and contrast management."""

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from .._errors import RosettaSecurityError

def build_contrast_matrix(colnames, contrast_str):
    """
    Build a contrast matrix for linear model comparisons.
    
    Args:
        colnames (list): Column names from the design matrix.
        contrast_str (str): Contrast formula (e.g., "conditionA - conditionB").
        
    Returns:
        rpy2.robjects.Matrix: Contrast matrix.
    """
    if not contrast_str:
        return None
    limma_pkg = importr("limma")
    
    try:
        # Convert list of strings to R character vector
        r_colnames = ro.StrVector(colnames)
        contrast_mat = limma_pkg.makeContrasts(contrast_str, levels=r_colnames)
        return contrast_mat
    except Exception as e:
        raise RosettaSecurityError(f"Failed to build contrast matrix '{contrast_str}': {e}")