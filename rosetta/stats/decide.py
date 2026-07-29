"""DecideTests integration for statistical significance."""

from .._bridge import ACTIVE_BACKEND

# Conditionally import rpy2 components based on the active backend
if ACTIVE_BACKEND == "rpy2":
    from rpy2.robjects.packages import importr
else:
    importr = None

def run_decide_tests(fit_object, method="separate", adj="BH", p_value=0.05):
    """
    Wrapper for decideTests() to determine significant genes.
    
    Args:
        fit_object: The fitted model object (from limma or edgeR).
        method: Method for multiple testing (default: 'separate').
        adj: Adjustment method (default: 'BH').
        p_value: Significance threshold.
        
    Returns:
        rpy2.robjects.Matrix: A matrix of 1, -1, 0 indicating significance.
    """
    # Determine the source based on the input object. Although both limma and edgeR 
    # provide decideTests, we provide a unified API here through this wrapper. 
    limma_pkg = importr("limma")
    
    try:
        results = limma_pkg.decideTests(
            fit_object, 
            method=method, 
            adjust_method=adj, 
            p_value=p_value
        )
        return results
    except Exception as e:
        raise RuntimeError(f"decideTests() failed for the provided fit object: {e}")

