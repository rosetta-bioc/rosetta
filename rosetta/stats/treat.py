# rosetta/stats/treat.py
from .._bridge import ACTIVE_BACKEND

# Conditionally import rpy2 components based on the active backend
if ACTIVE_BACKEND == "rpy2":
    from rpy2.robjects.packages import importr
else:
    importr = None

def run_treat(fit_object, lfc=1.0, trend=False):
    """
    Wrapper for limma::treat.
    Args:
        fit_object: The MArrayLM object from lmFit.
        lfc: The log-fold-change threshold (default: 1.0).
        trend: Logical, whether to allow for an empirical Bayes trend.
    """
    limma = importr("limma")
    return limma.treat(fit_object, lfc=lfc, trend=trend)