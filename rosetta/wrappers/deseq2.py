"""DESeq2 differential expression wrapper."""

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from .._bridge import _converter, to_r_matrix, to_pandas, to_r_df
from .._deps import ensure_installed
from .._errors import RDataError, RFormulaError

# --- Private Helper Functions ---
def _prepare_dds(counts: pd.DataFrame, metadata: pd.DataFrame, design: str):
    """Internal helper to setup DESeqDataSet with error handling."""
    if (counts < 0).any().any():
        raise RDataError("Count matrix contains negative values")
    if not set(counts.columns).issubset(set(metadata.index)):
        raise RDataError("Count matrix columns must match metadata row names")

    ensure_installed("DESeq2")
    deseq2_pkg = importr("DESeq2")
    stats_pkg = importr("stats")
    
    try:
        r_design = stats_pkg.as_formula(design)
    except Exception as e:
        raise RFormulaError(f"Invalid design formula '{design}': {e}") from e

    r_counts = to_r_matrix(counts)
    r_metadata = to_r_df(metadata)

    with localconverter(_converter):
        try:
            return deseq2_pkg.DESeqDataSetFromMatrix(
                countData=r_counts, colData=r_metadata, design=r_design
            )
        except Exception as e:
            # Here we wrap the R error into our custom RDataError
            raise RDataError(f"Failed to create DESeqDataSet: {e}") from e

# --- Public API ---

def preview_design(counts: pd.DataFrame, metadata: pd.DataFrame, design: str = "~ condition"):
    """Lightweight preview: Initialize DESeqDataSet without fitting the model."""
    return _prepare_dds(counts, metadata, design)

def get_results_names(counts: pd.DataFrame, metadata: pd.DataFrame, design: str = "~ condition") -> list:
    """Parameter check: Fit the model and return available results names."""
    deseq2_pkg = importr("DESeq2")
    dds = _prepare_dds(counts, metadata, design)
    
    with localconverter(_converter):
        try:
            dds = deseq2_pkg.DESeq(dds)
            return list(deseq2_pkg.resultsNames(dds))
        except Exception as e:
            raise RDataError(f"DESeq2 analysis failed: {e}") from e

def run_deseq2(counts: pd.DataFrame, metadata: pd.DataFrame, design: str = "~ condition"):
    """Core operation: Perform full model fitting and return the dds object."""
    deseq2_pkg = importr("DESeq2")
    dds = _prepare_dds(counts, metadata, design)
    with localconverter(_converter):
        return deseq2_pkg.DESeq(dds)

def deseq2(counts: pd.DataFrame, metadata: pd.DataFrame, design: str = "~ condition", **kwargs) -> pd.DataFrame:
    """Legacy API: Compatibility wrapper to run the full pipeline."""
    deseq2_pkg = importr("DESeq2")
    dds = run_deseq2(counts, metadata, design)
    
    with localconverter(_converter):
        res = deseq2_pkg.results(dds, **kwargs)
        
    return to_pandas(to_r_df(res))


def get_results(dds, contrast: list = None, lfc_threshold: float = 0.0, alpha: float = 0.1) -> pd.DataFrame:
    """
    Extract results from a fitted DESeqDataSet.
    
    Args:
        dds: The fitted R DESeqDataSet object (from run_deseq2).
        contrast: A list defining the comparison [factor, numerator, denominator].
        lfc_threshold: Log2 fold change threshold (absolute).
        alpha: Significance cutoff (FDR).
    """
    deseq2_pkg = importr("DESeq2")
    
    with localconverter(_converter):
        # Build arguments dictionary for results()
        args = {"alpha": alpha, "lfcThreshold": lfc_threshold}
        
        # Handle contrast if provided
        if contrast:
            # R expects a character vector for contrast: c(factor, num, den)
            args["contrast"] = ro.StrVector(contrast)
            
        try:
            res = deseq2_pkg.results(dds, **args)
            return to_pandas(to_r_df(res))
        except Exception as e:
            raise RDataError(f"Failed to extract results: {e}") from e