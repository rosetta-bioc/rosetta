"""DESeq2 differential expression wrapper."""

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from .._bridge import _converter, to_r_matrix, to_pandas, to_r_df
from .._deps import ensure_installed
from .._errors import RDataError, RFormulaError


_SHRINK_METHODS = {"apeglm", "ashr", "normal"}


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


def get_results(dds, contrast: list = None, lfc_threshold: float = 0.0, alpha: float = 0.1, shrink: str | None = None) -> pd.DataFrame:
    """
    Extract results from a fitted DESeqDataSet.

    Args:
        dds: The fitted R DESeqDataSet object (from run_deseq2).
        contrast: A list defining the comparison [factor, numerator, denominator].
        lfc_threshold: Log2 fold change threshold (absolute).
        alpha: Significance cutoff (FDR).
        shrink: LFC shrinkage method — one of 'apeglm', 'ashr', or 'normal'. None skips shrinkage.
    """
    if shrink is not None and shrink not in _SHRINK_METHODS:
        raise ValueError(f"shrink must be one of {sorted(_SHRINK_METHODS)}, got '{shrink}'")

    deseq2_pkg = importr("DESeq2")

    with localconverter(_converter):
        args = {"alpha": alpha, "lfcThreshold": lfc_threshold}

        if contrast:
            args["contrast"] = ro.StrVector(contrast)

        try:
            if shrink is None:
                res = deseq2_pkg.results(dds, **args)
            else:
                if shrink in ("apeglm", "ashr"):
                    ensure_installed(shrink)
                coef_names = deseq2_pkg.resultsNames(dds)
                coef_name = coef_names[len(coef_names) - 1]
                res = deseq2_pkg.lfcShrink(dds, coef=coef_name, type=shrink)
            return to_pandas(to_r_df(res))
        except Exception as e:
            raise RDataError(f"Failed to extract results: {e}") from e


def deseq2(counts: pd.DataFrame, metadata: pd.DataFrame, design: str = "~ condition", shrink: str | None = None, **kwargs) -> pd.DataFrame:
    """Run DESeq2 differential expression analysis.

    Args:
        counts: Gene count matrix (genes x samples) with non-negative integers.
        metadata: Sample metadata DataFrame with row names matching counts columns.
        design: R formula string for the experimental design.
        shrink: LFC shrinkage method — one of 'apeglm', 'ashr', or 'normal'. None skips shrinkage.
        **kwargs: Additional arguments passed to DESeq2::results().

    Returns:
        DataFrame with baseMean, log2FoldChange, lfcSE, stat, pvalue, padj.
    """
    if shrink is not None and shrink not in _SHRINK_METHODS:
        raise ValueError(f"shrink must be one of {sorted(_SHRINK_METHODS)}, got '{shrink}'")

    dds = run_deseq2(counts, metadata, design)
    return get_results(dds, shrink=shrink)
