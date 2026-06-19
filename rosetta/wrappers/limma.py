"""limma-voom differential expression wrapper."""

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from .._bridge import _converter, to_r_matrix, to_r_dataframe, to_pandas, r_nrow
from .._deps import ensure_installed
from .._errors import RDataError, RFormulaError
from ..stats.design import build_contrast_matrix
from ..stats.decide import run_decide_tests


def limma_voom(counts: pd.DataFrame, metadata: pd.DataFrame, design: str = "~ condition", contrast=None, decide_tests=False, **kwargs) -> pd.DataFrame:
    """Run limma-voom differential expression analysis.

    Args:
        counts: Gene count matrix (genes x samples) with non-negative integers.
        metadata: Sample metadata DataFrame with row names matching counts columns.
        design: R formula string for the experimental design.
        **kwargs: Additional arguments passed to limma::lmFit().

    Returns:
        DataFrame with logFC, AveExpr, t, P.Value, adj.P.Val, B.
    """
    if (counts < 0).any().any():
        raise RDataError("Count matrix contains negative values")
    if not set(counts.columns).issubset(set(metadata.index)):
        raise RDataError("Count matrix columns must match metadata row names")

    stats_pkg = importr("stats")
    with localconverter(_converter):
        try:
            stats_pkg.as_formula(design)
        except Exception as e:
            raise RFormulaError(f"Invalid design formula '{design}': {e}") from e

    ensure_installed("limma")
    ensure_installed("edgeR")

    limma_pkg = importr("limma")
    edger_pkg = importr("edgeR")

    r_counts = to_r_matrix(counts)
    r_metadata = to_r_dataframe(metadata)

    with localconverter(_converter):
        r_design_matrix = stats_pkg.model_matrix(ro.Formula(design), data=r_metadata)
        # Use edgeR::voomLmFit (edgeR v4) — combines voom + lmFit in one
        # optimized step. Recommended by Gordon Smyth over separate
        # voom() + lmFit() calls.
        dge = edger_pkg.DGEList(counts=r_counts)
        dge = edger_pkg.calcNormFactors(dge)
        fit = edger_pkg.voomLmFit(dge, r_design_matrix, **kwargs)

        # Handle contrast matrix using stats module
        if contrast:
            contrast_mat = build_contrast_matrix(r_design_matrix.colnames, contrast)
            fit = limma_pkg.contrasts_fit(fit, contrast_mat)
        
        # Empirical Bayes moderation
        fit = limma_pkg.eBayes(fit)
        
        # Determine significant genes using stats module
        if decide_tests:
            fit.results = run_decide_tests(fit)

        # Extract results
        r_df = limma_pkg.topTable(fit, number=r_nrow(r_counts))

    return to_pandas(r_df)
