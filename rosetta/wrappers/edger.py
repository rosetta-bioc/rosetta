"""edgeR differential expression wrapper."""

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from .._bridge import _converter, to_r_matrix, to_r_dataframe, to_pandas, to_r_df, r_nrow
from .._deps import ensure_installed
from .._errors import RDataError, RFormulaError


def edger(counts: pd.DataFrame, metadata: pd.DataFrame, design: str = "~ condition", contrast=None, lfc: float = 0, **kwargs) -> pd.DataFrame:
    """Run edgeR quasi-likelihood differential expression analysis.

    Args:
        counts: Gene count matrix (genes x samples) with non-negative integers.
        metadata: Sample metadata DataFrame with row names matching counts columns.
        design: R formula string for the experimental design.
        **kwargs: Additional arguments passed to edgeR::glmQLFit().

    Returns:
        DataFrame with logFC, logCPM, F, PValue, FDR.
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

    ensure_installed("edgeR")
    edger_pkg = importr("edgeR")

    r_counts = to_r_matrix(counts)
    r_metadata = to_r_dataframe(metadata)

    with localconverter(_converter):
        r_design_matrix = stats_pkg.model_matrix(ro.Formula(design), data=r_metadata)
        dge = edger_pkg.DGEList(counts=r_counts)
        dge = edger_pkg.calcNormFactors(dge)
        dge = edger_pkg.estimateDisp(dge, r_design_matrix)
        fit = edger_pkg.glmQLFit(dge, r_design_matrix, **kwargs)

        # Prepare contrast object if provided
        r_contrast = None
        if contrast is not None:
            # Handle both list (e.g., [0, 1]) and DataFrame (matrix) inputs
            r_contrast = ro.FloatVector(contrast) if isinstance(contrast, list) else to_r_matrix(contrast)

        # Determine testing method (TREAT vs. QL-test)
        if lfc > 0:
            # Prepare arguments for glmTreat, excluding 'fit' from kwargs
            treat_args = {"lfc": lfc}
            if r_contrast is not None:
                treat_args["contrast"] = r_contrast
            
            # Pass 'fit' as a positional argument explicitly
            res = edger_pkg.glmTreat(fit, **treat_args)
        else:
            # Prepare arguments for glmQLFTest, excluding 'fit' from kwargs
            qlf_args = {}
            if r_contrast is not None:
                qlf_args["contrast"] = r_contrast
            
            # Pass 'fit' as a positional argument explicitly
            res = edger_pkg.glmQLFTest(fit, **qlf_args)

        top = edger_pkg.topTags(res, n=r_nrow(r_counts))

    return to_pandas(to_r_df(top))
