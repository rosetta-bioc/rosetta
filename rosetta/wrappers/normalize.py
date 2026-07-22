"""Normalization and transformation wrappers (DESeq2 VST/rlog, edgeR TMM)."""

import pandas as pd
from .._bridge import ACTIVE_BACKEND, _converter, to_r_matrix, to_r_dataframe, to_r_df
from .._deps import ensure_installed
from .._errors import RDataError

# Conditionally import rpy2 components based on the active backend
if ACTIVE_BACKEND == "rpy2":
    import rpy2.robjects as ro
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr
else:
    ro = None
    localconverter = None
    importr = None


def vst(
    counts: pd.DataFrame,
    metadata: pd.DataFrame = None,
    design: str = "~ 1",
    blind: bool = True,
) -> pd.DataFrame:
    """Variance stabilizing transformation via DESeq2.

    Returns transformed count matrix as pandas DataFrame (genes x samples).
    If blind=True (default), transformation is blind to experimental design.

    Args:
        counts: Raw count matrix (genes x samples), non-negative integers.
        metadata: Sample metadata DataFrame. If None, a dummy frame is created.
        design: R formula string for the experimental design.
        blind: If True, transformation ignores the design (recommended for QC).

    Returns:
        pd.DataFrame of transformed values with same shape/index as input.
    """
    if counts.empty:
        raise RDataError("Count matrix is empty")
    if (counts < 0).any().any():
        raise RDataError("Count matrix contains negative values")

    ensure_installed("DESeq2")
    deseq2_pkg = importr("DESeq2")
    stats_pkg = importr("stats")
    base_pkg = importr("base")

    if metadata is None:
        metadata = pd.DataFrame(
            {"intercept": [1] * counts.shape[1]}, index=counts.columns
        )

    r_counts = to_r_matrix(counts)
    r_metadata = to_r_dataframe(metadata)

    with localconverter(_converter):
        r_design = stats_pkg.as_formula(design)
        dds = deseq2_pkg.DESeqDataSetFromMatrix(
            countData=r_counts, colData=r_metadata, design=r_design
        )
        vsd = deseq2_pkg.varianceStabilizingTransformation(dds, blind=blind)
        mat = ro.r["assay"](vsd)
        result = ro.conversion.get_conversion().rpy2py(base_pkg.as_data_frame(mat))

    result.index = counts.index
    result.columns = counts.columns
    return pd.DataFrame(result)


def rlog(
    counts: pd.DataFrame,
    metadata: pd.DataFrame = None,
    design: str = "~ 1",
    blind: bool = True,
) -> pd.DataFrame:
    """Regularized log transformation via DESeq2.

    Returns transformed count matrix as pandas DataFrame (genes x samples).
    Slower than vst() but more robust for small sample sizes (<20).

    Args:
        counts: Raw count matrix (genes x samples), non-negative integers.
        metadata: Sample metadata DataFrame. If None, a dummy frame is created.
        design: R formula string for the experimental design.
        blind: If True, transformation ignores the design (recommended for QC).

    Returns:
        pd.DataFrame of transformed values with same shape/index as input.
    """
    if counts.empty:
        raise RDataError("Count matrix is empty")
    if (counts < 0).any().any():
        raise RDataError("Count matrix contains negative values")

    ensure_installed("DESeq2")
    deseq2_pkg = importr("DESeq2")
    stats_pkg = importr("stats")
    base_pkg = importr("base")

    if metadata is None:
        metadata = pd.DataFrame(
            {"intercept": [1] * counts.shape[1]}, index=counts.columns
        )

    r_counts = to_r_matrix(counts)
    r_metadata = to_r_dataframe(metadata)

    with localconverter(_converter):
        r_design = stats_pkg.as_formula(design)
        dds = deseq2_pkg.DESeqDataSetFromMatrix(
            countData=r_counts, colData=r_metadata, design=r_design
        )
        rld = deseq2_pkg.rlog(dds, blind=blind)
        mat = ro.r["assay"](rld)
        result = ro.conversion.get_conversion().rpy2py(base_pkg.as_data_frame(mat))

    result.index = counts.index
    result.columns = counts.columns
    return pd.DataFrame(result)


def tmm_normalize(counts: pd.DataFrame) -> pd.DataFrame:
    """TMM normalization via edgeR.

    Returns log-CPM normalized matrix as pandas DataFrame (genes x samples).
    Uses calcNormFactors(method='TMM') then cpm(log=TRUE).

    Args:
        counts: Raw count matrix (genes x samples), non-negative values.

    Returns:
        pd.DataFrame of log-CPM normalized values with same shape/index as input.
    """
    if counts.empty:
        raise RDataError("Count matrix is empty")
    if (counts < 0).any().any():
        raise RDataError("Count matrix contains negative values")

    ensure_installed("edgeR")
    edger_pkg = importr("edgeR")
    base_pkg = importr("base")

    r_counts = to_r_matrix(counts)

    with localconverter(_converter):
        dge = edger_pkg.DGEList(counts=r_counts)
        dge = edger_pkg.calcNormFactors(dge, method="TMM")
        log_cpm = edger_pkg.cpm(dge, log=True)
        result = ro.conversion.get_conversion().rpy2py(base_pkg.as_data_frame(log_cpm))

    result.index = counts.index
    result.columns = counts.columns
    return pd.DataFrame(result)
