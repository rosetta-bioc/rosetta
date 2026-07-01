"""Column auto-detection for DESeq2, edgeR, and limma results."""

from __future__ import annotations

from typing import NamedTuple

import pandas as pd


class DetectedColumns(NamedTuple):
    """Detected column names for plotting."""

    pvalue: str  # adjusted p-value column
    lfc: str  # log fold change column
    mean_expr: str | None  # mean expression column (for MA plots)
    tool: str  # "deseq2", "edger", or "limma"


def detect_columns(df: pd.DataFrame) -> DetectedColumns:
    """Auto-detect which tool produced the results based on column names.

    Detection order:
      1. DESeq2: padj, log2FoldChange, baseMean
      2. edgeR: FDR, logFC, logCPM
      3. limma: adj.P.Val, logFC, AveExpr

    Parameters
    ----------
    df : pd.DataFrame
        Results DataFrame from a differential expression analysis.

    Returns
    -------
    DetectedColumns
        Named tuple with (pvalue, lfc, mean_expr, tool).

    Raises
    ------
    ValueError
        If columns cannot be identified as any supported tool output.
    """
    cols = set(df.columns)

    # DESeq2
    if "padj" in cols and "log2FoldChange" in cols:
        mean_expr = "baseMean" if "baseMean" in cols else None
        return DetectedColumns(
            pvalue="padj",
            lfc="log2FoldChange",
            mean_expr=mean_expr,
            tool="deseq2",
        )

    # edgeR
    if "FDR" in cols and "logFC" in cols:
        mean_expr = "logCPM" if "logCPM" in cols else None
        return DetectedColumns(
            pvalue="FDR",
            lfc="logFC",
            mean_expr=mean_expr,
            tool="edger",
        )

    # limma
    if "adj.P.Val" in cols and "logFC" in cols:
        mean_expr = "AveExpr" if "AveExpr" in cols else None
        return DetectedColumns(
            pvalue="adj.P.Val",
            lfc="logFC",
            mean_expr=mean_expr,
            tool="limma",
        )

    raise ValueError(
        "Cannot detect result type. Expected columns from DESeq2 "
        "(padj, log2FoldChange), edgeR (FDR, logFC), or limma (adj.P.Val, logFC). "
        f"Found columns: {sorted(df.columns.tolist())}"
    )
