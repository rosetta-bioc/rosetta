"""MA plot for differential expression results."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ._detect import detect_columns

if TYPE_CHECKING:
    from matplotlib.figure import Figure


def ma_plot(
    results: pd.DataFrame,
    alpha: float = 0.05,
    title: str | None = None,
    ax=None,
    **kwargs,
) -> "Figure":
    """Create an MA plot from differential expression results.

    Plots mean expression (A) vs log fold change (M). Significant genes
    (adjusted p-value < alpha) are colored red; others grey.

    Parameters
    ----------
    results : pd.DataFrame
        Differential expression results (DESeq2, edgeR, or limma format).
    alpha : float
        Significance threshold for adjusted p-value.
    title : str, optional
        Plot title. Auto-generated if None.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. Creates new figure if None.
    **kwargs
        Additional keyword arguments passed to ax.scatter().

    Returns
    -------
    matplotlib.figure.Figure
        The figure containing the MA plot.

    Raises
    ------
    ValueError
        If no mean expression column is found in the results.
    """
    import matplotlib.pyplot as plt

    detected = detect_columns(results)
    pval_col = detected.pvalue
    lfc_col = detected.lfc
    mean_col = detected.mean_expr

    if mean_col is None:
        raise ValueError(
            f"No mean expression column found for {detected.tool} results. "
            "MA plots require baseMean (DESeq2), logCPM (edgeR), or AveExpr (limma)."
        )

    # Work on a copy, drop NaN p-values
    df = results[[lfc_col, pval_col, mean_col]].dropna().copy()

    # For DESeq2, log-transform baseMean for better visualization
    if detected.tool == "deseq2":
        a_values = np.log10(df[mean_col].clip(lower=1e-1))
        x_label = "log10(baseMean)"
    else:
        a_values = df[mean_col]
        x_label = mean_col

    sig = df[pval_col] < alpha
    non_sig = ~sig

    # Create axes if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.get_figure()

    scatter_kw = {"s": 8, "alpha": 0.6, "edgecolors": "none"}
    scatter_kw.update(kwargs)

    ax.scatter(
        a_values[non_sig],
        df.loc[non_sig, lfc_col],
        c="grey",
        label="Not significant",
        **scatter_kw,
    )
    ax.scatter(
        a_values[sig],
        df.loc[sig, lfc_col],
        c="red",
        label=f"Significant (n={sig.sum()})",
        **scatter_kw,
    )

    # Horizontal line at zero
    ax.axhline(0, color="blue", linestyle="-", linewidth=0.8)

    # Labels
    ax.set_xlabel(x_label)
    ax.set_ylabel("Log2 Fold Change" if detected.tool == "deseq2" else "logFC")
    if title is None:
        title = f"MA Plot ({detected.tool})"
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.8)

    fig.tight_layout()
    return fig
