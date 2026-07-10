"""Volcano plot for differential expression results."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ._detect import detect_columns

if TYPE_CHECKING:
    from matplotlib.figure import Figure


def volcano(
    results: pd.DataFrame,
    alpha: float = 0.05,
    lfc_cutoff: float = 1.0,
    title: str | None = None,
    ax=None,
    highlight_genes: list[str] | None = None,
    **kwargs,
) -> "Figure":
    """Create a volcano plot from differential expression results.

    Plots -log10(adjusted p-value) vs log fold change. Points are colored:
      - Grey: not significant
      - Blue: significantly downregulated (lfc < -lfc_cutoff & padj < alpha)
      - Red: significantly upregulated (lfc > lfc_cutoff & padj < alpha)

    Parameters
    ----------
    results : pd.DataFrame
        Differential expression results (DESeq2, edgeR, or limma format).
    alpha : float
        Significance threshold for adjusted p-value.
    lfc_cutoff : float
        Log fold change cutoff for coloring.
    title : str, optional
        Plot title. Auto-generated if None.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. Creates new figure if None.
    highlight_genes : list of str, optional
        Gene names (index values) to label on the plot.
    **kwargs
        Additional keyword arguments passed to ax.scatter().

    Returns
    -------
    matplotlib.figure.Figure
        The figure containing the volcano plot.
    """
    import matplotlib.pyplot as plt

    detected = detect_columns(results)
    pval_col = detected.pvalue
    lfc_col = detected.lfc

    # Work on a copy, drop NaN p-values
    df = results[[lfc_col, pval_col]].dropna().copy()
    df["neg_log10_p"] = -np.log10(df[pval_col].clip(lower=1e-300))

    # Classify significance
    sig_up = (df[pval_col] < alpha) & (df[lfc_col] > lfc_cutoff)
    sig_down = (df[pval_col] < alpha) & (df[lfc_col] < -lfc_cutoff)
    non_sig = ~(sig_up | sig_down)

    # Create axes if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.get_figure()

    scatter_kw = {"s": 8, "alpha": 0.6, "edgecolors": "none"}
    scatter_kw.update(kwargs)

    # Plot each category
    ax.scatter(
        df.loc[non_sig, lfc_col],
        df.loc[non_sig, "neg_log10_p"],
        c="grey",
        label="Not significant",
        **scatter_kw,
    )
    ax.scatter(
        df.loc[sig_down, lfc_col],
        df.loc[sig_down, "neg_log10_p"],
        c="blue",
        label=f"Down (n={sig_down.sum()})",
        **scatter_kw,
    )
    ax.scatter(
        df.loc[sig_up, lfc_col],
        df.loc[sig_up, "neg_log10_p"],
        c="red",
        label=f"Up (n={sig_up.sum()})",
        **scatter_kw,
    )

    # Reference lines
    ax.axhline(-np.log10(alpha), color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(-lfc_cutoff, color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(lfc_cutoff, color="grey", linestyle="--", linewidth=0.8)

    # Labels
    ax.set_xlabel("Log2 Fold Change" if detected.tool == "deseq2" else "logFC")
    ax.set_ylabel(f"-log10({pval_col})")
    if title is None:
        title = f"Volcano Plot ({detected.tool})"
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.8)

    # Highlight specific genes
    if highlight_genes:
        for gene in highlight_genes:
            if gene in df.index:
                x = df.loc[gene, lfc_col]
                y = df.loc[gene, "neg_log10_p"]
                ax.annotate(
                    gene,
                    (x, y),
                    fontsize=7,
                    ha="center",
                    va="bottom",
                    arrowprops={"arrowstyle": "-", "color": "black", "lw": 0.5},
                )

    fig.tight_layout()
    return fig
