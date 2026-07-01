"""PCA plot for count matrices."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from matplotlib.figure import Figure


def _pca_numpy(data: np.ndarray, n_components: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """PCA via numpy SVD fallback.

    Parameters
    ----------
    data : np.ndarray
        Centered data matrix (samples x genes).
    n_components : int
        Number of components to return.

    Returns
    -------
    tuple of (scores, variance_explained_ratio)
    """
    # Center the data
    centered = data - data.mean(axis=0)
    # SVD
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)
    # Scores
    scores = U[:, :n_components] * S[:n_components]
    # Variance explained
    var_explained = (S**2) / (S**2).sum()
    return scores, var_explained[:n_components]


def _pca_sklearn(data: np.ndarray, n_components: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """PCA via sklearn."""
    from sklearn.decomposition import PCA as SklearnPCA

    model = SklearnPCA(n_components=n_components)
    scores = model.fit_transform(data)
    return scores, model.explained_variance_ratio_


def pca(
    counts: pd.DataFrame,
    metadata: pd.DataFrame | None = None,
    color_by: str | None = None,
    n_components: int = 2,
    title: str | None = None,
    ax=None,
    **kwargs,
) -> "Figure":
    """Create a PCA plot from a count matrix.

    Transforms counts with log2(counts + 1), then performs PCA and plots
    PC1 vs PC2 (or higher components if n_components > 2, but only the
    first two are plotted).

    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix with genes as rows and samples as columns.
    metadata : pd.DataFrame, optional
        Sample metadata with samples as rows. Index must match count columns.
    color_by : str, optional
        Column in metadata to color points by. Requires metadata.
    n_components : int
        Number of PCA components to compute (at least 2 for plotting).
    title : str, optional
        Plot title. Defaults to "PCA".
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. Creates new figure if None.
    **kwargs
        Additional keyword arguments passed to ax.scatter().

    Returns
    -------
    matplotlib.figure.Figure
        The figure containing the PCA plot.
    """
    import matplotlib
    import matplotlib.pyplot as plt

    # Log2 transform: samples as rows, genes as columns
    log_counts = np.log2(counts.values.T.astype(float) + 1)

    # Run PCA (prefer sklearn, fall back to numpy)
    try:
        scores, var_ratio = _pca_sklearn(log_counts, n_components=n_components)
    except ImportError:
        scores, var_ratio = _pca_numpy(log_counts, n_components=n_components)

    # Create axes if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.get_figure()

    scatter_kw = {"s": 50, "alpha": 0.8, "edgecolors": "black", "linewidths": 0.5}
    scatter_kw.update(kwargs)

    sample_names = counts.columns.tolist()

    if color_by and metadata is not None and color_by in metadata.columns:
        # Color by metadata column
        groups = metadata.loc[sample_names, color_by]
        unique_groups = groups.unique()
        cmap = matplotlib.colormaps.get_cmap("tab10")

        for i, group in enumerate(unique_groups):
            mask = groups == group
            ax.scatter(
                scores[mask, 0],
                scores[mask, 1],
                c=[cmap(i)],
                label=str(group),
                **scatter_kw,
            )
        ax.legend(title=color_by, fontsize=8, framealpha=0.8)
    else:
        ax.scatter(scores[:, 0], scores[:, 1], **scatter_kw)

    # Labels with variance explained
    ax.set_xlabel(f"PC1 ({var_ratio[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({var_ratio[1]*100:.1f}% variance)")
    if title is None:
        title = "PCA"
    ax.set_title(title)

    fig.tight_layout()
    return fig
