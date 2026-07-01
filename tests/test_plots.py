"""Tests for rosetta.plots module."""

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest
from matplotlib.figure import Figure

from rosetta.plots import volcano, ma_plot, pca
from rosetta.plots._detect import detect_columns, DetectedColumns


# --- Fixtures ---


@pytest.fixture
def deseq2_results():
    """Synthetic DESeq2-like results DataFrame."""
    np.random.seed(42)
    n = 500
    return pd.DataFrame(
        {
            "baseMean": np.random.exponential(500, n),
            "log2FoldChange": np.random.normal(0, 2, n),
            "lfcSE": np.random.uniform(0.1, 0.5, n),
            "stat": np.random.normal(0, 3, n),
            "pvalue": np.random.uniform(0, 1, n),
            "padj": np.random.uniform(0, 1, n),
        },
        index=[f"gene_{i}" for i in range(n)],
    )


@pytest.fixture
def edger_results():
    """Synthetic edgeR-like results DataFrame."""
    np.random.seed(43)
    n = 500
    return pd.DataFrame(
        {
            "logFC": np.random.normal(0, 2, n),
            "logCPM": np.random.normal(5, 2, n),
            "F": np.random.uniform(0, 20, n),
            "PValue": np.random.uniform(0, 1, n),
            "FDR": np.random.uniform(0, 1, n),
        },
        index=[f"gene_{i}" for i in range(n)],
    )


@pytest.fixture
def limma_results():
    """Synthetic limma-like results DataFrame."""
    np.random.seed(44)
    n = 500
    return pd.DataFrame(
        {
            "logFC": np.random.normal(0, 1.5, n),
            "AveExpr": np.random.normal(6, 2, n),
            "t": np.random.normal(0, 3, n),
            "P.Value": np.random.uniform(0, 1, n),
            "adj.P.Val": np.random.uniform(0, 1, n),
            "B": np.random.normal(0, 2, n),
        },
        index=[f"gene_{i}" for i in range(n)],
    )


@pytest.fixture
def count_matrix():
    """Synthetic count matrix (genes x samples)."""
    np.random.seed(45)
    return pd.DataFrame(
        np.random.negative_binomial(5, 0.1, size=(200, 6)),
        index=[f"gene_{i}" for i in range(200)],
        columns=["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"],
    )


@pytest.fixture
def sample_metadata():
    """Sample metadata for PCA coloring."""
    return pd.DataFrame(
        {"condition": ["control"] * 3 + ["treated"] * 3},
        index=["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"],
    )


# --- _detect tests ---


class TestDetectColumns:
    def test_detect_deseq2(self, deseq2_results):
        detected = detect_columns(deseq2_results)
        assert detected.pvalue == "padj"
        assert detected.lfc == "log2FoldChange"
        assert detected.mean_expr == "baseMean"
        assert detected.tool == "deseq2"

    def test_detect_edger(self, edger_results):
        detected = detect_columns(edger_results)
        assert detected.pvalue == "FDR"
        assert detected.lfc == "logFC"
        assert detected.mean_expr == "logCPM"
        assert detected.tool == "edger"

    def test_detect_limma(self, limma_results):
        detected = detect_columns(limma_results)
        assert detected.pvalue == "adj.P.Val"
        assert detected.lfc == "logFC"
        assert detected.mean_expr == "AveExpr"
        assert detected.tool == "limma"

    def test_detect_unknown_raises(self):
        df = pd.DataFrame({"x": [1, 2], "y": [3, 4]})
        with pytest.raises(ValueError, match="Cannot detect result type"):
            detect_columns(df)


# --- Volcano tests ---


class TestVolcano:
    def test_returns_figure(self, deseq2_results):
        fig = volcano(deseq2_results)
        assert isinstance(fig, Figure)
        matplotlib.pyplot.close(fig)

    def test_edger_results(self, edger_results):
        fig = volcano(edger_results, alpha=0.1, lfc_cutoff=0.5)
        assert isinstance(fig, Figure)
        matplotlib.pyplot.close(fig)

    def test_limma_results(self, limma_results):
        fig = volcano(limma_results)
        assert isinstance(fig, Figure)
        matplotlib.pyplot.close(fig)

    def test_custom_title(self, deseq2_results):
        fig = volcano(deseq2_results, title="My Volcano")
        ax = fig.axes[0]
        assert ax.get_title() == "My Volcano"
        matplotlib.pyplot.close(fig)

    def test_highlight_genes(self, deseq2_results):
        genes = ["gene_0", "gene_1", "gene_999"]  # gene_999 doesn't exist
        fig = volcano(deseq2_results, highlight_genes=genes)
        assert isinstance(fig, Figure)
        matplotlib.pyplot.close(fig)

    def test_custom_ax(self, deseq2_results):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        result_fig = volcano(deseq2_results, ax=ax)
        assert result_fig is fig
        plt.close(fig)


# --- MA plot tests ---


class TestMAPlot:
    def test_returns_figure(self, deseq2_results):
        fig = ma_plot(deseq2_results)
        assert isinstance(fig, Figure)
        matplotlib.pyplot.close(fig)

    def test_edger_results(self, edger_results):
        fig = ma_plot(edger_results, alpha=0.1)
        assert isinstance(fig, Figure)
        matplotlib.pyplot.close(fig)

    def test_limma_results(self, limma_results):
        fig = ma_plot(limma_results)
        assert isinstance(fig, Figure)
        matplotlib.pyplot.close(fig)

    def test_custom_title(self, deseq2_results):
        fig = ma_plot(deseq2_results, title="My MA")
        ax = fig.axes[0]
        assert ax.get_title() == "My MA"
        matplotlib.pyplot.close(fig)

    def test_no_mean_column_raises(self):
        df = pd.DataFrame({"padj": [0.01, 0.1], "log2FoldChange": [1.5, -0.3]})
        with pytest.raises(ValueError, match="No mean expression column"):
            ma_plot(df)


# --- PCA tests ---


class TestPCA:
    def test_returns_figure(self, count_matrix):
        fig = pca(count_matrix)
        assert isinstance(fig, Figure)
        matplotlib.pyplot.close(fig)

    def test_color_by_metadata(self, count_matrix, sample_metadata):
        fig = pca(count_matrix, metadata=sample_metadata, color_by="condition")
        assert isinstance(fig, Figure)
        matplotlib.pyplot.close(fig)

    def test_custom_title(self, count_matrix):
        fig = pca(count_matrix, title="My PCA")
        ax = fig.axes[0]
        assert ax.get_title() == "My PCA"
        matplotlib.pyplot.close(fig)

    def test_axes_labels_contain_variance(self, count_matrix):
        fig = pca(count_matrix)
        ax = fig.axes[0]
        assert "variance" in ax.get_xlabel().lower()
        assert "variance" in ax.get_ylabel().lower()
        matplotlib.pyplot.close(fig)

    def test_custom_ax(self, count_matrix):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        result_fig = pca(count_matrix, ax=ax)
        assert result_fig is fig
        plt.close(fig)
