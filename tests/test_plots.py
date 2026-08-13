"""Tests for rosetta.plots — MA, PCA, volcano, column detection."""
import numpy as np
import pandas as pd
import pytest

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for CI


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _deseq2_df(n=20, n_sig=5):
    rng = np.random.default_rng(0)
    df = pd.DataFrame({
        "baseMean": rng.uniform(10, 1000, n),
        "log2FoldChange": rng.uniform(-3, 3, n),
        "padj": np.concatenate([rng.uniform(0.0, 0.04, n_sig),
                                 rng.uniform(0.1, 1.0, n - n_sig)]),
        "pvalue": rng.uniform(0, 1, n),
    }, index=[f"gene{i}" for i in range(n)])
    return df


def _edger_df(n=20):
    rng = np.random.default_rng(1)
    return pd.DataFrame({
        "logFC": rng.uniform(-3, 3, n),
        "logCPM": rng.uniform(1, 10, n),
        "FDR": rng.uniform(0, 1, n),
        "PValue": rng.uniform(0, 1, n),
    }, index=[f"gene{i}" for i in range(n)])


def _limma_df(n=20):
    rng = np.random.default_rng(2)
    return pd.DataFrame({
        "logFC": rng.uniform(-3, 3, n),
        "AveExpr": rng.uniform(1, 10, n),
        "adj.P.Val": rng.uniform(0, 1, n),
        "P.Value": rng.uniform(0, 1, n),
    }, index=[f"gene{i}" for i in range(n)])


# ---------------------------------------------------------------------------
# plots/_detect.py
# ---------------------------------------------------------------------------

def test_detect_deseq2():
    from rosetta.plots._detect import detect_columns
    det = detect_columns(_deseq2_df())
    assert det.tool == "deseq2"
    assert det.lfc == "log2FoldChange"
    assert det.pvalue == "padj"
    assert det.mean_expr == "baseMean"


def test_detect_edger():
    from rosetta.plots._detect import detect_columns
    det = detect_columns(_edger_df())
    assert det.tool == "edger"
    assert det.lfc == "logFC"
    assert det.mean_expr == "logCPM"


def test_detect_limma():
    from rosetta.plots._detect import detect_columns
    det = detect_columns(_limma_df())
    assert det.tool == "limma"
    assert det.mean_expr == "AveExpr"


def test_detect_unknown_raises():
    from rosetta.plots._detect import detect_columns
    with pytest.raises(ValueError, match="Cannot detect"):
        detect_columns(pd.DataFrame({"x": [1], "y": [2]}))


# ---------------------------------------------------------------------------
# plots/ma.py
# ---------------------------------------------------------------------------

def test_ma_plot_deseq2_returns_figure():
    from rosetta.plots.ma import ma_plot
    fig = ma_plot(_deseq2_df())
    assert fig is not None
    import matplotlib.pyplot as plt
    plt.close("all")


def test_ma_plot_edger_returns_figure():
    from rosetta.plots.ma import ma_plot
    fig = ma_plot(_edger_df())
    assert fig is not None
    import matplotlib.pyplot as plt
    plt.close("all")


def test_ma_plot_limma_returns_figure():
    from rosetta.plots.ma import ma_plot
    fig = ma_plot(_limma_df())
    assert fig is not None
    import matplotlib.pyplot as plt
    plt.close("all")


def test_ma_plot_custom_title():
    from rosetta.plots.ma import ma_plot
    fig = ma_plot(_deseq2_df(), title="My MA")
    ax = fig.axes[0]
    assert ax.get_title() == "My MA"
    import matplotlib.pyplot as plt
    plt.close("all")


def test_ma_plot_no_mean_col_raises():
    from rosetta.plots.ma import ma_plot
    df = pd.DataFrame({"log2FoldChange": [1.0], "padj": [0.01]})
    with pytest.raises(ValueError, match="No mean expression column"):
        ma_plot(df)


def test_ma_plot_uses_provided_ax():
    import matplotlib.pyplot as plt
    from rosetta.plots.ma import ma_plot
    fig, ax = plt.subplots()
    returned_fig = ma_plot(_deseq2_df(), ax=ax)
    assert returned_fig is fig
    plt.close("all")


# ---------------------------------------------------------------------------
# plots/volcano.py
# ---------------------------------------------------------------------------

def test_volcano_deseq2_returns_figure():
    from rosetta.plots.volcano import volcano
    fig = volcano(_deseq2_df())
    assert fig is not None
    import matplotlib.pyplot as plt
    plt.close("all")


def test_volcano_custom_title():
    from rosetta.plots.volcano import volcano
    fig = volcano(_deseq2_df(), title="My Volcano")
    ax = fig.axes[0]
    assert "My Volcano" in ax.get_title()
    import matplotlib.pyplot as plt
    plt.close("all")


def test_volcano_highlight_genes():
    from rosetta.plots.volcano import volcano
    df = _deseq2_df(n=20, n_sig=10)
    highlight = [df.index[0], df.index[1]]
    fig = volcano(df, highlight_genes=highlight)
    assert fig is not None
    import matplotlib.pyplot as plt
    plt.close("all")


# ---------------------------------------------------------------------------
# plots/pca.py
# ---------------------------------------------------------------------------

def _count_matrix(n_genes=50, n_samples=8):
    rng = np.random.default_rng(3)
    counts = pd.DataFrame(
        rng.integers(0, 500, size=(n_genes, n_samples)),
        index=[f"gene{i}" for i in range(n_genes)],
        columns=[f"sample{i}" for i in range(n_samples)],
    )
    return counts


def test_pca_returns_figure():
    from rosetta.plots.pca import pca
    fig = pca(_count_matrix())
    assert fig is not None
    import matplotlib.pyplot as plt
    plt.close("all")


def test_pca_with_metadata_color_by():
    import matplotlib.pyplot as plt
    from rosetta.plots.pca import pca
    counts = _count_matrix(n_samples=6)
    meta = pd.DataFrame(
        {"condition": ["A", "A", "A", "B", "B", "B"]},
        index=counts.columns,
    )
    fig = pca(counts, metadata=meta, color_by="condition")
    assert fig is not None
    plt.close("all")


def test_pca_custom_title():
    import matplotlib.pyplot as plt
    from rosetta.plots.pca import pca
    fig = pca(_count_matrix(), title="My PCA")
    ax = fig.axes[0]
    assert "My PCA" in ax.get_title()
    plt.close("all")
