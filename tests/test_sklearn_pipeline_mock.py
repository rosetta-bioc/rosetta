"""Mock-based tests for sklearn_pipeline.py — no live R or sklearn needed."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from unittest.mock import MagicMock, patch, call


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _X(n_samples=6, n_genes=10, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        rng.integers(0, 500, (n_samples, n_genes)),
        index=[f"s{i}" for i in range(n_samples)],
        columns=[f"gene{i}" for i in range(n_genes)],
    )


def _meta(n_samples=6):
    return pd.DataFrame(
        {"condition": ["A"] * 3 + ["B"] * 3},
        index=[f"s{i}" for i in range(n_samples)],
    )


# =========================================================================
# DESeq2Transformer
# =========================================================================

def _make_deseq2_transformer(meta=None, **kwargs):
    import rosetta.sklearn_pipeline as mod
    meta = meta if meta is not None else _meta()

    mock_deseq2 = MagicMock()
    mock_stats  = MagicMock()
    mock_ensure = MagicMock()

    # importr("DESeq2") → mock_deseq2, importr("stats") → mock_stats
    def _importr(name):
        return mock_deseq2 if name == "DESeq2" else mock_stats

    # results() returns a mock; to_pandas turns it into a real df with padj
    sig_df = pd.DataFrame(
        {"padj": [0.01, 0.9, 0.02, 0.7, 0.03, 0.8, 0.05, 0.4, 0.9, 0.8]},
        index=[f"gene{i}" for i in range(10)],
    )

    ctx = [
        patch.object(mod, "importr", side_effect=_importr),
        patch("rosetta.sklearn_pipeline.ensure_installed", mock_ensure),
        patch("rosetta.sklearn_pipeline.to_r_matrix", return_value=MagicMock()),
        patch("rosetta.sklearn_pipeline.to_r_df", return_value=MagicMock()),
        patch("rosetta.sklearn_pipeline.to_pandas", return_value=sig_df),
        patch("rosetta.sklearn_pipeline.localconverter"),
        patch("rosetta.sklearn_pipeline._converter", MagicMock()),
    ]
    return mod.DESeq2Transformer(meta, **kwargs), ctx


def test_deseq2_transformer_stores_params():
    from rosetta.sklearn_pipeline import DESeq2Transformer
    meta = _meta()
    t = DESeq2Transformer(meta, design="~ batch", alpha=0.1, lfc_threshold=1.0)
    assert t.design == "~ batch"
    assert t.alpha == 0.1
    assert t.lfc_threshold == 1.0
    assert t.metadata is meta


def test_deseq2_transformer_fit_returns_self():
    t, ctx = _make_deseq2_transformer()
    with _stack(ctx):
        result = t.fit(_X())
    assert result is t


def test_deseq2_transformer_fit_sets_significant_genes():
    t, ctx = _make_deseq2_transformer(alpha=0.05)
    with _stack(ctx):
        t.fit(_X())
    # padj < 0.05: gene0 (0.01), gene2 (0.02), gene4 (0.03)
    assert set(t.significant_genes_) == {"gene0", "gene2", "gene4"}


def test_deseq2_transformer_fit_negative_counts_raises():
    from rosetta._errors import RDataError
    from rosetta.sklearn_pipeline import DESeq2Transformer
    t = DESeq2Transformer(_meta())
    X = _X()
    X.iloc[0, 0] = -1
    with patch("rosetta.sklearn_pipeline.ensure_installed"):
        with pytest.raises(RDataError, match="negative"):
            t.fit(X)


def test_deseq2_transformer_fit_non_df_raises():
    from rosetta.sklearn_pipeline import DESeq2Transformer
    t = DESeq2Transformer(_meta())
    with patch("rosetta.sklearn_pipeline.ensure_installed"):
        with pytest.raises(TypeError, match="pandas DataFrame"):
            t.fit(np.zeros((6, 10)))


def test_deseq2_transformer_transform_filters_genes():
    t, ctx = _make_deseq2_transformer(alpha=0.05)
    with _stack(ctx):
        t.fit(_X())
    X_new = _X()
    out = t.transform(X_new)
    assert list(out.columns) == sorted(out.columns.tolist())  # subset
    assert set(out.columns).issubset(set(X_new.columns))
    assert set(out.columns) == {"gene0", "gene2", "gene4"}


def test_deseq2_transformer_transform_no_overlap_raises():
    t, ctx = _make_deseq2_transformer(alpha=0.05)
    with _stack(ctx):
        t.fit(_X())
    X_unrelated = pd.DataFrame(
        np.zeros((6, 3)), columns=["other1", "other2", "other3"]
    )
    with pytest.raises(ValueError, match="None of the significant"):
        t.transform(X_unrelated)


def test_deseq2_transformer_not_fitted_raises():
    from sklearn.exceptions import NotFittedError
    from rosetta.sklearn_pipeline import DESeq2Transformer
    t = DESeq2Transformer(_meta())
    with pytest.raises(NotFittedError):
        t.transform(_X())


# =========================================================================
# EdgeRTransformer
# =========================================================================

def _make_edger_transformer(meta=None, fdr_col="FDR", fdr_values=None, **kwargs):
    import rosetta.sklearn_pipeline as mod
    meta = meta if meta is not None else _meta()

    mock_edger = MagicMock()
    mock_stats  = MagicMock()
    mock_ro     = MagicMock()

    def _importr(name):
        return mock_edger if name == "edgeR" else mock_stats

    if fdr_values is None:
        fdr_values = [0.01, 0.9, 0.03, 0.8, 0.5, 0.02, 0.7, 0.04, 0.6, 0.9]

    sig_df = pd.DataFrame(
        {fdr_col: fdr_values},
        index=[f"gene{i}" for i in range(10)],
    )

    ctx = [
        patch.object(mod, "importr", side_effect=_importr),
        patch("rosetta.sklearn_pipeline.ensure_installed"),
        patch("rosetta.sklearn_pipeline.to_r_matrix", return_value=MagicMock()),
        patch("rosetta.sklearn_pipeline.to_r_df", return_value=MagicMock()),
        patch("rosetta.sklearn_pipeline.to_pandas", return_value=sig_df),
        patch("rosetta.sklearn_pipeline.localconverter"),
        patch("rosetta.sklearn_pipeline._converter", MagicMock()),
        patch("rpy2.robjects", mock_ro),
    ]
    return mod.EdgeRTransformer(meta, **kwargs), ctx


def test_edger_transformer_stores_params():
    from rosetta.sklearn_pipeline import EdgeRTransformer
    meta = _meta()
    t = EdgeRTransformer(meta, alpha=0.1, filter_low_counts=False)
    assert t.alpha == 0.1
    assert t.filter_low_counts is False


def test_edger_transformer_fit_returns_self():
    t, ctx = _make_edger_transformer()
    with _stack(ctx):
        result = t.fit(_X())
    assert result is t


def test_edger_transformer_fit_sets_significant_genes():
    t, ctx = _make_edger_transformer(alpha=0.05)
    with _stack(ctx):
        t.fit(_X())
    # FDR < 0.05: gene0(0.01), gene2(0.03), gene5(0.02), gene7(0.04)
    assert set(t.significant_genes_) == {"gene0", "gene2", "gene5", "gene7"}


def test_edger_transformer_fit_negative_counts_raises():
    from rosetta._errors import RDataError
    from rosetta.sklearn_pipeline import EdgeRTransformer
    t = EdgeRTransformer(_meta())
    X = _X()
    X.iloc[0, 0] = -1
    with patch("rosetta.sklearn_pipeline.ensure_installed"):
        with pytest.raises(RDataError, match="negative"):
            t.fit(X)


def test_edger_transformer_transform_filters_genes():
    t, ctx = _make_edger_transformer(alpha=0.05)
    with _stack(ctx):
        t.fit(_X())
    out = t.transform(_X())
    assert set(out.columns) == {"gene0", "gene2", "gene5", "gene7"}


def test_edger_transformer_transform_no_overlap_raises():
    t, ctx = _make_edger_transformer(alpha=0.05)
    with _stack(ctx):
        t.fit(_X())
    X_other = pd.DataFrame(np.zeros((6, 2)), columns=["a", "b"])
    with pytest.raises(ValueError, match="None of the significant"):
        t.transform(X_other)


def test_edger_transformer_not_fitted_raises():
    from sklearn.exceptions import NotFittedError
    from rosetta.sklearn_pipeline import EdgeRTransformer
    t = EdgeRTransformer(_meta())
    with pytest.raises(NotFittedError):
        t.transform(_X())


def test_edger_transformer_fallback_pvalue_col():
    """When 'FDR' absent, falls back to 'PValue'."""
    t, ctx = _make_edger_transformer(
        fdr_col="PValue",
        fdr_values=[0.01, 0.9, 0.03, 0.8, 0.5, 0.02, 0.7, 0.9, 0.6, 0.9],
        alpha=0.05,
    )
    with _stack(ctx):
        t.fit(_X())
    assert set(t.significant_genes_) == {"gene0", "gene2", "gene5"}


# ---------------------------------------------------------------------------
# helper: stack context managers
# ---------------------------------------------------------------------------

from contextlib import contextmanager, ExitStack

@contextmanager
def _stack(ctx_managers):
    with ExitStack() as stack:
        for cm in ctx_managers:
            stack.enter_context(cm)
        yield
