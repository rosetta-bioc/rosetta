"""Tests for rosetta.sklearn_compat — TransformerMixin-compatible DE wrappers."""

import numpy as np
import pandas as pd
import pytest

from rosetta.results import RosettaDataFrame
from rosetta.sklearn_compat import (
    DESeq2Transformer,
    EdgeRTransformer,
    LimmaTransformer,
    _BaseRosettaTransformer,
)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def deseq2_result():
    np.random.seed(0)
    n = 200
    idx = [f"GENE{i}" for i in range(n)]
    df = RosettaDataFrame(
        {
            "baseMean": np.random.exponential(100, n),
            "log2FoldChange": np.random.normal(0, 1.5, n),
            "lfcSE": np.abs(np.random.normal(0.3, 0.1, n)),
            "stat": np.random.normal(0, 2, n),
            "pvalue": np.random.uniform(0, 1, n),
            "padj": np.concatenate(
                [np.random.uniform(0, 0.01, 30), np.random.uniform(0.1, 1, n - 30)]
            ),
        },
        index=idx,
    )
    df._rosetta_method = "deseq2"
    return df


@pytest.fixture
def edger_result():
    np.random.seed(1)
    n = 150
    idx = [f"GENE{i}" for i in range(n)]
    df = RosettaDataFrame(
        {
            "logFC": np.random.normal(0, 2, n),
            "logCPM": np.random.normal(5, 1, n),
            "F": np.abs(np.random.normal(3, 2, n)),
            "PValue": np.random.uniform(0, 1, n),
            "FDR": np.concatenate(
                [np.random.uniform(0, 0.02, 20), np.random.uniform(0.1, 1, n - 20)]
            ),
        },
        index=idx,
    )
    df._rosetta_method = "edger"
    return df


@pytest.fixture
def sample_counts():
    np.random.seed(42)
    genes = [f"G{i}" for i in range(50)]
    samples = [f"S{i}" for i in range(6)]
    return pd.DataFrame(
        np.random.poisson(100, (50, 6)), index=genes, columns=samples
    )


@pytest.fixture
def sample_metadata():
    return pd.DataFrame(
        {"condition": ["ctrl", "ctrl", "ctrl", "treat", "treat", "treat"]},
        index=[f"S{i}" for i in range(6)],
    )


# ---------------------------------------------------------------------------
# Helper: a concrete transformer that uses pre-built results (no R required)
# ---------------------------------------------------------------------------

class _MockDETransformer(_BaseRosettaTransformer):
    """Concrete subclass for testing without calling R."""

    def __init__(self, mock_result, **kwargs):
        # We don't need real counts/metadata for unit tests
        kwargs.setdefault("counts", pd.DataFrame())
        kwargs.setdefault("metadata", pd.DataFrame())
        super().__init__(**kwargs)
        self._mock_result = mock_result

    def _run_de(self):
        return self._mock_result

    @property
    def _method_name(self):
        return "mock"


# ---------------------------------------------------------------------------
# _BaseRosettaTransformer unit tests
# ---------------------------------------------------------------------------

class TestBaseTransformer:
    def test_fit_sets_result(self, deseq2_result):
        tr = _MockDETransformer(deseq2_result)
        tr.fit()
        assert tr.result_ is not None
        assert isinstance(tr.result_, RosettaDataFrame)

    def test_fit_sets_feature_names(self, deseq2_result):
        tr = _MockDETransformer(deseq2_result)
        tr.fit()
        assert tr.feature_names_ is not None
        assert len(tr.feature_names_) == 30  # 30 sig genes in fixture

    def test_transform_shape(self, deseq2_result):
        # 6 samples in metadata
        meta = pd.DataFrame(index=[f"S{i}" for i in range(6)])
        tr = _MockDETransformer(deseq2_result, metadata=meta)
        tr.fit()
        out = tr.transform()
        assert out.shape == (6, 30)

    def test_transform_before_fit_raises(self, deseq2_result):
        tr = _MockDETransformer(deseq2_result)
        with pytest.raises(RuntimeError, match="fit()"):
            tr.transform()

    def test_get_feature_names_before_fit_raises(self, deseq2_result):
        tr = _MockDETransformer(deseq2_result)
        with pytest.raises(RuntimeError, match="fit()"):
            tr.get_feature_names_out()

    def test_get_feature_names_after_fit(self, deseq2_result):
        tr = _MockDETransformer(deseq2_result)
        tr.fit()
        names = tr.get_feature_names_out()
        assert isinstance(names, np.ndarray)
        assert all(n.startswith("GENE") for n in names)

    def test_select_all(self, deseq2_result):
        tr = _MockDETransformer(deseq2_result, select="all")
        tr.fit()
        assert len(tr.feature_names_) == 200  # all genes

    def test_lfc_threshold_filters(self, deseq2_result):
        tr_no_thresh = _MockDETransformer(deseq2_result)
        tr_thresh = _MockDETransformer(deseq2_result, lfc_threshold=1.0)
        tr_no_thresh.fit()
        tr_thresh.fit()
        assert len(tr_thresh.feature_names_) <= len(tr_no_thresh.feature_names_)

    def test_no_genes_raises(self, deseq2_result):
        tr = _MockDETransformer(deseq2_result, alpha=1e-99)
        with pytest.raises(ValueError, match="No genes passed filters"):
            tr.fit()

    def test_fit_transform_equivalent(self, deseq2_result):
        meta = pd.DataFrame(index=[f"S{i}" for i in range(4)])
        tr1 = _MockDETransformer(deseq2_result, metadata=meta)
        tr2 = _MockDETransformer(deseq2_result, metadata=meta)
        out1 = tr1.fit_transform()
        tr2.fit()
        out2 = tr2.transform()
        np.testing.assert_array_equal(out1, out2)

    def test_method_metadata_set_on_result(self, deseq2_result):
        tr = _MockDETransformer(deseq2_result)
        tr.fit()
        assert tr.result_._rosetta_method == "mock"

    def test_padj_col_deseq2(self, deseq2_result):
        tr = _MockDETransformer(deseq2_result)
        tr.fit()
        assert tr._padj_col() == "padj"

    def test_padj_col_edger(self, edger_result):
        tr = _MockDETransformer(edger_result)
        tr.fit()
        assert tr._padj_col() == "FDR"

    def test_value_col_fallback(self, edger_result):
        # EdgeR uses 'logFC'; asking for 'log2FoldChange' should fall back
        tr = _MockDETransformer(edger_result, value="log2FoldChange")
        tr.fit()
        assert tr._value_col() == "logFC"

    def test_stat_value_col(self, deseq2_result):
        tr = _MockDETransformer(deseq2_result, value="stat")
        tr.fit()
        assert tr._value_col() == "stat"
        out = tr.transform()
        assert out.shape[1] == 30


# ---------------------------------------------------------------------------
# get_params / set_params (sklearn clone compatibility)
# ---------------------------------------------------------------------------

class TestSklearnParamInterface:
    def test_get_params_keys(self, deseq2_result):
        tr = _MockDETransformer(deseq2_result)
        params = tr.get_params()
        for key in ("alpha", "lfc_threshold", "select", "value", "design"):
            assert key in params

    def test_set_params_updates_alpha(self, deseq2_result):
        tr = _MockDETransformer(deseq2_result)
        tr.set_params(alpha=0.01)
        assert tr.alpha == 0.01

    def test_set_params_returns_self(self, deseq2_result):
        tr = _MockDETransformer(deseq2_result)
        result = tr.set_params(alpha=0.01)
        assert result is tr


# ---------------------------------------------------------------------------
# DESeq2Transformer extra params
# ---------------------------------------------------------------------------

class TestDESeq2TransformerParams:
    def test_extra_params_in_get_params(self, sample_counts, sample_metadata):
        tr = DESeq2Transformer(sample_counts, sample_metadata, shrinkage="apeglm")
        params = tr.get_params()
        assert "shrinkage" in params
        assert params["shrinkage"] == "apeglm"
        assert "contrast" in params

    def test_default_params(self, sample_counts, sample_metadata):
        tr = DESeq2Transformer(sample_counts, sample_metadata)
        assert tr.alpha == 0.05
        assert tr.shrinkage is None
        assert tr.contrast is None
        assert tr.select == "significant"
        assert tr.value == "log2FoldChange"


# ---------------------------------------------------------------------------
# EdgeRTransformer / LimmaTransformer construction
# ---------------------------------------------------------------------------

class TestOtherTransformers:
    def test_edger_transformer_init(self, sample_counts, sample_metadata):
        tr = EdgeRTransformer(sample_counts, sample_metadata)
        assert tr._method_name == "edger"
        assert tr.alpha == 0.05

    def test_limma_transformer_init(self, sample_counts, sample_metadata):
        tr = LimmaTransformer(sample_counts, sample_metadata)
        assert tr._method_name == "limma"

    def test_edger_mock_fit_transform(self, edger_result):
        class _MockEdgeR(_MockDETransformer):
            @property
            def _method_name(self):
                return "edger"

        meta = pd.DataFrame(index=[f"S{i}" for i in range(3)])
        tr = _MockEdgeR(edger_result, metadata=meta)
        out = tr.fit_transform()
        assert out.shape[0] == 3
        assert out.shape[1] == 20  # 20 sig genes in edger fixture


# ---------------------------------------------------------------------------
# sklearn Pipeline smoke test (no R, uses mock)
# ---------------------------------------------------------------------------

class TestSklearnPipelineIntegration:
    def test_pipeline_with_pca(self, deseq2_result):
        pytest.importorskip("sklearn")
        from sklearn.decomposition import PCA
        from sklearn.pipeline import Pipeline

        meta = pd.DataFrame(index=[f"S{i}" for i in range(10)])
        tr = _MockDETransformer(deseq2_result, metadata=meta)

        pipe = Pipeline([
            ("de", tr),
            ("pca", PCA(n_components=5)),
        ])
        out = pipe.fit_transform(None)
        assert out.shape == (10, 5)

    def test_pipeline_fit_then_transform(self, deseq2_result):
        pytest.importorskip("sklearn")
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler

        meta = pd.DataFrame(index=[f"S{i}" for i in range(8)])
        tr = _MockDETransformer(deseq2_result, metadata=meta)

        pipe = Pipeline([
            ("de", tr),
            ("scale", StandardScaler()),
        ])
        pipe.fit(None)
        out = pipe.transform(None)
        assert out.shape[0] == 8

    def test_feature_names_propagate(self, deseq2_result):
        pytest.importorskip("sklearn")
        meta = pd.DataFrame(index=[f"S{i}" for i in range(5)])
        tr = _MockDETransformer(deseq2_result, metadata=meta)
        tr.fit()
        names = tr.get_feature_names_out()
        assert len(names) == 30
        assert all(isinstance(n, str) for n in names)
