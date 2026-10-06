"""Tests for rosetta.pipelines — high-level one-call workflows.

Differential-expression tests patch existing wrapper classes so stale imports
or calls to removed methods cannot be hidden by fake wrapper modules.
"""

import sys
import types
from unittest.mock import MagicMock, call, patch, sentinel

import pandas as pd
import pytest

from rosetta import _bridge
from rosetta.pipelines import compare, diff_expr, enrichment
from rosetta.results import RosettaDataFrame

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fake_df(sig_col="padj"):
    df = pd.DataFrame(
        {
            "log2FoldChange": [1.5, -0.3, 2.1, 0.1],
            sig_col: [0.01, 0.8, 0.001, 0.5],
            "baseMean": [100, 50, 200, 30],
        },
        index=["GeneA", "GeneB", "GeneC", "GeneD"],
    )
    return RosettaDataFrame(df)


@pytest.fixture
def counts():
    return pd.DataFrame(
        {"S1": [10, 20, 5], "S2": [12, 18, 7], "S3": [9, 22, 4], "S4": [11, 19, 6]},
        index=["G1", "G2", "G3"],
    )


@pytest.fixture
def metadata():
    return pd.DataFrame(
        {"condition": ["A", "A", "B", "B"]}, index=["S1", "S2", "S3", "S4"]
    )


# ---------------------------------------------------------------------------
# Wrapper fixtures — patch classes on their real modules
# ---------------------------------------------------------------------------

@pytest.fixture
def deseq2_class():
    with patch("rosetta.wrappers.deseq2.DESeq2", autospec=True) as wrapper:
        wrapper.return_value.get_results.return_value = _fake_df()
        wrapper.return_value.lfc_shrink.return_value = _fake_df()
        yield wrapper


@pytest.fixture
def edger_class():
    with patch("rosetta.wrappers.edger.EdgeR", autospec=True) as wrapper:
        wrapper.return_value.get_results.return_value = _fake_df("FDR")
        yield wrapper


@pytest.fixture
def limma_class():
    with patch("rosetta.wrappers.limma.Limma", autospec=True) as wrapper:
        wrapper.return_value.get_results.return_value = _fake_df("adj.P.Val")
        yield wrapper


def _cp_mod():
    mod = types.ModuleType("rosetta.wrappers.clusterprofiler")
    mod.enrich_go = MagicMock(return_value=RosettaDataFrame(pd.DataFrame({"ID": ["GO:1"]})))
    mod.enrich_kegg = MagicMock(return_value=RosettaDataFrame(pd.DataFrame({"ID": ["hsa1"]})))
    mod.enrich_pathway = MagicMock(return_value=RosettaDataFrame(pd.DataFrame({"ID": ["R-1"]})))
    return mod


# ---------------------------------------------------------------------------
# diff_expr — invalid method
# ---------------------------------------------------------------------------

def test_diff_expr_invalid_method_raises(counts, metadata):
    with pytest.raises(ValueError, match="Unknown method"):
        diff_expr(counts, metadata, method="magic")


# ---------------------------------------------------------------------------
# diff_expr — DESeq2
# ---------------------------------------------------------------------------

def test_diff_expr_deseq2_calls_run_and_get(counts, metadata, deseq2_class):
    model = deseq2_class.return_value

    result = diff_expr(counts, metadata)

    deseq2_class.assert_called_once_with(counts, metadata, "~ condition")
    assert model.method_calls == [
        call.run_deseq(),
        call.get_results(alpha=0.05, lfcThreshold=0.0),
    ]
    assert result is model.get_results.return_value
    assert result._rosetta_method == "deseq2"


def test_diff_expr_deseq2_passes_alpha_and_lfc(counts, metadata, deseq2_class):
    diff_expr(
        counts, metadata, design="~ batch + condition", method="deseq2",
        alpha=0.1, lfc_threshold=1.0,
    )

    deseq2_class.assert_called_once_with(counts, metadata, "~ batch + condition")
    deseq2_class.return_value.get_results.assert_called_once_with(
        alpha=0.1, lfcThreshold=1.0,
    )


@pytest.mark.parametrize("contrast", [None, []])
def test_diff_expr_deseq2_omits_absent_contrast(counts, metadata, deseq2_class, contrast):
    with patch.object(_bridge, "ro") as ro:
        diff_expr(counts, metadata, contrast=contrast)

    ro.StrVector.assert_not_called()
    deseq2_class.return_value.get_results.assert_called_once_with(
        alpha=0.05, lfcThreshold=0.0,
    )


def test_diff_expr_deseq2_converts_contrast(counts, metadata, deseq2_class):
    contrast = ["condition", "B", "A"]
    with patch.object(_bridge, "ro") as ro:
        ro.StrVector.return_value = sentinel.r_contrast
        diff_expr(counts, metadata, contrast=contrast)

    ro.StrVector.assert_called_once_with(contrast)
    deseq2_class.return_value.get_results.assert_called_once_with(
        alpha=0.05, lfcThreshold=0.0, contrast=sentinel.r_contrast,
    )


@pytest.mark.parametrize("shrinkage", ["apeglm", "ashr", "normal"])
def test_diff_expr_deseq2_with_shrinkage(counts, metadata, deseq2_class, shrinkage):
    model = deseq2_class.return_value
    model.r_obj = sentinel.unfitted_dds
    model.deseq_pkg = MagicMock()
    model.deseq_pkg.resultsNames.return_value = ["Intercept", "condition_B_vs_A"]

    def fit():
        model.r_obj = sentinel.fitted_dds

    model.run_deseq.side_effect = fit
    with patch.object(_bridge, "localconverter") as localconverter:
        result = diff_expr(counts, metadata, method="deseq2", shrinkage=shrinkage)

    localconverter.assert_called_once_with(_bridge._converter)
    assert model.method_calls == [
        call.run_deseq(),
        call.deseq_pkg.resultsNames(sentinel.fitted_dds),
        call.lfc_shrink(coef="condition_B_vs_A", type=shrinkage),
    ]
    model.get_results.assert_not_called()
    assert result is model.lfc_shrink.return_value
    assert result._rosetta_method == "deseq2"


def test_diff_expr_deseq2_shrinkage_without_coefficient(counts, metadata, deseq2_class):
    model = deseq2_class.return_value
    model.r_obj = sentinel.fitted_dds
    model.deseq_pkg = MagicMock()
    model.deseq_pkg.resultsNames.return_value = []

    with patch.object(_bridge, "localconverter"):
        result = diff_expr(
            counts, metadata, shrinkage="normal", alpha=0.1, lfc_threshold=1.0,
        )

    model.lfc_shrink.assert_not_called()
    assert model.method_calls == [
        call.run_deseq(),
        call.deseq_pkg.resultsNames(sentinel.fitted_dds),
        call.get_results(alpha=0.1, lfcThreshold=1.0),
    ]
    assert result is model.get_results.return_value
    assert result._rosetta_method == "deseq2"


# ---------------------------------------------------------------------------
# diff_expr — edgeR and limma
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("lfc_threshold", [0.0, 1.0])
def test_diff_expr_edger(counts, metadata, edger_class, lfc_threshold):
    model = edger_class.return_value
    model.run_test.return_value = sentinel.test_result

    result = diff_expr(
        counts, metadata, design="~ batch + condition", method="edger",
        lfc_threshold=lfc_threshold,
    )

    edger_class.assert_called_once_with(counts, metadata, "~ batch + condition")
    assert model.method_calls == [
        call.run_test(lfc=lfc_threshold),
        call.get_results(sentinel.test_result),
    ]
    assert result is model.get_results.return_value
    assert result._rosetta_method == "edger"


def test_diff_expr_limma(counts, metadata, limma_class):
    model = limma_class.return_value

    result = diff_expr(counts, metadata, design="~ batch + condition", method="limma")

    limma_class.assert_called_once_with(counts, metadata, "~ batch + condition")
    assert model.method_calls == [call.run_ebayes(), call.get_results()]
    assert result is model.get_results.return_value
    assert result._rosetta_method == "limma"


@pytest.mark.parametrize(
    ("method", "wrapper_path", "analysis_method"),
    [
        ("deseq2", "rosetta.wrappers.deseq2.DESeq2", "run_deseq"),
        ("edger", "rosetta.wrappers.edger.EdgeR", "run_test"),
        ("limma", "rosetta.wrappers.limma.Limma", "run_ebayes"),
    ],
)
def test_diff_expr_propagates_analysis_error(counts, metadata, method,
                                            wrapper_path, analysis_method):
    error = RuntimeError("Analysis failed")
    with patch(wrapper_path, autospec=True) as wrapper:
        model = wrapper.return_value
        getattr(model, analysis_method).side_effect = error

        with pytest.raises(RuntimeError) as exc:
            diff_expr(counts, metadata, method=method)

    assert exc.value is error
    model.get_results.assert_not_called()


# ---------------------------------------------------------------------------
# enrichment — invalid method
# ---------------------------------------------------------------------------

def test_enrichment_invalid_method_raises():
    with pytest.raises(ValueError, match="Unknown method"):
        enrichment(["1234"], method="pathway_x")


# ---------------------------------------------------------------------------
# enrichment — delegates to correct wrapper
# ---------------------------------------------------------------------------

def test_enrichment_go():
    mod = _cp_mod()
    with patch.dict(sys.modules, {"rosetta.wrappers.clusterprofiler": mod}):
        result = enrichment(["1234", "5678"], method="go")
    mod.enrich_go.assert_called_once()
    assert isinstance(result, pd.DataFrame)


def test_enrichment_kegg():
    mod = _cp_mod()
    with patch.dict(sys.modules, {"rosetta.wrappers.clusterprofiler": mod}):
        result = enrichment(["1234"], method="kegg", organism="hsa")
    mod.enrich_kegg.assert_called_once()
    assert isinstance(result, pd.DataFrame)


def test_enrichment_reactome():
    mod = _cp_mod()
    with patch.dict(sys.modules, {"rosetta.wrappers.clusterprofiler": mod}):
        result = enrichment(["1234"], method="reactome")
    mod.enrich_pathway.assert_called_once()
    assert isinstance(result, pd.DataFrame)


# ---------------------------------------------------------------------------
# compare — uses diff_expr internally, so we patch there
# ---------------------------------------------------------------------------

def test_compare_consensus(counts, metadata):
    with patch("rosetta.pipelines.diff_expr") as mock_de:
        mock_de.side_effect = [_fake_df("padj"), _fake_df("FDR"), _fake_df("adj.P.Val")]
        result = compare(counts, metadata, methods=["deseq2", "edger", "limma"])
    assert "n_methods" in result.columns
    assert mock_de.call_count == 3
    assert result["n_methods"].max() <= 3


def test_compare_partial_failure(counts, metadata):
    with patch("rosetta.pipelines.diff_expr") as mock_de:
        mock_de.side_effect = [
            _fake_df("padj"),
            RuntimeError("edgeR exploded"),
            _fake_df("adj.P.Val"),
        ]
        result = compare(counts, metadata, methods=["deseq2", "edger", "limma"])
    assert "n_methods" in result.columns
    assert result["n_methods"].max() <= 2


def test_compare_all_fail_raises(counts, metadata):
    with patch("rosetta.pipelines.diff_expr", side_effect=RuntimeError("broken")):
        with pytest.raises(RuntimeError, match="All methods failed"):
            compare(counts, metadata, methods=["deseq2", "edger", "limma"])


def test_compare_default_methods(counts, metadata):
    with patch("rosetta.pipelines.diff_expr", return_value=_fake_df()) as mock_de:
        compare(counts, metadata)
    assert mock_de.call_count == 3


def test_compare_returns_rosetta_dataframe(counts, metadata):
    with patch("rosetta.pipelines.diff_expr", return_value=_fake_df()):
        result = compare(counts, metadata)
    assert isinstance(result, RosettaDataFrame)


def test_compare_n_methods_is_row_sum(counts, metadata):
    with patch("rosetta.pipelines.diff_expr") as mock_de:
        mock_de.side_effect = [_fake_df("padj"), _fake_df("FDR")]
        result = compare(counts, metadata, methods=["deseq2", "edger"])
    assert (result["n_methods"] == result[["deseq2", "edger"]].sum(axis=1)).all()


def test_compare_uses_real_diff_expr_dispatch(counts, metadata, deseq2_class,
                                            edger_class, limma_class):
    result = compare(counts, metadata)

    assert result.columns.tolist() == ["deseq2", "edger", "limma", "n_methods"]
    assert result["n_methods"].to_dict() == {
        "GeneA": 3, "GeneC": 3, "GeneB": 0, "GeneD": 0,
    }
    for wrapper in (deseq2_class, edger_class, limma_class):
        wrapper.assert_called_once_with(counts, metadata, "~ condition")
        wrapper.return_value.get_results.assert_called_once()
