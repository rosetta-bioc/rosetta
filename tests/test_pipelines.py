"""Tests for rosetta.pipelines — high-level one-call workflows.

Strategy: pipelines.py uses lazy imports inside each function body.
We inject fake modules into sys.modules so those names resolve correctly.
"""

import sys
import types
import pandas as pd
import pytest
from unittest.mock import patch, MagicMock

from rosetta.pipelines import diff_expr, enrichment, compare
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
# Module factories — simulate lazy imports inside pipelines.py
# ---------------------------------------------------------------------------

def _deseq2_mod(dds=None, get_rv=None, shrink_rv=None):
    mod = types.ModuleType("rosetta.wrappers.deseq2")
    mod.run_deseq2 = MagicMock(return_value=dds or MagicMock())
    mod.get_results = MagicMock(return_value=get_rv or _fake_df())
    mod.lfc_shrink = MagicMock(return_value=shrink_rv or _fake_df())
    return mod


def _edger_mod(rv=None):
    mod = types.ModuleType("rosetta.wrappers.edger")
    mod.edger = MagicMock(return_value=rv or _fake_df("FDR"))
    return mod


def _limma_mod(rv=None):
    mod = types.ModuleType("rosetta.wrappers.limma")
    mod.limma_voom = MagicMock(return_value=rv or _fake_df("adj.P.Val"))
    return mod


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

def test_diff_expr_deseq2_calls_run_and_get(counts, metadata):
    mod = _deseq2_mod()
    with patch.dict(sys.modules, {"rosetta.wrappers.deseq2": mod}):
        result = diff_expr(counts, metadata, method="deseq2")
    mod.run_deseq2.assert_called_once()
    mod.get_results.assert_called_once()
    assert isinstance(result, pd.DataFrame)


def test_diff_expr_deseq2_passes_alpha_and_lfc(counts, metadata):
    mod = _deseq2_mod()
    with patch.dict(sys.modules, {"rosetta.wrappers.deseq2": mod}):
        diff_expr(counts, metadata, method="deseq2", alpha=0.1, lfc_threshold=1.0)
    _, kwargs = mod.get_results.call_args
    assert kwargs.get("alpha") == 0.1
    assert kwargs.get("lfc_threshold") == 1.0


def test_diff_expr_deseq2_with_shrinkage(counts, metadata):
    mod = _deseq2_mod()
    coef_list = ["Intercept", "condition_B_vs_A"]

    pkg_mock = MagicMock()
    pkg_mock.resultsNames.return_value = coef_list

    rpy2_pkg = types.ModuleType("rpy2.robjects.packages")
    rpy2_pkg.importr = MagicMock(return_value=pkg_mock)

    rpy2_conv = types.ModuleType("rpy2.robjects.conversion")
    rpy2_conv.localconverter = MagicMock()

    bridge = types.ModuleType("rosetta._bridge")
    bridge._converter = MagicMock()

    with patch.dict(sys.modules, {
        "rosetta.wrappers.deseq2": mod,
        "rpy2.robjects.packages": rpy2_pkg,
        "rpy2.robjects.conversion": rpy2_conv,
        "rosetta._bridge": bridge,
    }):
        result = diff_expr(counts, metadata, method="deseq2", shrinkage="normal")

    mod.lfc_shrink.assert_called_once()
    assert isinstance(result, pd.DataFrame)


# ---------------------------------------------------------------------------
# diff_expr — edgeR and limma
# ---------------------------------------------------------------------------

def test_diff_expr_edger(counts, metadata):
    mod = _edger_mod()
    with patch.dict(sys.modules, {"rosetta.wrappers.edger": mod}):
        result = diff_expr(counts, metadata, method="edger")
    mod.edger.assert_called_once()
    assert isinstance(result, pd.DataFrame)


def test_diff_expr_limma(counts, metadata):
    mod = _limma_mod()
    with patch.dict(sys.modules, {"rosetta.wrappers.limma": mod}):
        result = diff_expr(counts, metadata, method="limma")
    mod.limma_voom.assert_called_once()
    assert isinstance(result, pd.DataFrame)


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
