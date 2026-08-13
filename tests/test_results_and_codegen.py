"""Tests for results.py (RosettaDataFrame + _build_report) and codegen.py."""
import numpy as np
import pandas as pd
import pytest

import matplotlib
matplotlib.use("Agg")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _deseq2_df(n=30, n_sig=8):
    rng = np.random.default_rng(42)
    return pd.DataFrame({
        "baseMean": rng.uniform(10, 1000, n),
        "log2FoldChange": np.concatenate([rng.uniform(0.5, 3, n_sig // 2),
                                          rng.uniform(-3, -0.5, n_sig // 2),
                                          rng.uniform(-0.3, 0.3, n - n_sig)]),
        "padj": np.concatenate([rng.uniform(0.001, 0.04, n_sig),
                                  rng.uniform(0.1, 1.0, n - n_sig)]),
    }, index=[f"gene{i}" for i in range(n)])


def _edger_df(n=20, n_sig=5):
    rng = np.random.default_rng(1)
    return pd.DataFrame({
        "logFC": rng.uniform(-3, 3, n),
        "logCPM": rng.uniform(1, 10, n),
        "FDR": np.concatenate([rng.uniform(0.001, 0.04, n_sig),
                                rng.uniform(0.1, 1.0, n - n_sig)]),
    }, index=[f"gene{i}" for i in range(n)])


def _limma_df(n=20, n_sig=5):
    rng = np.random.default_rng(2)
    return pd.DataFrame({
        "logFC": rng.uniform(-3, 3, n),
        "AveExpr": rng.uniform(1, 10, n),
        "adj.P.Val": np.concatenate([rng.uniform(0.001, 0.04, n_sig),
                                      rng.uniform(0.1, 1.0, n - n_sig)]),
    }, index=[f"gene{i}" for i in range(n)])


def _enrichment_df():
    return pd.DataFrame({
        "Description": [f"GO term {i}" for i in range(10)],
        "GeneRatio": ["3/100"] * 10,
        "p.adjust": [0.001, 0.01, 0.03, 0.04, 0.06, 0.1, 0.2, 0.3, 0.5, 0.9],
    })


# ---------------------------------------------------------------------------
# RosettaDataFrame basics
# ---------------------------------------------------------------------------

def test_rosetta_df_is_dataframe():
    from rosetta.results import RosettaDataFrame
    rdf = RosettaDataFrame(_deseq2_df())
    assert isinstance(rdf, pd.DataFrame)


def test_rosetta_df_constructor_preserves_type():
    from rosetta.results import RosettaDataFrame
    rdf = RosettaDataFrame(_deseq2_df())
    sliced = rdf.iloc[:5]
    assert isinstance(sliced, RosettaDataFrame)


def test_rosetta_df_report_returns_string(capsys):
    from rosetta.results import RosettaDataFrame
    rdf = RosettaDataFrame(_deseq2_df())
    result = rdf.report()
    assert isinstance(result, str)
    assert "DESeq2" in result


def test_rosetta_df_report_with_method_prefix():
    from rosetta.results import RosettaDataFrame
    rdf = RosettaDataFrame(_deseq2_df())
    rdf._rosetta_method = "deseq2"
    text = rdf.report()
    assert "Method: deseq2" in text


# ---------------------------------------------------------------------------
# _build_report dispatch
# ---------------------------------------------------------------------------

def test_build_report_deseq2():
    from rosetta.results import _build_report
    text = _build_report(_deseq2_df())
    assert "DESeq2" in text
    assert "Upregulated" in text
    assert "Downregulated" in text


def test_build_report_edger():
    from rosetta.results import _build_report
    text = _build_report(_edger_df())
    assert "edgeR" in text


def test_build_report_limma():
    from rosetta.results import _build_report
    text = _build_report(_limma_df())
    assert "limma" in text


def test_build_report_enrichment():
    from rosetta.results import _build_report
    text = _build_report(_enrichment_df())
    assert "Enrichment" in text
    assert "GO term" in text  # top terms listed


def test_build_report_unknown_fallback():
    from rosetta.results import _build_report
    df = pd.DataFrame({"x": [1, 2], "y": [3, 4]})
    text = _build_report(df)
    assert "2 rows" in text


def test_build_report_deseq2_no_sig():
    """All padj > alpha → 0 significant, no LFC range line."""
    from rosetta.results import _build_report
    df = pd.DataFrame({
        "log2FoldChange": [0.1, -0.1],
        "padj": [0.5, 0.8],
    })
    text = _build_report(df, alpha=0.05)
    assert "DESeq2" in text
    assert "LFC range" not in text


def test_build_report_deseq2_with_nans():
    """NaN padj values are dropped gracefully."""
    from rosetta.results import _build_report
    df = pd.DataFrame({
        "log2FoldChange": [1.0, -1.0, 0.5],
        "padj": [0.01, float("nan"), 0.9],
    })
    text = _build_report(df, alpha=0.05)
    assert "DESeq2" in text


# ---------------------------------------------------------------------------
# RosettaDataFrame.volcano / .ma_plot (delegates to plots/)
# ---------------------------------------------------------------------------

def test_rosetta_df_volcano():
    import matplotlib.pyplot as plt
    from rosetta.results import RosettaDataFrame
    rdf = RosettaDataFrame(_deseq2_df())
    fig = rdf.volcano()
    assert fig is not None
    plt.close("all")


def test_rosetta_df_ma_plot():
    import matplotlib.pyplot as plt
    from rosetta.results import RosettaDataFrame
    rdf = RosettaDataFrame(_deseq2_df())
    fig = rdf.ma_plot()
    assert fig is not None
    plt.close("all")


# ---------------------------------------------------------------------------
# QuickResult
# ---------------------------------------------------------------------------

def test_quick_result_dict_access():
    from rosetta.quick_result import QuickResult
    qr = QuickResult({"a": 1, "b": 2}, method="test")
    assert qr["a"] == 1
    assert "b" in qr
    assert len(qr) == 2


def test_quick_result_iteration():
    from rosetta.quick_result import QuickResult
    qr = QuickResult({"x": 10, "y": 20}, method="test")
    assert set(qr) == {"x", "y"}
    assert set(qr.keys()) == {"x", "y"}
    assert set(qr.values()) == {10, 20}
    assert dict(qr.items()) == {"x": 10, "y": 20}


def test_quick_result_get_default():
    from rosetta.quick_result import QuickResult
    qr = QuickResult({}, method="test")
    assert qr.get("missing", 99) == 99


def test_quick_result_properties():
    from rosetta.quick_result import QuickResult
    qr = QuickResult({"k": "v"}, method="mymethod", metadata={"n": 5})
    assert qr.method == "mymethod"
    assert qr.metadata == {"n": 5}
    assert qr.data == {"k": "v"}


def test_quick_result_repr():
    from rosetta.quick_result import QuickResult
    qr = QuickResult({"a": 1}, method="seurat")
    assert "seurat" in repr(qr)


# --- report dispatchers ---

def test_quick_result_report_seurat(capsys):
    from rosetta.quick_result import QuickResult
    clusters = pd.Series([0, 0, 1, 1, 2], name="cluster")
    umap = pd.DataFrame({"UMAP_1": [1, 2, 3, 4, 5], "UMAP_2": [5, 4, 3, 2, 1]})
    qr = QuickResult(
        {"clusters": clusters, "umap": umap, "variable_features": list(range(200))},
        method="seurat",
    )
    text = qr.report()
    assert "Seurat" in text
    assert "cells" in text
    assert "Clusters" in text


def test_quick_result_report_phyloseq(capsys):
    from rosetta.quick_result import QuickResult
    diversity = pd.DataFrame({
        "Shannon": [2.1, 1.8, 2.5, 2.0],
        "Simpson": [0.8, 0.7, 0.9, 0.85],
    })
    qr = QuickResult({"diversity": diversity}, method="phyloseq")
    text = qr.report()
    assert "Phyloseq" in text
    assert "Shannon" in text


def test_quick_result_report_deseq2(capsys):
    from rosetta.quick_result import QuickResult
    qr = QuickResult({"results": _deseq2_df()}, method="deseq2")
    text = qr.report()
    assert "DESeq2" in text


def test_quick_result_report_edger(capsys):
    from rosetta.quick_result import QuickResult
    qr = QuickResult({"results": _edger_df()}, method="edger")
    text = qr.report()
    assert "edgeR" in text


def test_quick_result_report_deseq2_no_results_df():
    from rosetta.quick_result import QuickResult
    qr = QuickResult({}, method="deseq2")
    text = qr.report()
    assert "no results DataFrame" in text


def test_quick_result_report_edger_no_results_df():
    from rosetta.quick_result import QuickResult
    qr = QuickResult({}, method="edger")
    text = qr.report()
    assert "no results DataFrame" in text


def test_quick_result_report_generic(capsys):
    from rosetta.quick_result import QuickResult
    qr = QuickResult(
        {"df": pd.DataFrame({"a": [1]}), "lst": [1, 2, 3], "val": 42},
        method="unknown_tool",
    )
    text = qr.report()
    assert "QuickResult" in text
    assert "df" in text


# ---------------------------------------------------------------------------
# codegen.py
# ---------------------------------------------------------------------------

def test_codegen_enable_disable():
    from rosetta import codegen
    codegen.enable()
    assert codegen._is_enabled()
    codegen.disable()
    assert not codegen._is_enabled()


def test_codegen_emit_and_last():
    from rosetta import codegen
    codegen.enable()
    codegen._emit("dds <- DESeqDataSetFromMatrix(counts, meta, design=~condition)")
    log = codegen.last()
    assert "DESeqDataSetFromMatrix" in log


def test_codegen_clear():
    from rosetta import codegen
    codegen.enable()
    codegen._emit("some_r_code()")
    codegen.clear()
    assert codegen.last() == ""


def test_codegen_block():
    from rosetta import codegen
    codegen.enable()
    codegen.clear()
    codegen._block(["line1 <- foo()", "line2 <- bar()"])
    log = codegen.last()
    assert "line1" in log
    assert "line2" in log


def test_codegen_disabled_emit_still_records():
    """_emit records even when disabled (for later retrieval)."""
    from rosetta import codegen
    codegen.disable()
    codegen.clear()
    codegen._emit("silent_call()")
    assert "silent_call" in codegen.last()


def test_codegen_enable_clears_log():
    from rosetta import codegen
    codegen.disable()
    codegen.clear()
    codegen._emit("old_line()")
    codegen.enable()  # enable() calls clear()
    assert codegen.last() == ""
