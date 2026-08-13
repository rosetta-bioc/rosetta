"""Targeted mock tests to close coverage gaps in phyloseq, limma, deseq2, decide, design."""
from __future__ import annotations

import pandas as pd
import numpy as np
import pytest
from unittest.mock import MagicMock, patch, call
from contextlib import ExitStack, contextmanager


# =========================================================================
# Helpers
# =========================================================================

def _otu(n_taxa=4, n_samples=3):
    rng = np.random.default_rng(0)
    return pd.DataFrame(
        rng.integers(1, 200, (n_taxa, n_samples)),
        index=[f"OTU{i}" for i in range(n_taxa)],
        columns=[f"S{i}" for i in range(n_samples)],
    )


def _sample_data(n=3):
    return pd.DataFrame({"group": ["A", "A", "B"]}, index=[f"S{i}" for i in range(n)])


def _tax(n=4):
    return pd.DataFrame(
        {"Phylum": [f"P{i}" for i in range(n)]},
        index=[f"OTU{i}" for i in range(n)],
    )


@contextmanager
def _stack(cms):
    with ExitStack() as stack:
        for cm in cms:
            stack.enter_context(cm)
        yield


# =========================================================================
# phyloseq — lines 35–85 (entire __init__ / _create_object / methods)
# =========================================================================

def _make_phyloseq(otu=None, sample_data=None, tax_table=None):
    """Build a Phyloseq with all R mocked out."""
    import rosetta.wrappers.phyloseq as mod
    otu = otu if otu is not None else _otu()

    ps_pkg  = MagicMock(name="phyloseq")
    base_pkg = MagicMock(name="base")
    ps_obj  = MagicMock(name="ps_obj")
    ps_pkg.phyloseq.return_value = ps_obj

    ctx = [
        patch("rosetta.wrappers.phyloseq.ensure_installed"),
        patch("rosetta.wrappers.phyloseq.importr", side_effect=[ps_pkg, base_pkg]),
        patch("rosetta.wrappers.phyloseq.localconverter"),
        patch("rosetta.wrappers.phyloseq._converter", MagicMock()),
        patch("rosetta.wrappers.phyloseq.to_r_matrix", return_value=MagicMock()),
        patch("rosetta.wrappers.phyloseq.to_r_dataframe", return_value=MagicMock()),
        patch("rosetta.wrappers.phyloseq.to_pandas", return_value=pd.DataFrame({"x": [1.0]})),
    ]
    with _stack(ctx):
        obj = mod.Phyloseq(otu, sample_data=sample_data, tax_table=tax_table)
    obj.ps_pkg  = ps_pkg
    obj.base_pkg = base_pkg
    obj._ps_obj = ps_obj
    return obj, ps_pkg, base_pkg


def test_phyloseq_empty_otu_raises():
    from rosetta.wrappers.phyloseq import Phyloseq
    from rosetta._errors import RDataError
    with pytest.raises(RDataError, match="empty"):
        Phyloseq(pd.DataFrame())


def test_phyloseq_negative_otu_raises():
    from rosetta.wrappers.phyloseq import Phyloseq
    from rosetta._errors import RDataError
    bad = _otu()
    bad.iloc[0, 0] = -1
    with pytest.raises(RDataError, match="negative"):
        with patch("rosetta.wrappers.phyloseq.ensure_installed"):
            Phyloseq(bad)


def test_phyloseq_init_calls_importr():
    import rosetta.wrappers.phyloseq as mod
    ps_pkg = MagicMock(name="phyloseq")
    base_pkg = MagicMock(name="base")
    ps_pkg.phyloseq.return_value = MagicMock()

    importr_mock = MagicMock(side_effect=[ps_pkg, base_pkg])
    ctx = [
        patch("rosetta.wrappers.phyloseq.ensure_installed"),
        patch("rosetta.wrappers.phyloseq.importr", importr_mock),
        patch("rosetta.wrappers.phyloseq.localconverter"),
        patch("rosetta.wrappers.phyloseq._converter", MagicMock()),
        patch("rosetta.wrappers.phyloseq.to_r_matrix", return_value=MagicMock()),
    ]
    with _stack(ctx):
        mod.Phyloseq(_otu())

    importr_mock.assert_any_call("phyloseq")
    importr_mock.assert_any_call("base")


def test_phyloseq_create_object_with_sample_and_tax():
    """_create_object appends sample_data and tax_table components."""
    import rosetta.wrappers.phyloseq as mod
    ps_pkg = MagicMock(name="phyloseq")
    base_pkg = MagicMock(name="base")
    ps_pkg.phyloseq.return_value = MagicMock()

    ctx = [
        patch("rosetta.wrappers.phyloseq.ensure_installed"),
        patch("rosetta.wrappers.phyloseq.importr", side_effect=[ps_pkg, base_pkg]),
        patch("rosetta.wrappers.phyloseq.localconverter"),
        patch("rosetta.wrappers.phyloseq._converter", MagicMock()),
        patch("rosetta.wrappers.phyloseq.to_r_matrix", return_value=MagicMock()),
        patch("rosetta.wrappers.phyloseq.to_r_dataframe", return_value=MagicMock()),
    ]
    with _stack(ctx):
        mod.Phyloseq(_otu(), sample_data=_sample_data(), tax_table=_tax())

    ps_pkg.sample_data.assert_called_once()
    ps_pkg.tax_table.assert_called_once()
    # phyloseq(*components) should be called with 3 args
    assert ps_pkg.phyloseq.call_count == 1


def test_phyloseq_create_object_otu_only():
    """Without optional args, only otu_table component is added."""
    import rosetta.wrappers.phyloseq as mod
    ps_pkg = MagicMock()
    base_pkg = MagicMock()
    ps_pkg.phyloseq.return_value = MagicMock()

    ctx = [
        patch("rosetta.wrappers.phyloseq.ensure_installed"),
        patch("rosetta.wrappers.phyloseq.importr", side_effect=[ps_pkg, base_pkg]),
        patch("rosetta.wrappers.phyloseq.localconverter"),
        patch("rosetta.wrappers.phyloseq._converter", MagicMock()),
        patch("rosetta.wrappers.phyloseq.to_r_matrix", return_value=MagicMock()),
    ]
    with _stack(ctx):
        mod.Phyloseq(_otu())

    ps_pkg.sample_data.assert_not_called()
    ps_pkg.tax_table.assert_not_called()


def test_phyloseq_estimate_richness_no_measures():
    obj, ps_pkg, _ = _make_phyloseq()
    result_df = pd.DataFrame({"Shannon": [2.1, 1.9, 2.4]})

    with patch("rosetta.wrappers.phyloseq.localconverter"), \
         patch("rosetta.wrappers.phyloseq.to_pandas", return_value=result_df):
        result = obj.estimate_richness()

    ps_pkg.estimate_richness.assert_called_once()
    assert list(result.columns) == ["Shannon"]


def test_phyloseq_estimate_richness_with_measures():
    obj, ps_pkg, _ = _make_phyloseq()
    result_df = pd.DataFrame({"Shannon": [2.1, 1.9, 2.4]})

    with patch("rosetta.wrappers.phyloseq.localconverter"), \
         patch("rosetta.wrappers.phyloseq.to_pandas", return_value=result_df), \
         patch("rosetta.wrappers.phyloseq.ro") as mock_ro:
        mock_ro.StrVector.return_value = MagicMock()
        result = obj.estimate_richness(measures=["Shannon"])

    mock_ro.StrVector.assert_called_once_with(["Shannon"])
    ps_pkg.estimate_richness.assert_called_once()


def test_phyloseq_run_ordination():
    obj, ps_pkg, base_pkg = _make_phyloseq()
    coord_df = pd.DataFrame({"Axis.1": [0.1, -0.2, 0.3], "Axis.2": [0.4, 0.5, -0.1]})
    ord_obj = MagicMock()
    ps_pkg.ordinate.return_value = ord_obj

    with patch("rosetta.wrappers.phyloseq.localconverter"), \
         patch("rosetta.wrappers.phyloseq.to_pandas", return_value=coord_df), \
         patch("rosetta.wrappers.phyloseq.filter_kwargs", return_value={"method": "NMDS"}):
        result = obj.run_ordination(method="NMDS", distance="bray")

    ps_pkg.ordinate.assert_called_once()
    ord_obj.rx2.assert_called_with("vectors")
    base_pkg.as_data_frame.assert_called_once()
    assert "Axis.1" in result.columns


def test_phyloseq_run_ordination_filters_kwargs():
    """Only PARAMS_ORDINATION keys are forwarded to R."""
    obj, ps_pkg, base_pkg = _make_phyloseq()
    coord_df = pd.DataFrame({"PC1": [0.1]})
    ord_obj = MagicMock()
    ps_pkg.ordinate.return_value = ord_obj

    with patch("rosetta.wrappers.phyloseq.localconverter"), \
         patch("rosetta.wrappers.phyloseq.to_pandas", return_value=coord_df):
        # 'invalid_param' should be stripped by filter_kwargs
        result = obj.run_ordination(method="PCoA", distance="unifrac", invalid_param="foo")

    call_kwargs = ps_pkg.ordinate.call_args
    assert "invalid_param" not in str(call_kwargs)


# =========================================================================
# limma — lines 53–80 (init internals, apply_contrasts, run_ebayes, get_results)
# =========================================================================

def _make_limma():
    import rosetta.wrappers.limma as mod
    counts = pd.DataFrame(
        {"S1": [100, 200], "S2": [150, 180], "S3": [90, 210]},
        index=["G1", "G2"],
    )
    meta = pd.DataFrame({"condition": ["ctrl", "ctrl", "treated"]}, index=["S1", "S2", "S3"])

    limma_pkg = MagicMock(name="limma")
    edger_pkg = MagicMock(name="edgeR")
    stats_pkg = MagicMock(name="stats")
    fit_obj   = MagicMock(name="fit")
    edger_pkg.voomLmFit.return_value = fit_obj

    ctx = [
        patch("rosetta.wrappers.limma.ensure_installed"),
        patch("rosetta.wrappers.limma.importr", side_effect=[limma_pkg, edger_pkg, stats_pkg]),
        patch("rosetta.wrappers.limma.ro", MagicMock()),
        patch("rosetta.wrappers.limma.localconverter"),
        patch("rosetta.wrappers.limma._converter", MagicMock()),
        patch("rosetta.wrappers.limma.to_r_matrix", return_value=MagicMock()),
        patch("rosetta.wrappers.limma.to_r_dataframe", return_value=MagicMock()),
    ]
    with _stack(ctx):
        obj = mod.Limma(counts, meta)
    obj.limma_pkg = limma_pkg
    obj.edger_pkg = edger_pkg
    obj.obj = fit_obj
    return obj, limma_pkg, edger_pkg, fit_obj


def test_limma_init_calls_voomlmfit():
    _, _, edger_pkg, _ = _make_limma()
    edger_pkg.voomLmFit.assert_called_once()


def test_limma_init_dge_chain():
    _, _, edger_pkg, _ = _make_limma()
    edger_pkg.DGEList.assert_called_once()
    edger_pkg.calcNormFactors.assert_called_once()


def test_limma_apply_contrasts_calls_contrasts_fit():
    obj, limma_pkg, _, fit_obj = _make_limma()
    design_matrix = MagicMock()
    design_matrix.colnames = ["(Intercept)", "conditiontreated"]
    fit_obj.rx2.return_value = design_matrix
    new_fit = MagicMock()
    limma_pkg.contrasts_fit.return_value = new_fit

    with patch("rosetta.wrappers.limma.localconverter"), \
         patch("rosetta.wrappers.limma.build_contrast_matrix", return_value=MagicMock()) as mock_bcm:
        result = obj.apply_contrasts(["conditiontreated - (Intercept)"])

    mock_bcm.assert_called_once()
    limma_pkg.contrasts_fit.assert_called_once()
    assert result is obj  # returns self


def test_limma_run_ebayes_delegates():
    obj, limma_pkg, _, fit_obj = _make_limma()
    ebayes_result = MagicMock()
    limma_pkg.eBayes.return_value = ebayes_result

    with patch.object(obj, "_call_r", return_value=ebayes_result) as mock_call:
        result = obj.run_ebayes(trend=True)

    mock_call.assert_called_once_with("eBayes", obj._PARAMS_EBAYES, trend=True)


def test_limma_get_results_calls_toptable():
    obj, limma_pkg, _, fit_obj = _make_limma()
    result_df = pd.DataFrame({"logFC": [1.2, -0.8], "adj.P.Val": [0.01, 0.04]})
    limma_pkg.topTable.return_value = MagicMock()

    with patch("rosetta.wrappers.limma.localconverter"), \
         patch("rosetta.wrappers.limma.to_pandas", return_value=result_df), \
         patch("rosetta.wrappers.limma.r_nrow", return_value=2):
        result = obj.get_results()

    limma_pkg.topTable.assert_called_once()
    assert "logFC" in result.columns


def test_limma_get_results_default_number():
    """number kwarg is auto-set to nrow when not provided."""
    obj, limma_pkg, _, _ = _make_limma()
    result_df = pd.DataFrame({"logFC": [1.0]})

    with patch("rosetta.wrappers.limma.localconverter"), \
         patch("rosetta.wrappers.limma.to_pandas", return_value=result_df), \
         patch("rosetta.wrappers.limma.r_nrow", return_value=42) as mock_nrow:
        obj.get_results()

    call_kwargs = limma_pkg.topTable.call_args[1]
    assert call_kwargs.get("number") == 42


def test_limma_get_results_respects_explicit_number():
    obj, limma_pkg, _, _ = _make_limma()
    result_df = pd.DataFrame({"logFC": [1.0]})

    with patch("rosetta.wrappers.limma.localconverter"), \
         patch("rosetta.wrappers.limma.to_pandas", return_value=result_df), \
         patch("rosetta.wrappers.limma.r_nrow", return_value=99):
        obj.get_results(number=10)

    call_kwargs = limma_pkg.topTable.call_args[1]
    assert call_kwargs.get("number") == 10


# =========================================================================
# deseq2 — lines 72,74,86,87,113,126,128,130,133,135,137–141
# =========================================================================

def _make_deseq2():
    import rosetta.wrappers.deseq2 as mod
    counts = pd.DataFrame(
        {"S1": [10, 20], "S2": [15, 18], "S3": [9, 21]},
        index=["G1", "G2"],
    )
    meta = pd.DataFrame({"condition": ["A", "A", "B"]}, index=["S1", "S2", "S3"])

    deseq_pkg = MagicMock(name="DESeq2")
    stats_pkg = MagicMock(name="stats")
    dds_obj   = MagicMock(name="dds")
    deseq_pkg.DESeqDataSetFromMatrix.return_value = dds_obj

    ctx = [
        patch("rosetta.wrappers.deseq2.ensure_installed"),
        patch("rosetta.wrappers.deseq2.importr", side_effect=[deseq_pkg, stats_pkg]),
        patch("rosetta.wrappers.deseq2.localconverter"),
        patch("rosetta.wrappers.deseq2._converter", MagicMock()),
        patch("rosetta.wrappers.deseq2.to_r_matrix", return_value=MagicMock()),
        patch("rosetta.wrappers.deseq2.to_r_df", return_value=MagicMock()),
        patch("rosetta.wrappers.deseq2.codegen", MagicMock()),
    ]
    with _stack(ctx):
        obj = mod.DESeq2(counts, meta, "~ condition")
    obj.deseq_pkg = deseq_pkg
    obj.obj = dds_obj
    return obj, deseq_pkg, dds_obj


def _deseq2_ctx(deseq_pkg=None, stats_pkg=None):
    deseq_pkg = deseq_pkg or MagicMock()
    stats_pkg = stats_pkg or MagicMock()
    return [
        patch("rosetta.wrappers.deseq2.ensure_installed"),
        patch("rosetta.wrappers.deseq2.importr", side_effect=[deseq_pkg, stats_pkg]),
        patch("rosetta.wrappers.deseq2.localconverter"),
        patch("rosetta.wrappers.deseq2._converter", MagicMock()),
        patch("rosetta.wrappers.deseq2.to_r_matrix", return_value=MagicMock()),
        patch("rosetta.wrappers.deseq2.to_r_df", return_value=MagicMock()),
        patch("rosetta.wrappers.deseq2.codegen", MagicMock()),
    ]


def test_deseq2_negative_counts_raises():
    from rosetta.wrappers.deseq2 import DESeq2
    from rosetta._errors import RDataError
    counts = pd.DataFrame({"S1": [-1, 20], "S2": [15, 18]}, index=["G1", "G2"])
    meta = pd.DataFrame({"condition": ["A", "B"]}, index=["S1", "S2"])
    with _stack(_deseq2_ctx()):
        with pytest.raises(RDataError, match="negative"):
            DESeq2(counts, meta, "~ condition")


def test_deseq2_mismatched_samples_raises():
    from rosetta.wrappers.deseq2 import DESeq2
    from rosetta._errors import RDataError
    counts = pd.DataFrame({"S1": [10, 20], "S2": [15, 18]}, index=["G1", "G2"])
    meta = pd.DataFrame({"condition": ["A", "B"]}, index=["X1", "X2"])
    with _stack(_deseq2_ctx()):
        with pytest.raises(RDataError, match="columns must match"):
            DESeq2(counts, meta, "~ condition")


def test_deseq2_bad_formula_raises():
    from rosetta.wrappers.deseq2 import DESeq2
    from rosetta._errors import RDataError, RFormulaError
    counts = pd.DataFrame({"S1": [10], "S2": [15]}, index=["G1"])
    meta = pd.DataFrame({"condition": ["A", "B"]}, index=["S1", "S2"])

    stats_pkg = MagicMock()
    stats_pkg.as_formula.side_effect = Exception("bad formula")
    with _stack(_deseq2_ctx(stats_pkg=stats_pkg)):
        with pytest.raises(RFormulaError):
            DESeq2(counts, meta, "~~ bad")


def test_deseq2_run_deseq_calls_deseq():
    obj, deseq_pkg, dds_obj = _make_deseq2()
    deseq_pkg.DESeq.return_value = MagicMock()

    with patch("rosetta.wrappers.deseq2.localconverter"), \
         patch("rosetta.wrappers.deseq2.codegen", MagicMock()), \
         patch.object(obj, "_call_r", return_value=obj) as mock_call:
        obj.run_deseq()

    mock_call.assert_called_once_with("DESeq", obj._PARAMS_DESEQ)


def test_deseq2_lfc_shrink_invalid_type_raises():
    obj, _, _ = _make_deseq2()
    with pytest.raises(ValueError, match="Invalid shrinkage type"):
        obj.lfc_shrink(type="badtype")


def test_deseq2_lfc_shrink_valid_types():
    obj, deseq_pkg, _ = _make_deseq2()
    result_df = pd.DataFrame({"log2FoldChange": [1.0]})
    deseq_pkg.lfcShrink.return_value = MagicMock()

    for t in ("apeglm", "ashr", "normal"):
        with patch("rosetta.wrappers.deseq2.localconverter"), \
             patch("rosetta.wrappers.deseq2.to_pandas", return_value=result_df), \
             patch("rosetta.wrappers.deseq2.to_r_df", return_value=MagicMock()):
            result = obj.lfc_shrink(type=t)
        assert "log2FoldChange" in result.columns


def test_deseq2_get_results_calls_results():
    obj, deseq_pkg, _ = _make_deseq2()
    result_df = pd.DataFrame({"padj": [0.01, 0.5]})
    deseq_pkg.results.return_value = MagicMock()

    with patch("rosetta.wrappers.deseq2.localconverter"), \
         patch("rosetta.wrappers.deseq2.to_pandas", return_value=result_df), \
         patch("rosetta.wrappers.deseq2.to_r_df", return_value=MagicMock()):
        result = obj.get_results()

    deseq_pkg.results.assert_called_once()
    assert "padj" in result.columns


def test_deseq2_legacy_run_deseq2():
    """Legacy bridge function constructs DESeq2 and calls run_deseq."""
    import rosetta.wrappers.deseq2 as mod
    counts = pd.DataFrame({"S1": [10], "S2": [15]}, index=["G1"])
    meta = pd.DataFrame({"condition": ["A", "B"]}, index=["S1", "S2"])

    mock_model = MagicMock()
    mock_model.run_deseq.return_value = mock_model
    mock_model.obj = MagicMock()

    with patch.object(mod, "DESeq2", return_value=mock_model):
        result = mod.run_deseq2(counts, meta, "~ condition")

    mock_model.run_deseq.assert_called_once()
    assert result is mock_model.obj


# =========================================================================
# stats/decide — lines 30,31 (RuntimeError on failure)
# =========================================================================

def test_decide_tests_runtime_error_on_failure():
    from rosetta.stats.decide import run_decide_tests
    limma_pkg = MagicMock()
    limma_pkg.decideTests.side_effect = Exception("R error")

    with patch("rosetta.stats.decide.importr", return_value=limma_pkg):
        with pytest.raises(RuntimeError, match="decideTests\\(\\) failed"):
            run_decide_tests(MagicMock())


def test_decide_tests_success():
    from rosetta.stats.decide import run_decide_tests
    limma_pkg = MagicMock()
    result_mat = MagicMock()
    limma_pkg.decideTests.return_value = result_mat
    fit_obj = MagicMock()

    with patch("rosetta.stats.decide.importr", return_value=limma_pkg):
        result = run_decide_tests(fit_obj, method="global", adj="fdr", p_value=0.01)

    limma_pkg.decideTests.assert_called_once_with(
        fit_obj,
        method="global",
        adjust_method="fdr",
        p_value=0.01,
    )
    assert result is result_mat


# =========================================================================
# stats/design — lines 19,27,28 (empty contrast → None, rpy2 unavailable, failure)
# =========================================================================

def test_build_contrast_matrix_empty_returns_none():
    from rosetta.stats.design import build_contrast_matrix
    assert build_contrast_matrix(["A", "B"], "") is None
    assert build_contrast_matrix(["A", "B"], None) is None


def test_build_contrast_matrix_no_rpy2_raises():
    import rosetta.stats.design as mod
    original = mod.importr
    mod.importr = None
    try:
        with pytest.raises(RuntimeError, match="rpy2 is required"):
            mod.build_contrast_matrix(["A", "B"], "A - B")
    finally:
        mod.importr = original


def test_build_contrast_matrix_r_failure_raises():
    from rosetta.stats.design import build_contrast_matrix
    from rosetta._errors import RosettaSecurityError
    limma_pkg = MagicMock()
    limma_pkg.makeContrasts.side_effect = Exception("R crash")

    with patch("rosetta.stats.design.importr", return_value=limma_pkg), \
         patch("rosetta.stats.design.ro") as mock_ro:
        mock_ro.StrVector.return_value = MagicMock()
        with pytest.raises(RosettaSecurityError, match="Failed to build contrast matrix"):
            build_contrast_matrix(["condA", "condB"], "condA - condB")
