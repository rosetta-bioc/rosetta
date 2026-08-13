"""Mock-based tests for wrappers/* and sklearn_compat — no live R session."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from unittest.mock import MagicMock, patch, PropertyMock


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

def _counts(n_genes=20, n_samples=6, seed=0):
    rng = np.random.default_rng(seed)
    df = pd.DataFrame(
        rng.integers(0, 500, (n_genes, n_samples)),
        index=[f"gene{i}" for i in range(n_genes)],
        columns=[f"s{i}" for i in range(n_samples)],
    )
    return df


def _meta(n_samples=6):
    return pd.DataFrame(
        {"condition": ["A"] * 3 + ["B"] * 3},
        index=[f"s{i}" for i in range(n_samples)],
    )


# ---------------------------------------------------------------------------
# _errors.py
# ---------------------------------------------------------------------------

def test_r_data_error_is_exception():
    from rosetta._errors import RDataError
    with pytest.raises(RDataError, match="bad data"):
        raise RDataError("bad data")


def test_r_formula_error_is_exception():
    from rosetta._errors import RFormulaError
    with pytest.raises(RFormulaError):
        raise RFormulaError("bad formula")


# ---------------------------------------------------------------------------
# _deps.py
# ---------------------------------------------------------------------------

def test_is_installed_builtin():
    from rosetta._deps import is_installed
    # os is always available — uses the R-package check path only if rpy2 exists
    # but is_installed should not raise for any string
    result = is_installed("nonexistent_r_pkg_xyz")
    assert isinstance(result, bool)


def test_ensure_installed_missing_raises(monkeypatch):
    from rosetta import _deps
    from rosetta._errors import RPackageMissing
    monkeypatch.setattr(_deps, "is_installed", lambda _: False)
    with pytest.raises(RPackageMissing, match="nonexistent_pkg"):
        _deps.ensure_installed("nonexistent_pkg")


def test_ensure_installed_present_does_not_raise(monkeypatch):
    from rosetta import _deps
    monkeypatch.setattr(_deps, "is_installed", lambda _: True)
    _deps.ensure_installed("any_pkg")  # should not raise


# ---------------------------------------------------------------------------
# wrappers/normalize.py — input validation (no R needed)
# ---------------------------------------------------------------------------

def test_vst_empty_counts_raises():
    from rosetta._errors import RDataError
    import rosetta.wrappers.normalize as norm
    with patch.object(norm, "importr", side_effect=ImportError("no rpy2")):
        with pytest.raises(RDataError, match="empty"):
            norm.vst(pd.DataFrame())


def test_vst_negative_counts_raises():
    from rosetta._errors import RDataError
    import rosetta.wrappers.normalize as norm
    df = _counts()
    df.iloc[0, 0] = -1
    with patch.object(norm, "importr", side_effect=ImportError("no rpy2")):
        with pytest.raises(RDataError, match="negative"):
            norm.vst(df, _meta())


def test_rlog_empty_raises():
    from rosetta._errors import RDataError
    import rosetta.wrappers.normalize as norm
    with patch.object(norm, "importr", side_effect=ImportError("no rpy2")):
        with pytest.raises(RDataError, match="empty"):
            norm.rlog(pd.DataFrame())


def test_tmm_empty_raises():
    from rosetta._errors import RDataError
    import rosetta.wrappers.normalize as norm
    with patch.object(norm, "importr", side_effect=ImportError("no rpy2")):
        with pytest.raises(RDataError, match="empty"):
            norm.tmm_normalize(pd.DataFrame())


def test_tmm_negative_raises():
    from rosetta._errors import RDataError
    import rosetta.wrappers.normalize as norm
    df = _counts()
    df.iloc[0, 0] = -5
    with patch.object(norm, "importr", side_effect=ImportError("no rpy2")):
        with pytest.raises(RDataError, match="negative"):
            norm.tmm_normalize(df)


# ---------------------------------------------------------------------------
# wrappers/deseq2.py — input validation
# ---------------------------------------------------------------------------

def test_deseq2_negative_counts_raises():
    from rosetta._errors import RDataError
    import rosetta.wrappers.deseq2 as mod
    counts = _counts()
    counts.iloc[0, 0] = -1

    mock_importr = MagicMock()
    mock_ensure = MagicMock()
    with patch.object(mod, "importr", mock_importr), \
         patch("rosetta.wrappers.deseq2.ensure_installed", mock_ensure):
        with pytest.raises(RDataError, match="negative"):
            mod.DESeq2(counts, _meta(), "~ condition")


def test_deseq2_mismatched_samples_raises():
    from rosetta._errors import RDataError
    import rosetta.wrappers.deseq2 as mod
    counts = _counts()
    meta = pd.DataFrame({"condition": ["A", "B"]}, index=["other1", "other2"])

    mock_importr = MagicMock()
    with patch.object(mod, "importr", mock_importr), \
         patch("rosetta.wrappers.deseq2.ensure_installed", MagicMock()):
        with pytest.raises(RDataError, match="columns must match"):
            mod.DESeq2(counts, meta, "~ condition")


# ---------------------------------------------------------------------------
# wrappers/edger.py — input validation
# ---------------------------------------------------------------------------

def test_edger_negative_counts_raises():
    from rosetta._errors import RDataError
    import rosetta.wrappers.edger as mod
    counts = _counts()
    counts.iloc[0, 0] = -1

    with patch.object(mod, "importr", MagicMock()), \
         patch("rosetta.wrappers.edger.ensure_installed", MagicMock()):
        with pytest.raises(RDataError, match="negative"):
            mod.EdgeR(counts, _meta())


def test_edger_mismatched_samples_raises():
    from rosetta._errors import RDataError
    import rosetta.wrappers.edger as mod
    meta = pd.DataFrame({"condition": ["A"]}, index=["other"])

    with patch.object(mod, "importr", MagicMock()), \
         patch("rosetta.wrappers.edger.ensure_installed", MagicMock()):
        with pytest.raises(RDataError, match="columns must match"):
            mod.EdgeR(_counts(), meta)


# ---------------------------------------------------------------------------
# wrappers/clusterprofiler.py — input validation
# ---------------------------------------------------------------------------

def test_clusterprofiler_empty_gene_list_raises():
    from rosetta._errors import RDataError
    import rosetta.wrappers.clusterprofiler as mod

    mock_pkg = MagicMock()
    with patch.object(mod, "importr", return_value=mock_pkg), \
         patch("rosetta.wrappers.clusterprofiler.ensure_installed", MagicMock()):
        cp = mod.ClusterProfiler.__new__(mod.ClusterProfiler)
        cp.cp_pkg = mock_pkg
        cp.obj = None

        with pytest.raises(RDataError, match="empty"):
            cp._run_enrich("enrichGO", [])


def test_clusterprofiler_none_gene_list_raises():
    from rosetta._errors import RDataError
    import rosetta.wrappers.clusterprofiler as mod

    mock_pkg = MagicMock()
    cp = mod.ClusterProfiler.__new__(mod.ClusterProfiler)
    cp.cp_pkg = mock_pkg
    cp.obj = None

    with pytest.raises(RDataError, match="empty"):
        cp._run_enrich("enrichGO", None)


def test_clusterprofiler_param_mapping():
    """Python snake_case kwargs get remapped to R CamelCase."""
    import rosetta.wrappers.clusterprofiler as mod

    mock_pkg = MagicMock()
    mock_pkg.enrichGO.return_value = MagicMock()

    cp = mod.ClusterProfiler.__new__(mod.ClusterProfiler)
    cp.cp_pkg = mock_pkg
    cp.obj = None

    with patch("rosetta.wrappers.clusterprofiler.localconverter"), \
         patch("rosetta.wrappers.clusterprofiler.to_pandas", return_value=pd.DataFrame()), \
         patch("rosetta.wrappers.clusterprofiler.to_r_df", return_value=MagicMock()):
        cp._run_enrich("enrichGO", ["BRCA1", "TP53"], pvalue_cutoff=0.05, min_gs_size=10)

    call_kwargs = mock_pkg.enrichGO.call_args[1]
    assert "pvalueCutoff" in call_kwargs
    assert "minGSSize" in call_kwargs
    assert "pvalue_cutoff" not in call_kwargs


# ---------------------------------------------------------------------------
# sklearn_compat.py — no R needed for structural checks
# ---------------------------------------------------------------------------

def test_sklearn_check_raises_without_sklearn(monkeypatch):
    import rosetta.sklearn_compat as sc
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "sklearn":
            raise ImportError("no sklearn")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="scikit-learn"):
        sc._check_sklearn()


def test_base_transformer_stores_params():
    import rosetta.sklearn_compat as sc

    counts = _counts()
    meta = _meta()

    # Patch _check_sklearn so we don't need sklearn installed
    with patch.object(sc, "_check_sklearn"):
        t = sc._BaseRosettaTransformer(
            counts, meta, design="~ condition", alpha=0.01, lfc_threshold=1.0,
            select="all", value="stat"
        )

    assert t.alpha == 0.01
    assert t.lfc_threshold == 1.0
    assert t.select == "all"
    assert t.value == "stat"
    assert t.counts is counts
    assert t.metadata is meta


def test_base_transformer_select_significant():
    """_select_genes filters to significant rows via result_ attribute."""
    import rosetta.sklearn_compat as sc

    df = pd.DataFrame({
        "log2FoldChange": [2.0, -1.5, 0.1, 3.0],
        "padj": [0.01, 0.03, 0.5, 0.001],
    }, index=[f"g{i}" for i in range(4)])

    with patch.object(sc, "_check_sklearn"):
        t = sc._BaseRosettaTransformer(_counts(), _meta(), alpha=0.05, lfc_threshold=0.0)
    t.result_ = df  # inject fitted result

    genes = t._select_genes()
    assert set(genes) == {"g0", "g1", "g3"}


def test_base_transformer_select_lfc_threshold():
    """lfc_threshold filters out small-effect genes."""
    import rosetta.sklearn_compat as sc

    df = pd.DataFrame({
        "log2FoldChange": [2.0, 0.3, -1.5],
        "padj": [0.01, 0.02, 0.03],
    }, index=["g0", "g1", "g2"])

    with patch.object(sc, "_check_sklearn"):
        t = sc._BaseRosettaTransformer(_counts(), _meta(), alpha=0.05, lfc_threshold=1.0)
    t.result_ = df

    genes = t._select_genes()
    assert "g1" not in genes  # |0.3| < 1.0
    assert "g0" in genes
    assert "g2" in genes


def test_base_transformer_select_all():
    """select='all' returns every gene."""
    import rosetta.sklearn_compat as sc

    df = pd.DataFrame({
        "log2FoldChange": [2.0, 0.1],
        "padj": [0.01, 0.9],
    }, index=["g0", "g1"])

    with patch.object(sc, "_check_sklearn"):
        t = sc._BaseRosettaTransformer(_counts(), _meta(), alpha=0.05, select="all")
    t.result_ = df

    genes = t._select_genes()
    assert set(genes) == {"g0", "g1"}


def test_base_transformer_no_genes_raises():
    """All genes filtered out → ValueError with helpful message."""
    import rosetta.sklearn_compat as sc

    df = pd.DataFrame({
        "log2FoldChange": [0.1, 0.2],
        "padj": [0.9, 0.8],
    }, index=["g0", "g1"])

    with patch.object(sc, "_check_sklearn"):
        t = sc._BaseRosettaTransformer(_counts(), _meta(), alpha=0.05, select="significant")
    t.result_ = df

    with pytest.raises(ValueError, match="No genes passed filters"):
        t._select_genes()


def test_base_transformer_padj_col_detection():
    """_padj_col detects FDR and adj.P.Val column names."""
    import rosetta.sklearn_compat as sc

    with patch.object(sc, "_check_sklearn"):
        t = sc._BaseRosettaTransformer(_counts(), _meta())

    for col in ("padj", "FDR", "adj.P.Val", "p.adjust"):
        t.result_ = pd.DataFrame({col: [0.01], "logFC": [1.0]})
        assert t._padj_col() == col
