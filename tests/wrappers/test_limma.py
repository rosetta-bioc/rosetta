"""Tests for rosetta.wrappers.limma."""

import pandas as pd
import pytest
from unittest.mock import MagicMock, patch, PropertyMock

from rosetta._errors import RDataError, RFormulaError


# ---------------------------------------------------------------------------
# Fixtures (shared via conftest but redefined locally for mock tests)
# ---------------------------------------------------------------------------

@pytest.fixture
def counts():
    return pd.DataFrame(
        {"S1": [100, 200, 50], "S2": [120, 180, 60], "S3": [90, 210, 45]},
        index=["G1", "G2", "G3"],
    )


@pytest.fixture
def metadata():
    return pd.DataFrame(
        {"condition": ["ctrl", "ctrl", "treated"]}, index=["S1", "S2", "S3"]
    )


# ---------------------------------------------------------------------------
# Helper: build a fully-mocked Limma object without touching R
# ---------------------------------------------------------------------------

def _make_limma(counts, metadata, design="~ condition"):
    """Construct a Limma instance with all R calls mocked."""
    limma_pkg = MagicMock(name="limma")
    edger_pkg = MagicMock(name="edgeR")
    stats_pkg = MagicMock(name="stats")

    fit_obj = MagicMock(name="fit")
    edger_pkg.DGEList.return_value = MagicMock()
    edger_pkg.calcNormFactors.return_value = MagicMock()
    edger_pkg.voomLmFit.return_value = fit_obj

    with patch("rosetta.wrappers.limma.ensure_installed"), \
         patch("rosetta.wrappers.limma.importr", side_effect=[limma_pkg, edger_pkg, stats_pkg]), \
         patch("rosetta.wrappers.limma.ro") as mock_ro, \
         patch("rosetta.wrappers.limma.localconverter") as mock_lc, \
         patch("rosetta.wrappers.limma.to_r_matrix", return_value=MagicMock()), \
         patch("rosetta.wrappers.limma.to_r_dataframe", return_value=MagicMock()):
        mock_lc.return_value.__enter__ = MagicMock(return_value=None)
        mock_lc.return_value.__exit__ = MagicMock(return_value=False)
        mock_ro.Formula.return_value = MagicMock()

        from rosetta.wrappers.limma import Limma
        obj = Limma(counts, metadata, design=design)

    obj._limma_pkg = limma_pkg
    obj._edger_pkg = edger_pkg
    obj._fit_obj = fit_obj
    return obj


# ---------------------------------------------------------------------------
# Input validation — no R required
# ---------------------------------------------------------------------------

def test_negative_counts_raises(counts, metadata):
    from rosetta.wrappers.limma import Limma
    bad = counts.copy()
    bad.iloc[0, 0] = -1
    with pytest.raises(RDataError, match="negative"):
        with patch("rosetta.wrappers.limma.ensure_installed"), \
             patch("rosetta.wrappers.limma.importr", return_value=MagicMock()):
            Limma(bad, metadata, design="~ condition")


def test_mismatched_samples_raises(counts, metadata):
    from rosetta.wrappers.limma import Limma
    bad_meta = metadata.rename(index={"S1": "X999"})
    with pytest.raises(RDataError, match="columns must match"):
        with patch("rosetta.wrappers.limma.ensure_installed"), \
             patch("rosetta.wrappers.limma.importr", return_value=MagicMock()):
            Limma(counts, bad_meta, design="~ condition")


def test_bad_formula_raises(counts, metadata):
    from rosetta.wrappers.limma import Limma
    with pytest.raises(RFormulaError):
        with patch("rosetta.wrappers.limma.ensure_installed"), \
             patch("rosetta.wrappers.limma.importr", side_effect=[MagicMock(), MagicMock(), MagicMock()]), \
             patch("rosetta.wrappers.limma.ro") as mock_ro, \
             patch("rosetta.wrappers.limma.to_r_matrix", return_value=MagicMock()), \
             patch("rosetta.wrappers.limma.to_r_dataframe", return_value=MagicMock()):
            mock_ro.Formula.side_effect = ValueError("bad formula")
            Limma(counts, metadata, design="not a formula ~~~")


# ---------------------------------------------------------------------------
# run_fit — mock-based, covers lines 53–80
# ---------------------------------------------------------------------------

def test_run_fit_calls_voomlmfit(counts, metadata):
    lm = _make_limma(counts, metadata)
    # obj should be the voomLmFit return value
    assert lm.obj is not None


def test_run_fit_dge_chain(counts, metadata):
    """DGEList → calcNormFactors → voomLmFit are all called."""
    limma_pkg = MagicMock()
    edger_pkg = MagicMock()
    stats_pkg = MagicMock()
    fit_obj = MagicMock()
    edger_pkg.DGEList.return_value = MagicMock()
    edger_pkg.calcNormFactors.return_value = MagicMock()
    edger_pkg.voomLmFit.return_value = fit_obj

    with patch("rosetta.wrappers.limma.ensure_installed"), \
         patch("rosetta.wrappers.limma.importr", side_effect=[limma_pkg, edger_pkg, stats_pkg]), \
         patch("rosetta.wrappers.limma.ro") as mock_ro, \
         patch("rosetta.wrappers.limma.localconverter") as mock_lc, \
         patch("rosetta.wrappers.limma.to_r_matrix", return_value=MagicMock()), \
         patch("rosetta.wrappers.limma.to_r_dataframe", return_value=MagicMock()):
        mock_lc.return_value.__enter__ = MagicMock(return_value=None)
        mock_lc.return_value.__exit__ = MagicMock(return_value=False)
        mock_ro.Formula.return_value = MagicMock()
        from rosetta.wrappers.limma import Limma
        Limma(counts, metadata, design="~ condition")

    edger_pkg.DGEList.assert_called_once()
    edger_pkg.calcNormFactors.assert_called_once()
    edger_pkg.voomLmFit.assert_called_once()


# ---------------------------------------------------------------------------
# apply_contrasts — covers lines 62–68
# ---------------------------------------------------------------------------

def test_apply_contrasts_calls_contrasts_fit(counts, metadata):
    lm = _make_limma(counts, metadata)

    design_matrix = MagicMock()
    design_matrix.colnames = ["(Intercept)", "conditiontreated"]
    lm.obj = MagicMock()
    lm.obj.rx2.return_value = design_matrix
    lm._limma_pkg.contrasts_fit.return_value = MagicMock(name="contrasts_fit_result")

    with patch("rosetta.wrappers.limma.localconverter") as mock_lc, \
         patch("rosetta.wrappers.limma.build_contrast_matrix", return_value=MagicMock()):
        mock_lc.return_value.__enter__ = MagicMock(return_value=None)
        mock_lc.return_value.__exit__ = MagicMock(return_value=False)
        result = lm.apply_contrasts(["conditiontreated"])

    lm._limma_pkg.contrasts_fit.assert_called_once()
    assert result is lm  # fluent API returns self


# ---------------------------------------------------------------------------
# run_ebayes — covers line 71
# ---------------------------------------------------------------------------

def test_run_ebayes_delegates_to_call_r(counts, metadata):
    lm = _make_limma(counts, metadata)
    lm._call_r = MagicMock(return_value=lm)
    result = lm.run_ebayes(robust=True)
    lm._call_r.assert_called_once_with("eBayes", lm._PARAMS_EBAYES, robust=True)


# ---------------------------------------------------------------------------
# get_results — covers lines 74–80
# ---------------------------------------------------------------------------

def test_get_results_returns_dataframe(counts, metadata):
    lm = _make_limma(counts, metadata)
    expected = pd.DataFrame(
        {"logFC": [1.2, -0.5, 0.3], "adj.P.Val": [0.01, 0.05, 0.2]},
        index=["G1", "G2", "G3"],
    )

    with patch("rosetta.wrappers.limma.localconverter") as mock_lc, \
         patch("rosetta.wrappers.limma.to_pandas", return_value=expected), \
         patch("rosetta.wrappers.limma.r_nrow", return_value=3), \
         patch("rosetta.wrappers.limma.filter_kwargs", return_value={}):
        mock_lc.return_value.__enter__ = MagicMock(return_value=None)
        mock_lc.return_value.__exit__ = MagicMock(return_value=False)
        lm._limma_pkg.topTable.return_value = MagicMock()
        result = lm.get_results()

    assert isinstance(result, pd.DataFrame)
    assert "logFC" in result.columns


def test_get_results_injects_number_when_missing(counts, metadata):
    """number defaults to nrow(fit) when not supplied."""
    lm = _make_limma(counts, metadata)
    expected = pd.DataFrame({"logFC": [1.0]})

    with patch("rosetta.wrappers.limma.localconverter") as mock_lc, \
         patch("rosetta.wrappers.limma.to_pandas", return_value=expected), \
         patch("rosetta.wrappers.limma.r_nrow", return_value=42) as mock_nrow, \
         patch("rosetta.wrappers.limma.filter_kwargs", return_value={}):
        mock_lc.return_value.__enter__ = MagicMock(return_value=None)
        mock_lc.return_value.__exit__ = MagicMock(return_value=False)
        lm._limma_pkg.topTable.return_value = MagicMock()
        lm.get_results()

    _, call_kwargs = lm._limma_pkg.topTable.call_args
    assert call_kwargs.get("number") == 42


# ---------------------------------------------------------------------------
# Full R-based pipeline (skipped when limma not installed)
# ---------------------------------------------------------------------------

try:
    from rosetta._deps import is_installed as _is_installed
    _has_limma = _is_installed("limma")
except Exception:
    _has_limma = False

_limma_mark = pytest.mark.skipif(not _has_limma, reason="limma R package not installed")


@_limma_mark
def test_limma_pipeline(sample_counts, sample_metadata):
    from rosetta.wrappers.limma import Limma
    model = Limma(sample_counts, sample_metadata, design="~ condition")
    model.run_ebayes(verbose=False)
    result = model.get_results(coef=1)
    assert isinstance(result, pd.DataFrame)
    assert "logFC" in result.columns


@_limma_mark
def test_apply_contrasts_real(sample_counts, sample_metadata):
    from rosetta.wrappers.limma import Limma
    model = Limma(sample_counts, sample_metadata, design="~ condition")
    model.apply_contrasts("conditiontreated")
    model.run_ebayes()
    result = model.get_results()
    assert isinstance(result, pd.DataFrame)
