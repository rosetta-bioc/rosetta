"""Tests for rosetta.wrappers.phyloseq."""

import sys
import types
import numpy as np
import pandas as pd
import pytest
from unittest.mock import MagicMock, patch

from rosetta._errors import RDataError


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def otu():
    return pd.DataFrame(
        {"S1": [10, 5, 0], "S2": [20, 8, 3], "S3": [15, 2, 1]},
        index=["OTU1", "OTU2", "OTU3"],
    )


@pytest.fixture
def sample_data():
    return pd.DataFrame({"condition": ["A", "A", "B"]}, index=["S1", "S2", "S3"])


@pytest.fixture
def tax_table():
    return pd.DataFrame(
        {"Phylum": ["Firmicutes", "Bacteroidetes", "Proteobacteria"]},
        index=["OTU1", "OTU2", "OTU3"],
    )


# ---------------------------------------------------------------------------
# Helper: build a fully-mocked Phyloseq object without touching R
# ---------------------------------------------------------------------------

def _make_phyloseq(otu, sample_data=None, tax_table=None):
    """Construct a Phyloseq with all R calls mocked out."""
    ps_pkg = MagicMock(name="ps_pkg")
    base_pkg = MagicMock(name="base_pkg")

    # phyloseq() constructor returns a fake R object
    r_obj = MagicMock(name="r_phyloseq_obj")
    ps_pkg.phyloseq.return_value = r_obj
    ps_pkg.otu_table.return_value = MagicMock()
    ps_pkg.sample_data.return_value = MagicMock()
    ps_pkg.tax_table.return_value = MagicMock()

    with patch("rosetta.wrappers.phyloseq.ensure_installed"), \
         patch("rosetta.wrappers.phyloseq.importr", side_effect=[ps_pkg, base_pkg]), \
         patch("rosetta.wrappers.phyloseq.localconverter") as mock_lc, \
         patch("rosetta.wrappers.phyloseq.to_r_matrix", return_value=MagicMock()), \
         patch("rosetta.wrappers.phyloseq.to_r_dataframe", return_value=MagicMock()):
        mock_lc.return_value.__enter__ = MagicMock(return_value=None)
        mock_lc.return_value.__exit__ = MagicMock(return_value=False)

        from rosetta.wrappers.phyloseq import Phyloseq
        obj = Phyloseq(otu, sample_data=sample_data, tax_table=tax_table)

    # Stash mocks for assertions
    obj._ps_pkg = ps_pkg
    obj._base_pkg = base_pkg
    return obj


# ---------------------------------------------------------------------------
# Input validation (no R required)
# ---------------------------------------------------------------------------

def test_empty_otu_raises(sample_data):
    from rosetta.wrappers.phyloseq import Phyloseq
    with pytest.raises(RDataError, match="empty"):
        with patch("rosetta.wrappers.phyloseq.ensure_installed"), \
             patch("rosetta.wrappers.phyloseq.importr", return_value=MagicMock()):
            Phyloseq(pd.DataFrame(), sample_data=sample_data)


def test_negative_otu_raises(otu, sample_data):
    from rosetta.wrappers.phyloseq import Phyloseq
    bad = otu.copy()
    bad.iloc[0, 0] = -1
    with pytest.raises(RDataError, match="negative"):
        with patch("rosetta.wrappers.phyloseq.ensure_installed"), \
             patch("rosetta.wrappers.phyloseq.importr", return_value=MagicMock()):
            Phyloseq(bad, sample_data=sample_data)


# ---------------------------------------------------------------------------
# Construction — various input combinations
# ---------------------------------------------------------------------------

def test_construct_otu_only(otu):
    ps = _make_phyloseq(otu)
    assert ps is not None


def test_construct_with_sample_data(otu, sample_data):
    ps = _make_phyloseq(otu, sample_data=sample_data)
    assert ps is not None


def test_construct_with_all_tables(otu, sample_data, tax_table):
    ps = _make_phyloseq(otu, sample_data=sample_data, tax_table=tax_table)
    assert ps is not None


def test_sample_data_component_added_when_provided(otu, sample_data):
    """sample_data component is passed to phyloseq() when provided."""
    ps_pkg = MagicMock()
    base_pkg = MagicMock()
    ps_pkg.phyloseq.return_value = MagicMock()
    ps_pkg.otu_table.return_value = MagicMock()
    sample_data_component = MagicMock()
    ps_pkg.sample_data.return_value = sample_data_component

    with patch("rosetta.wrappers.phyloseq.ensure_installed"), \
         patch("rosetta.wrappers.phyloseq.importr", side_effect=[ps_pkg, base_pkg]), \
         patch("rosetta.wrappers.phyloseq.localconverter") as mock_lc, \
         patch("rosetta.wrappers.phyloseq.to_r_matrix", return_value=MagicMock()), \
         patch("rosetta.wrappers.phyloseq.to_r_dataframe", return_value=MagicMock()):
        mock_lc.return_value.__enter__ = MagicMock(return_value=None)
        mock_lc.return_value.__exit__ = MagicMock(return_value=False)
        from rosetta.wrappers.phyloseq import Phyloseq
        Phyloseq(otu, sample_data=sample_data)

    ps_pkg.sample_data.assert_called_once()


def test_tax_table_component_added_when_provided(otu, tax_table):
    ps_pkg = MagicMock()
    base_pkg = MagicMock()
    ps_pkg.phyloseq.return_value = MagicMock()
    ps_pkg.otu_table.return_value = MagicMock()
    tax_component = MagicMock()
    ps_pkg.tax_table.return_value = tax_component

    with patch("rosetta.wrappers.phyloseq.ensure_installed"), \
         patch("rosetta.wrappers.phyloseq.importr", side_effect=[ps_pkg, base_pkg]), \
         patch("rosetta.wrappers.phyloseq.localconverter") as mock_lc, \
         patch("rosetta.wrappers.phyloseq.to_r_matrix", return_value=MagicMock()), \
         patch("rosetta.wrappers.phyloseq.to_r_dataframe", return_value=MagicMock()):
        mock_lc.return_value.__enter__ = MagicMock(return_value=None)
        mock_lc.return_value.__exit__ = MagicMock(return_value=False)
        from rosetta.wrappers.phyloseq import Phyloseq
        Phyloseq(otu, tax_table=tax_table)

    ps_pkg.tax_table.assert_called_once()


# ---------------------------------------------------------------------------
# estimate_richness
# ---------------------------------------------------------------------------

def test_estimate_richness_returns_dataframe(otu, sample_data):
    ps = _make_phyloseq(otu, sample_data=sample_data)
    expected = pd.DataFrame({"Shannon": [2.1, 1.8, 1.5]}, index=["S1", "S2", "S3"])

    with patch("rosetta.wrappers.phyloseq.localconverter") as mock_lc, \
         patch("rosetta.wrappers.phyloseq.to_pandas", return_value=expected):
        mock_lc.return_value.__enter__ = MagicMock(return_value=None)
        mock_lc.return_value.__exit__ = MagicMock(return_value=False)
        result = ps.estimate_richness()

    assert isinstance(result, pd.DataFrame)
    assert "Shannon" in result.columns


def test_estimate_richness_with_measures(otu, sample_data):
    ps = _make_phyloseq(otu, sample_data=sample_data)
    expected = pd.DataFrame({"Shannon": [2.1, 1.8, 1.5]})

    with patch("rosetta.wrappers.phyloseq.localconverter") as mock_lc, \
         patch("rosetta.wrappers.phyloseq.to_pandas", return_value=expected), \
         patch("rosetta.wrappers.phyloseq.ro") as mock_ro:
        mock_lc.return_value.__enter__ = MagicMock(return_value=None)
        mock_lc.return_value.__exit__ = MagicMock(return_value=False)
        mock_ro.StrVector.return_value = MagicMock()
        result = ps.estimate_richness(measures=["Shannon"])

    ps._ps_pkg.estimate_richness.assert_called_once()
    assert isinstance(result, pd.DataFrame)


# ---------------------------------------------------------------------------
# run_ordination
# ---------------------------------------------------------------------------

def test_run_ordination_returns_dataframe(otu, sample_data):
    ps = _make_phyloseq(otu, sample_data=sample_data)
    expected = pd.DataFrame({"Axis.1": [0.1, -0.2, 0.3], "Axis.2": [0.4, 0.5, -0.1]})

    ord_mock = MagicMock()
    ord_mock.rx2.return_value = MagicMock()
    ps._ps_pkg.ordinate.return_value = ord_mock
    ps._base_pkg.as_data_frame.return_value = MagicMock()

    with patch("rosetta.wrappers.phyloseq.localconverter") as mock_lc, \
         patch("rosetta.wrappers.phyloseq.to_pandas", return_value=expected), \
         patch("rosetta.wrappers.phyloseq.filter_kwargs", return_value={}):
        mock_lc.return_value.__enter__ = MagicMock(return_value=None)
        mock_lc.return_value.__exit__ = MagicMock(return_value=False)
        result = ps.run_ordination(method="PCoA", distance="bray")

    assert isinstance(result, pd.DataFrame)
    assert "Axis.1" in result.columns


def test_run_ordination_passes_kwargs(otu):
    ps = _make_phyloseq(otu)
    expected = pd.DataFrame({"Axis.1": [0.1]})

    ord_mock = MagicMock()
    ord_mock.rx2.return_value = MagicMock()
    ps._ps_pkg.ordinate.return_value = ord_mock
    ps._base_pkg.as_data_frame.return_value = MagicMock()

    filtered = {"method": "NMDS", "distance": "jaccard"}

    with patch("rosetta.wrappers.phyloseq.localconverter") as mock_lc, \
         patch("rosetta.wrappers.phyloseq.to_pandas", return_value=expected), \
         patch("rosetta.wrappers.phyloseq.filter_kwargs", return_value=filtered):
        mock_lc.return_value.__enter__ = MagicMock(return_value=None)
        mock_lc.return_value.__exit__ = MagicMock(return_value=False)
        ps.run_ordination(method="NMDS", distance="jaccard")

    ps._ps_pkg.ordinate.assert_called_once_with(ps.obj, **filtered)


# ---------------------------------------------------------------------------
# _phyloseq_available helper
# ---------------------------------------------------------------------------

def test_phyloseq_available_true():
    with patch("rosetta._deps.is_installed", return_value=True):
        from rosetta.wrappers.phyloseq import _phyloseq_available
        assert _phyloseq_available() is True


def test_phyloseq_available_false():
    with patch("rosetta._deps.is_installed", return_value=False):
        from rosetta.wrappers.phyloseq import _phyloseq_available
        assert _phyloseq_available() is False


def test_phyloseq_available_exception_returns_false():
    with patch("rosetta._deps.is_installed", side_effect=Exception("boom")):
        from rosetta.wrappers.phyloseq import _phyloseq_available
        assert _phyloseq_available() is False
