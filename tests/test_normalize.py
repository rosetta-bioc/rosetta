"""Tests for normalization wrappers (VST, rlog, TMM)."""

import numpy as np
import pandas as pd
import pytest

from rosetta import vst, rlog, tmm_normalize
from rosetta._deps import is_installed
from rosetta._errors import RDataError

needs_deseq2 = pytest.mark.skipif(
    not is_installed("DESeq2"), reason="DESeq2 R package not installed"
)
needs_edger = pytest.mark.skipif(
    not is_installed("edgeR"), reason="edgeR R package not installed"
)


@pytest.fixture
def norm_counts():
    """Larger count matrix suitable for normalization (10 genes x 6 samples).

    DESeq2 VST/rlog need more genes and samples than minimal 3x4 to run
    without warnings or errors about insufficient replication.
    """
    np.random.seed(123)
    data = np.random.negative_binomial(5, 0.1, size=(100, 6))
    return pd.DataFrame(
        data,
        index=[f"Gene{i}" for i in range(100)],
        columns=["S1", "S2", "S3", "S4", "S5", "S6"],
    )


@pytest.fixture
def norm_metadata():
    """Metadata for norm_counts (6 samples, 2 conditions)."""
    return pd.DataFrame(
        {"condition": ["control", "control", "control", "treated", "treated", "treated"]},
        index=["S1", "S2", "S3", "S4", "S5", "S6"],
    )


class TestVST:
    @needs_deseq2
    def test_vst_basic(self, norm_counts):
        """VST returns DataFrame with same shape, values are transformed."""
        result = vst(norm_counts)

        assert isinstance(result, pd.DataFrame)
        assert result.shape == norm_counts.shape
        assert list(result.index) == list(norm_counts.index)
        assert list(result.columns) == list(norm_counts.columns)
        # Transformed values should not be raw integers (they're continuous)
        assert not (result.values == result.values.astype(int)).all()

    @needs_deseq2
    def test_vst_blind_false(self, norm_counts, norm_metadata):
        """VST with blind=False and metadata works without error."""
        result = vst(norm_counts, metadata=norm_metadata, design="~ condition", blind=False)

        assert isinstance(result, pd.DataFrame)
        assert result.shape == norm_counts.shape

    def test_vst_empty_error(self):
        """VST raises RDataError on empty input."""
        empty = pd.DataFrame()
        with pytest.raises(RDataError, match="empty"):
            vst(empty)


class TestRlog:
    @needs_deseq2
    def test_rlog_basic(self, norm_counts):
        """rlog returns DataFrame with same shape, values are transformed."""
        result = rlog(norm_counts)

        assert isinstance(result, pd.DataFrame)
        assert result.shape == norm_counts.shape
        assert list(result.index) == list(norm_counts.index)
        assert list(result.columns) == list(norm_counts.columns)
        # Transformed values should not be raw integers
        assert not (result.values == result.values.astype(int)).all()


class TestTMM:
    @needs_edger
    def test_tmm_basic(self, norm_counts):
        """TMM returns DataFrame with same shape, values are log-scale."""
        result = tmm_normalize(norm_counts)

        assert isinstance(result, pd.DataFrame)
        assert result.shape == norm_counts.shape
        assert list(result.index) == list(norm_counts.index)
        assert list(result.columns) == list(norm_counts.columns)
        # log-CPM values should contain negatives (log2 of fractions < 1)
        # and generally be much smaller than raw counts
        assert result.values.max() < norm_counts.values.max()

    def test_tmm_negative_error(self):
        """TMM raises RDataError on negative values."""
        bad = pd.DataFrame(
            {"S1": [-1, 5, 10], "S2": [3, 4, 5]},
            index=["G1", "G2", "G3"],
        )
        with pytest.raises(RDataError, match="negative"):
            tmm_normalize(bad)
