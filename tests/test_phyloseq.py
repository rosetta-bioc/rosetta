"""Tests for rosetta.wrappers.phyloseq."""

import numpy as np
import pandas as pd
import pytest

from rosetta._errors import RDataError
from rosetta.wrappers.phyloseq import Phyloseq


def _phyloseq_available():
    """Check if phyloseq is installed in the R environment."""
    try:
        from rosetta._deps import is_installed
        return is_installed("phyloseq")
    except Exception:
        return False


@pytest.fixture
def otu_table():
    """Simulated OTU count matrix."""
    np.random.seed(42)
    data = np.random.randint(0, 500, size=(5, 4))
    return pd.DataFrame(data, index=["OTU1", "OTU2", "OTU3", "OTU4", "OTU5"], columns=["S1", "S2", "S3", "S4"])


@pytest.fixture
def sample_meta():
    """Simulated sample metadata."""
    return pd.DataFrame({"site": ["gut", "gut", "skin", "skin"]}, index=["S1", "S2", "S3", "S4"])


def test_empty_otu_raises():
    """Ensure Phyloseq initialization fails with empty OTU table."""
    with pytest.raises(RDataError, match="empty"):
        Phyloseq(pd.DataFrame())


def test_negative_otu_raises(otu_table):
    """Ensure Phyloseq initialization fails with negative values."""
    bad = otu_table.copy()
    bad.iloc[0, 0] = -1
    with pytest.raises(RDataError, match="negative"):
        Phyloseq(bad)


@pytest.mark.skipif(not _phyloseq_available(), reason="phyloseq not installed in R")
def test_phyloseq_pipeline(otu_table, sample_meta):
    """Test Phyloseq initialization and alpha diversity estimation."""
    ps = Phyloseq(otu_table, sample_meta)
    
    # Test richness estimation
    richness = ps.estimate_richness(measures=["Shannon", "Simpson"])
    assert isinstance(richness, pd.DataFrame)
    assert "Shannon" in richness.columns
    assert len(richness) == 4


@pytest.mark.skipif(not _phyloseq_available(), reason="phyloseq not installed in R")
def test_run_ordination(otu_table, sample_meta):
    """Test ordination analysis integration."""
    ps = Phyloseq(otu_table, sample_meta)
    
    # Test ordination coordinates extraction
    ordination_df = ps.run_ordination(method="PCoA", distance="bray")
    assert isinstance(ordination_df, pd.DataFrame)
    assert len(ordination_df) == 4  # Should match sample count