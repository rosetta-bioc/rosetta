"""Tests for rosetta.wrappers.phyloseq."""

"""Tests for rosetta.wrappers.phyloseq."""

import numpy as np
import pandas as pd
import pytest

from rosetta._errors import RDataError
from rosetta.wrappers.phyloseq import Phyloseq

# --- Fixtures ---

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

@pytest.fixture
def ps_model(otu_table, sample_meta):
    """Fixture to provide an initialized Phyloseq model."""
    return Phyloseq(otu_table, sample_data=sample_meta)

# --- Initialization Tests ---

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

# --- Analysis Tests ---

def test_phyloseq_pipeline(ps_model):
    """Test Phyloseq alpha diversity estimation."""
    richness = ps_model.estimate_richness(measures=["Shannon", "Simpson"])
    assert isinstance(richness, pd.DataFrame)
    assert "Shannon" in richness.columns
    assert len(richness) == 4

def test_run_ordination(ps_model):
    """Test ordination analysis integration.""" 
    ordination_df = ps_model.run_ordination(method="PCoA", distance="bray")
    assert isinstance(ordination_df, pd.DataFrame)
    assert len(ordination_df) == 4