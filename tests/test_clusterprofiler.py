"""Tests for rosetta.wrappers.clusterprofiler."""

import pandas as pd
import pytest
from rosetta._errors import RDataError
from rosetta.wrappers.clusterprofiler import (
    enrich_go, enrich_kegg, enrich_pathway, enrich_custom
)

# Helper functions for skipif logic
def _cp_available():
    try:
        from rosetta._deps import is_installed
        return is_installed("clusterProfiler") and is_installed("org.Hs.eg.db")
    except Exception:
        return False

def _kegg_available():
    try:
        from rosetta._deps import is_installed
        return is_installed("clusterProfiler")
    except Exception:
        return False

# 1. GO Enrichment (Using sample_genes fixture)
@pytest.mark.skipif(not _cp_available(), reason="clusterProfiler/org.Hs.eg.db not installed")
def test_enrich_go_returns_dataframe(sample_genes):
    result = enrich_go(sample_genes, pvalue_cutoff=0.5)
    assert isinstance(result, pd.DataFrame)
    if not result.empty:
        assert "Description" in result.columns

@pytest.mark.skipif(not _cp_available(), reason="clusterProfiler/org.Hs.eg.db not installed")
def test_enrich_go_with_custom_params(sample_genes):
    """Verify that GO enrichment accepts custom GS size parameters."""
    result = enrich_go(
        sample_genes, 
        pvalue_cutoff=1.0, 
        min_gs_size=5,  # Test custom threshold
        max_gs_size=1000
    )
    assert isinstance(result, pd.DataFrame)


# 2. KEGG Enrichment
@pytest.mark.skipif(not _kegg_available(), reason="clusterProfiler not installed")
def test_enrich_kegg_returns_dataframe():
    # KEGG usually requires ENTREZIDs, so we don't use sample_genes fixture here
    genes = ["7157", "672", "1956", "4609", "5728", "5925", "207", "3845", "673", "5290"]
    result = enrich_kegg(genes, organism="hsa", pvalue_cutoff=0.5)
    assert isinstance(result, pd.DataFrame)
    assert "p.adjust" in result.columns

# 3. Reactome Enrichment
@pytest.mark.skipif(not _kegg_available(), reason="Needs ReactomePA")
def test_enrich_pathway_returns_dataframe():
    genes = ["7157", "672", "1956"] 
    result = enrich_pathway(genes, pvalueCutoff=0.5)
    assert isinstance(result, pd.DataFrame)
    assert not result.empty

# 4. Custom Enrichment (Using sample_genes and custom_term2gene fixtures)
def test_enrich_custom_returns_dataframe(sample_genes, custom_term2gene):
    """
    Test the custom enrichment wrapper with a small-scale mapping.
    
    Adjusts minGSSize to 1 to accommodate the limited number of genes 
    in the custom test fixture.
    """
    result = enrich_custom(
        sample_genes, 
        term2gene=custom_term2gene, 
        pvalueCutoff=1.0, 
        minGSSize=1,    # Adjust threshold to allow small gene sets for testing
        maxGSSize=500
    )
    
    assert isinstance(result, pd.DataFrame)
    assert not result.empty
    assert "Pathway_A" in result["ID"].values

# 5. Error Handling
def test_empty_gene_list_raises():
    with pytest.raises(RDataError):
        enrich_go([])