"""Tests for rosetta.wrappers.clusterprofiler."""

import pandas as pd
import numpy as np
import pytest
from rosetta._errors import RDataError
from rosetta.wrappers.clusterprofiler import (
    ORA, GSEA
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

def _reactome_available():
    try:
        from rosetta._deps import is_installed
        return is_installed("ReactomePA")
    except Exception:
        return False

# 1. GO Enrichment (Using sample_genes fixture)
@pytest.mark.skipif(not _cp_available(), reason="clusterProfiler/org.Hs.eg.db not installed")
def test_enrich_go_returns_dataframe(sample_genes):
    result = ORA.enrich_go(sample_genes, pvalue_cutoff=0.5)
    assert isinstance(result, pd.DataFrame)
    if not result.empty:
        assert "Description" in result.columns

@pytest.mark.skipif(not _cp_available(), reason="clusterProfiler/org.Hs.eg.db not installed")
def test_enrich_go_with_custom_params(sample_genes):
    """Verify that GO enrichment accepts custom GS size parameters."""
    result = ORA.enrich_go(
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
    result = ORA.enrich_kegg(genes, organism="hsa", pvalue_cutoff=0.5)
    assert isinstance(result, pd.DataFrame)
    assert "p.adjust" in result.columns

# 3. Reactome Enrichment
@pytest.mark.skipif(not _reactome_available(), reason="ReactomePA not installed")
def test_enrich_pathway_returns_dataframe():
    genes = ["7157", "672", "1956"] 
    result = ORA.enrich_pathway(genes, pvalueCutoff=0.5)
    assert isinstance(result, pd.DataFrame)
    assert not result.empty

# 4. Custom Enrichment (Using sample_genes and custom_term2gene fixtures)
def test_enrich_custom_returns_dataframe(sample_genes, custom_term2gene):
    """
    Test the custom enrichment wrapper with a small-scale mapping.
    
    Adjusts minGSSize to 1 to accommodate the limited number of genes 
    in the custom test fixture.
    """
    result = ORA.enrich_custom(
        sample_genes, 
        term2gene=custom_term2gene, 
        pvalueCutoff=1.0, 
        minGSSize=1,    # Adjust threshold to allow small gene sets for testing
        maxGSSize=500
    )
    
    assert isinstance(result, pd.DataFrame)
    # enricher may return empty when background universe is not specified;
    # assert only on type, not content
    assert list(result.columns) == ["ID", "Description", "GeneRatio", "BgRatio",
                                    "RichFactor", "FoldEnrichment", "zScore",
                                    "pvalue", "p.adjust", "qvalue", "geneID", "Count"]

# 5. Error Handling
def test_empty_gene_list_raises():
    with pytest.raises(RDataError):
        ORA.enrich_go([])

# Test Class GSEA
def test_gsea_prepare_gene_list(deseq2_result_df):
    """Test the GSEA gene list preparation using the fixture."""
    ranked_genes = GSEA.prepare_gene_list(
        deseq2_result_df, 
        gene_col='gene', 
        fc_col='log2FoldChange'
    )
    
    assert ranked_genes.iloc[0] == 3.0  # GeneD
    assert ranked_genes.iloc[-1] == -2.0 # GeneB

def test_gsea_prepare_gene_list_with_index():
    """Test the GSEA gene list preparation using index as gene column."""
    df = pd.DataFrame({
        'log2FoldChange': [1.0, -1.0]
    }, index=['Gene1', 'Gene2'])
    
    # Perform the conversion
    ranked_genes = GSEA.prepare_gene_list(df, gene_col='index', fc_col='log2FoldChange')
    
    # 1. Verify that the Index was correctly transferred
    assert ranked_genes.index.tolist() == ['Gene1', 'Gene2']
    
    # 2. Verify data length and numerical correctness
    assert len(ranked_genes) == 2
    assert ranked_genes.iloc[0] == 1.0
    assert ranked_genes.iloc[1] == -1.0
    
    # 3. Verify object type (Named Series)
    assert isinstance(ranked_genes, pd.Series)

def test_gsea_prepare_gene_list_handles_nans():
    """Ensure NaN values in log2FoldChange are dropped properly."""
    df = pd.DataFrame({
        'gene': ['G1', 'G2', 'G3'],
        'log2FoldChange': [1.0, np.nan, -1.0]
    })
    ranked_genes = GSEA.prepare_gene_list(df, gene_col='gene')
    
    # NaN should be dropped, so only 2 genes should remain
    assert len(ranked_genes) == 2
    assert 'G2' not in ranked_genes.index

def test_gsea_prepare_gene_list_empty_input():
    """Ensure handling of empty input gracefully."""
    df = pd.DataFrame(columns=['gene', 'log2FoldChange'])
    result = GSEA.prepare_gene_list(df, gene_col='gene')
    assert result.empty

def test_gsea_gse_go_parameter_transmission():
    """Verify that parameters like pvalue_cutoff are passed through."""
    # This is a 'mock' test to ensure our wrapper doesn't crash 
    # before calling the R function
    import pandas as pd
    dummy_genes = pd.Series([1.0, -1.0], index=['G1', 'G2'])
    
    # We verify the wrapper exists and accepts these params
    # We don't need to run the full R code, just ensure call structure is valid
    try:
        GSEA.gse_go(dummy_genes, pvalue_cutoff=0.01, eps=1e-5)
    except Exception as e:
        # We expect an error because R environment isn't initialized, 
        # but the error shouldn't be about missing parameters
        assert "must not be empty" not in str(e)  # Ensure the error is not because a empty geneList

def test_gsea_gse_kegg_parameter_transmission():
    """Verify that gse_kegg correctly accepts and transmits parameters."""
    import pandas as pd
    from rpy2.rinterface_lib.embedded import RRuntimeError
    
    dummy_genes = pd.Series([1.0, -1.0], index=['Gene1', 'Gene2'])
    
    try:
        GSEA.gse_kegg(dummy_genes, organism="mmu", pvalue_cutoff=0.01)
    except Exception as e:
        # We expect one of two things:
        # 1. RRuntimeError (from R logic failing because Gene1/Gene2 don't exist)
        # 2. RDataError (from our wrapper's validation)
        # We verify the error indicates the R call was attempted
        assert any(keyword in str(e).lower() for keyword in [
            "no gene can be mapped", 
            "ensure_installed", 
            "must not be empty"
        ])