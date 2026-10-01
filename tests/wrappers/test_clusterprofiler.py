"""Tests for rosetta.wrappers.clusterprofiler."""

import numpy as np
import pandas as pd
import pytest

from rosetta._deps import is_installed
from rosetta._errors import RDataError
from rosetta.wrappers.clusterprofiler import ClusterProfiler

# --- Fixtures ---

@pytest.fixture
def cp_model():
    """Fixture to provide an initialized ClusterProfiler instance."""
    if not is_installed("clusterProfiler"):
        pytest.skip("clusterProfiler R package not installed")
    return ClusterProfiler()

# --- GO Enrichment ---

def test_enrich_go_returns_dataframe(cp_model, sample_genes):
    result = cp_model.enrich_go(sample_genes, pvalue_cutoff=0.5)
    assert isinstance(result, pd.DataFrame)
    if not result.empty:
        assert "Description" in result.columns

def test_enrich_go_with_custom_params(cp_model, sample_genes):
    """Verify that GO enrichment accepts custom GS size parameters."""
    result = cp_model.enrich_go(
        sample_genes,
        pvalue_cutoff=1.0,
        minGSSize=5,
        maxGSSize=1000
    )
    assert isinstance(result, pd.DataFrame)

# --- KEGG Enrichment ---

def test_enrich_kegg_returns_dataframe(cp_model):
    # KEGG usually requires ENTREZIDs
    genes = ["7157", "672", "1956", "4609", "5728", "5925", "207", "3845", "673", "5290"]
    result = cp_model.enrich_kegg(genes, organism="hsa", pvalue_cutoff=0.5)
    assert isinstance(result, pd.DataFrame)
    assert "p.adjust" in result.columns

# --- Error Handling ---

def test_empty_gene_list_raises(cp_model):
    with pytest.raises(RDataError, match="empty"):
        cp_model.enrich_go([])

# --- GSEA & Data Preparation (Static Methods) ---

def test_gsea_prepare_gene_list(deseq2_result_df):
    """Test the GSEA gene list preparation."""
    ranked_genes = ClusterProfiler.prepare_gene_list(
        deseq2_result_df,
        gene_col='gene',
        fc_col='log2FoldChange'
    )

    assert ranked_genes.iloc[0] == 3.0  # GeneD
    assert ranked_genes.iloc[-1] == -2.0 # GeneB
    assert isinstance(ranked_genes, pd.Series)

def test_gsea_prepare_gene_list_with_index():
    df = pd.DataFrame({'log2FoldChange': [1.0, -1.0]}, index=['Gene1', 'Gene2'])
    ranked_genes = ClusterProfiler.prepare_gene_list(df, gene_col='index', fc_col='log2FoldChange')

    assert ranked_genes.index.tolist() == ['Gene1', 'Gene2']
    assert len(ranked_genes) == 2


@pytest.mark.parametrize("index_name", [None, "gene_id", "index", "log2FoldChange"])
@pytest.mark.parametrize("has_index_column", [False, True])
def test_gsea_prepare_gene_list_preserves_row_ids(index_name, has_index_column):
    df = pd.DataFrame(
        {"log2FoldChange": [1.0, np.nan, 2.0, -1.0]},
        index=pd.Index(["GeneA", "GeneB", "GeneC", "GeneD"], name=index_name),
    )
    if has_index_column:
        df["index"] = ["not_a_gene_1", "not_a_gene_2", "not_a_gene_3", "not_a_gene_4"]
    original = df.copy(deep=True)

    ranked_genes = ClusterProfiler.prepare_gene_list(df, gene_col="index")

    expected = pd.Series(
        [2.0, 1.0, -1.0],
        index=pd.Index(["GeneC", "GeneA", "GeneD"], name=index_name),
        name="log2FoldChange",
    )
    pd.testing.assert_series_equal(ranked_genes, expected)
    pd.testing.assert_frame_equal(df, original)


def test_gsea_prepare_gene_list_handles_nans():
    df = pd.DataFrame({
        'gene': ['G1', 'G2', 'G3'],
        'log2FoldChange': [1.0, np.nan, -1.0]
    })
    ranked_genes = ClusterProfiler.prepare_gene_list(df, gene_col='gene')
    assert len(ranked_genes) == 2
    assert 'G2' not in ranked_genes.index
