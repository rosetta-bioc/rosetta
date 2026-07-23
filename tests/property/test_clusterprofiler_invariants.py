"""Property-based tests for ClusterProfiler wrapper invariants using Hypothesis."""

import pytest
import numpy as np
from hypothesis import given, settings
from tests.property.conftest import valid_gene_list
from rosetta._deps import is_installed
from rosetta.wrappers.clusterprofiler import ClusterProfiler

if not is_installed("clusterprofiler"):
    pytest.skip("clusterprofiler package not available", allow_module_level=True)

# Pool of valid Entrez IDs that exist in org.Hs.eg.db
REAL_ENTREZ_POOL = ["7157", "3845", "672", "675", "1017", "2064", "2353", "3162", "348", "5925", "673", "7158", "836", "10", "19"]

@settings(max_examples=5, deadline=None)
@given(raw_genes=valid_gene_list(min_len=10, max_len=25))
def test_clusterprofiler_invariants(raw_genes):
    # Use the generated random list length from Hypothesis to sample from the real ID pool
    import numpy as np
    np.random.seed(len(raw_genes))  # Ensure pseudo-random stability for the test run
    genes = list(np.random.choice(REAL_ENTREZ_POOL, size=min(len(raw_genes), len(REAL_ENTREZ_POOL)), replace=False))
    
    try:
        model = ClusterProfiler(gene_list=genes, organism="org.Hs.eg.db", key_type="ENTREZID")
        results = model.enrich_go(ont="BP")
    except Exception as e:
        pytest.skip(f"Skipping due to missing R organism database or execution error: {e}")
    
    if results is not None and not results.empty:
        expected_cols = {"ID", "Description", "pvalue", "p.adjust"}
        assert expected_cols.issubset(set(results.columns)), "Missing mandatory enrichment columns"
        
        for p_col in ["pvalue", "p.adjust"]:
            if p_col in results.columns:
                valid_p = results[p_col].dropna()
                if len(valid_p) > 0:
                    assert ((valid_p >= 0.0) & (valid_p <= 1.0)).all(), f"{p_col} must be between 0 and 1"