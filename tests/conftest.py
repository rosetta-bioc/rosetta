"""Shared test fixtures."""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def sample_counts():
    """3 genes x 4 samples count matrix."""
    np.random.seed(42)
    data = np.random.randint(0, 1000, size=(3, 4))
    return pd.DataFrame(data, index=["GeneA", "GeneB", "GeneC"], columns=["S1", "S2", "S3", "S4"])


@pytest.fixture
def sample_metadata():
    """Sample metadata with condition column."""
    return pd.DataFrame({"condition": ["control", "control", "treated", "treated"]}, index=["S1", "S2", "S3", "S4"])


@pytest.fixture
def limma_fit_object():
    """Provides a standard fitted limma object with contrast applied."""
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
    except Exception as e:
        pytest.skip(f"rpy2 not available: {e}")
    from rosetta.stats.design import build_contrast_matrix

    limma = importr("limma")
    
    # 1. Setup Data & Design
    data = ro.r.matrix(ro.FloatVector([10.0, 2.0, 5.0, 8.0, 3.0, 9.0, 1.0, 4.0, 6.0, 7.0, 2.0, 8.0]), 
                       nrow=3, ncol=4)
    design = ro.r.matrix(ro.FloatVector([1, 1, 1, 1, 0, 0, 1, 1]), nrow=4, ncol=2)
    colnames = ["Intercept", "Condition"]
    design.colnames = ro.StrVector(colnames)
    
    # 2. Pipeline with Contrast
    fit = limma.lmFit(data, design)
    
    # Apply contrast here to ensure the object only has 1 column of results
    contrast_mat = build_contrast_matrix(colnames, "Condition")
    fit = limma.contrasts_fit(fit, contrast_mat)
    fit = limma.eBayes(fit)
    
    return fit

@pytest.fixture
def custom_term2gene():
    """Provides a standard custom gene set mapping table for testing."""
    return pd.DataFrame({
        "term": ["Pathway_A", "Pathway_A", "Pathway_A", "Pathway_B"],
        "gene": ["TP53", "BRCA1", "EGFR", "PTEN"]
    })

@pytest.fixture
def sample_genes():
    """Provides a standard list of gene identifiers for testing."""
    return ["TP53", "BRCA1", "EGFR"]

@pytest.fixture
def deseq2_result_df():
    """Provides a standard DESeq2-like results DataFrame for GSEA testing."""
    return pd.DataFrame({
        'gene': ['GeneA', 'GeneB', 'GeneC', 'GeneD'],
        'log2FoldChange': [1.5, -2.0, 0.5, 3.0],
        'padj': [0.01, 0.05, 0.2, 0.001]
    })