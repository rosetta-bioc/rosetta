"""Hypothesis strategies for property-based testing in Rosetta."""

import pandas as pd
from hypothesis import strategies as st

@st.composite
def _raw_count_matrix(draw, min_genes, max_genes, min_samples, max_samples):
    n_genes = draw(st.integers(min_value=min_genes, max_value=max_genes))
    n_samples = draw(st.integers(min_value=min_samples, max_value=max_samples))
    
    gene_names = [f"Gene{i}" for i in range(n_genes)]
    sample_names = [f"S{j}" for j in range(n_samples)]
    
    data = draw(
        st.lists(
            st.lists(st.integers(min_value=0, max_value=2000), min_size=n_samples, max_size=n_samples),
            min_size=n_genes,
            max_size=n_genes
        )
    )
    return pd.DataFrame(data, index=gene_names, columns=sample_names)

def valid_count_matrix(min_genes=5, max_genes=20, min_samples=4, max_samples=10):
    """
    Returns a strategy for a count matrix with zeros, guaranteed to have no all-zero genes.
    """
    return _raw_count_matrix(min_genes, max_genes, min_samples, max_samples).filter(
        lambda df: (df.sum(axis=1) > 0).all()
    )

@st.composite
def valid_metadata(draw, sample_names):
    """
    Generates a valid sample metadata DataFrame with at least 2 distinct groups (condition).
    """
    # Assign each sample to one of two conditions (control / treated)
    conditions = draw(
        st.lists(
            st.sampled_from(["control", "treated"]),
            min_size=len(sample_names),
            max_size=len(sample_names)
        )
    )
    
    # Ensure both levels are present to avoid design matrix singularity errors
    if len(set(conditions)) < 2:
        conditions[0] = "control"
        conditions[1] = "treated"
        
    meta = pd.DataFrame({"condition": conditions}, index=sample_names)
    return meta

@st.composite
def valid_gene_list(draw, min_len=1, max_len=50):
    """
    Generates a list of unique gene identifiers.
    """
    n = draw(st.integers(min_value=min_len, max_value=max_len))
    return [f"Gene{i}" for i in draw(st.lists(st.integers(min_value=0, max_value=10000), min_size=n, max_size=n, unique=True))]