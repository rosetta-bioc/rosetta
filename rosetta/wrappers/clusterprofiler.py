"""clusterProfiler gene set enrichment wrapper."""

from typing import List, Any
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from .._bridge import _converter, to_pandas, to_r_df
from .._deps import ensure_installed
from .._errors import RDataError

def _run_r_enrichment(func_name: str, library: str = "clusterProfiler", **kwargs) -> pd.DataFrame:
    """
    Internal generic wrapper to execute clusterProfiler R functions.

    Args:
        func_name: The name of the clusterProfiler function to execute.
        **kwargs: Arguments to pass to the R function.

    Returns:
        A pandas DataFrame containing the enrichment analysis results.
    """
    ensure_installed(library)
    pkg = importr(library)

    # check if gene is empty
    gene_arg = kwargs.get("gene")
    if not gene_arg:
        raise RDataError("gene_list must not be empty")
    
    # Dynamically retrieve the R function from the clusterProfiler package
    r_func = getattr(pkg, func_name)

    with localconverter(_converter):
        # Handle automatic conversion of Python types (List, dict, etc.) to R types
        result = r_func(**kwargs)

    return to_pandas(to_r_df(result))


def enrich_go(
    gene_list: List[str], 
    organism: str = "org.Hs.eg.db", 
    ont: str = "BP", 
    pvalue_cutoff: float = 0.05, 
    min_gs_size: int = 10, 
    max_gs_size: int = 500,
    **kwargs
) -> pd.DataFrame:
    """Run Gene Ontology (GO) enrichment analysis."""
    kwargs.update({
        "OrgDb": organism,
        "ont": ont,
        "pvalueCutoff": pvalue_cutoff,
        "minGSSize": min_gs_size,
        "maxGSSize": max_gs_size
    })
    return _run_r_enrichment("enrichGO", gene=gene_list, **kwargs)


def enrich_kegg(
    gene_list: List[str], 
    organism: str = "hsa", 
    pvalue_cutoff: float = 0.05, 
    min_gs_size: int = 10, 
    max_gs_size: int = 500,
    **kwargs
) -> pd.DataFrame:
    """Run KEGG pathway enrichment analysis."""
    kwargs.update({
        "organism": organism,
        "pvalueCutoff": pvalue_cutoff,
        "minGSSize": min_gs_size,
        "maxGSSize": max_gs_size
    })
    return _run_r_enrichment("enrichKEGG", gene=gene_list, **kwargs)


def enrich_pathway(gene_list: List[str], **kwargs) -> pd.DataFrame:
    """Run Reactome pathway enrichment analysis."""
    return _run_r_enrichment("enrichPathway", library="ReactomePA", gene=gene_list, **kwargs)


def enrich_custom(
    gene_list: List[str], 
    term2gene: pd.DataFrame, 
    min_gs_size: int = 10, 
    max_gs_size: int = 500, 
    **kwargs
) -> pd.DataFrame:
    """Run custom enrichment analysis."""
    # Explicitly pack parameters into a dictionary to ensure they are not missed when passed to R
    r_params = {
        "gene": gene_list,
        "TERM2GENE": term2gene,
        "minGSSize": min_gs_size,
        "maxGSSize": max_gs_size,
    }
    # Merge any additional parameters
    r_params.update(kwargs)
    
    return _run_r_enrichment("enricher", **r_params)