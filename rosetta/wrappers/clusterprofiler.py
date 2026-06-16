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

    # Define a mapping to automatically handle parameter name differences in R.
    # ORA functions typically use 'gene', while GSEA functions use 'geneList'.
    if func_name.startswith("gse"):
        # For GSEA, ensure the input is passed as 'geneList'
        if "gene" in kwargs:
            kwargs["geneList"] = kwargs.pop("gene")
    
    # check if gene is empty
    gene_arg = kwargs.get("gene") or kwargs.get("geneList")
    if gene_arg is None:
        raise RDataError("gene_list must not be empty")
    
    # Dynamically retrieve the R function from the clusterProfiler package
    r_func = getattr(pkg, func_name)

    with localconverter(_converter):
        # Handle automatic conversion of Python types (List, dict, etc.) to R types
        result = r_func(**kwargs)

    return to_pandas(to_r_df(result))

class ORA:
    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def enrich_pathway(gene_list: List[str], **kwargs) -> pd.DataFrame:
        """Run Reactome pathway enrichment analysis."""
        return _run_r_enrichment("enrichPathway", library="ReactomePA", gene=gene_list, **kwargs)

    @staticmethod
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


class GSEA:
    @staticmethod
    def prepare_gene_list(df: pd.DataFrame, gene_col: str, fc_col: str = "log2FoldChange") -> pd.Series:
        """
        Convert DESeq2 results into a sorted named numeric vector for GSEA.
        
        Args:
            df: DataFrame from get_results() (e.g., index contains gene symbols).
            gene_col: Column name for gene IDs (or 'index' if using the dataframe index).
            fc_col: The log2FoldChange column name.
        """
        # 1. Handle if gene IDs are in the index
        if gene_col == "index":
            df = df.reset_index()
            # After resetting, the index column is now named 'index' 
            # (or whatever the original index name was)
            gene_col = "index" 
            
        # 2. Filter out missing values
        df = df.dropna(subset=[fc_col])
        
        # 3. Sort by Fold Change descending (required for GSEA)
        df = df.sort_values(by=fc_col, ascending=False)
        
        # 4. Create named series
        gene_list = df.set_index(gene_col)[fc_col]
        return gene_list
    
    @staticmethod
    def gse_go(
        gene_list: pd.Series, 
        organism: str = "org.Hs.eg.db", 
        ont: str = "BP", 
        pvalue_cutoff: float = 0.05, 
        eps: float = 1e-10, 
        **kwargs
    ) -> pd.DataFrame:
        """Run GSEA using Gene Ontology with significance filtering."""
        kwargs.update({
            "OrgDb": organism,
            "ont": ont,
            "pvalueCutoff": pvalue_cutoff,
            "eps": eps  # Adding eps to handle numerical stability
        })
        return _run_r_enrichment("gseGO", geneList=gene_list, **kwargs)

    @staticmethod
    def gse_kegg(
        gene_list: pd.Series, 
        organism: str = "hsa", 
        pvalue_cutoff: float = 0.05, 
        eps: float = 1e-10, 
        **kwargs
    ) -> pd.DataFrame:
        """Run GSEA using KEGG pathways with significance filtering."""
        kwargs.update({
            "organism": organism,
            "pvalueCutoff": pvalue_cutoff,
            "eps": eps
        })
        return _run_r_enrichment("gseKEGG", geneList=gene_list, **kwargs)