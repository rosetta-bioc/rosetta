"""clusterProfiler gene set enrichment wrapper."""

import pandas as pd
from typing import Any
from rosetta._bridge import ACTIVE_BACKEND, BaseWrapper, _converter, to_pandas, to_r_df
from rosetta.utils.kwargs import filter_kwargs
from rosetta._deps import ensure_installed
from rosetta._errors import RDataError

# Conditionally import rpy2 components based on the active backend
if ACTIVE_BACKEND == "rpy2":
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter
else:
    importr = None
    localconverter = None

class ClusterProfiler(BaseWrapper):
    """
    Class-based wrapper for clusterProfiler enrichment analysis.
    Provides a standardized interface for ORA and GSEA workflows.
    """

    def __init__(self, library: str = "clusterProfiler"):
        """Initialize the wrapper with the required R library."""
        ensure_installed(library)
        self.cp_pkg = importr(library)
        # Inherit from BaseWrapper; obj is set to None for stateless enrichment tools
        super().__init__(None, self.cp_pkg)

    def _run_enrich(self, func_name: str, gene_data: Any, **kwargs) -> pd.DataFrame:
        """Internal execution method with automatic parameter mapping."""
        if gene_data is None or (isinstance(gene_data, (list, pd.Series)) and len(gene_data) == 0):
            raise RDataError("Input gene list is empty")
        
        # 1. Parameter Mapping (Python snake_case -> R CamelCase)
        # This fixes the unused argument errors
        param_map = {
            "pvalue_cutoff": "pvalueCutoff",
            "qvalue_cutoff": "qvalueCutoff",
            "min_gs_size": "minGSSize",
            "max_gs_size": "maxGSSize"
        }
        for snake, camel in param_map.items():
            if snake in kwargs:
                kwargs[camel] = kwargs.pop(snake)

        # 2. Key Mapping
        key = "geneList" if func_name.startswith("gse") else "gene"
        kwargs[key] = gene_data
        
        # 3. Execution
        with localconverter(_converter):
            r_func = getattr(self.cp_pkg, func_name)
            res = r_func(**kwargs)
            return to_pandas(to_r_df(res))

    # --- ORA (Over-Representation Analysis) Interface ---

    def enrich_go(self, gene_list: Any, organism: str = "org.Hs.eg.db", **kwargs) -> pd.DataFrame:
        """Perform Gene Ontology (GO) enrichment analysis."""
        kwargs.update({"OrgDb": organism})
        return self._run_enrich("enrichGO", gene_list, **kwargs)

    def enrich_kegg(self, gene_list: Any, organism: str = "hsa", **kwargs) -> pd.DataFrame:
        """Perform KEGG pathway enrichment analysis."""
        kwargs.update({"organism": organism})
        return self._run_enrich("enrichKEGG", gene_list, **kwargs)

    # --- GSEA (Gene Set Enrichment Analysis) Interface ---

    def gse_go(self, gene_list: Any, organism: str = "org.Hs.eg.db", **kwargs) -> pd.DataFrame:
        """Perform GSEA using Gene Ontology."""
        kwargs.update({"OrgDb": organism})
        return self._run_enrich("gseGO", gene_list, **kwargs)

    @staticmethod
    def prepare_gene_list(df: pd.DataFrame, gene_col: str, fc_col: str = "log2FoldChange") -> pd.Series:
        """
        Convert differential expression results into a sorted named numeric vector for GSEA.
        
        Args:
            df: DataFrame containing log2FoldChange values.
            gene_col: Column name containing gene IDs.
            fc_col: Column name for fold change values.
        """
        if gene_col == "index":
            df = df.reset_index()
            gene_col = "index"
            
        df = df.dropna(subset=[fc_col])
        df = df.sort_values(by=fc_col, ascending=False)
        return df.set_index(gene_col)[fc_col]