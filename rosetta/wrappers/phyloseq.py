import pandas as pd
from typing import Optional, List
from rosetta._bridge import ACTIVE_BACKEND, BaseWrapper, _converter, to_r_matrix, to_r_dataframe, to_pandas
from rosetta.utils.kwargs import filter_kwargs
from rosetta._deps import ensure_installed
from rosetta._errors import RDataError
from .. import codegen

# Conditionally import rpy2 components based on the active backend
if ACTIVE_BACKEND == "rpy2":
    import rpy2.robjects as ro
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr
else:
    ro = None
    localconverter = None
    importr = None

def _phyloseq_available():
    """Check if phyloseq is installed in the R environment."""
    try:
        from rosetta._deps import is_installed
        return is_installed("phyloseq")
    except Exception:
        return False

class Phyloseq(BaseWrapper):
    """
    Class-based wrapper for phyloseq microbiome analysis.
    
    This class wraps the S4 phyloseq object, providing methods for alpha diversity
    estimation and ordination analysis.
    """

    _PARAMS_ORDINATION = {"method", "distance", "formula", "first", "last", "trymax", "k"}

    def __init__(self, otu_table: pd.DataFrame, sample_data: Optional[pd.DataFrame] = None, tax_table: Optional[pd.DataFrame] = None):
        if otu_table.empty:
            raise RDataError("otu_table must not be empty")
        if (otu_table < 0).any().any():
            raise RDataError("OTU table contains negative values")
            
        ensure_installed("phyloseq")
        self.ps_pkg = importr("phyloseq")
        self.base_pkg = importr("base")
        
        obj = self._create_object(otu_table, sample_data, tax_table)
        super().__init__(obj, self.ps_pkg)

    def _create_object(self, otu_table, sample_data, tax_table):
        """Build the R phyloseq object from validated inputs."""
        with localconverter(_converter):
            codegen._emit("otu <- otu_table(otu_table, taxa_are_rows=TRUE)")
            components = [self.ps_pkg.otu_table(to_r_matrix(otu_table), taxa_are_rows=True)]
            if sample_data is not None:
                codegen._emit("samples <- sample_data(sample_data)")
                components.append(self.ps_pkg.sample_data(to_r_dataframe(sample_data)))
            if tax_table is not None:
                codegen._emit("taxonomy <- tax_table(tax_table)")
                components.append(self.ps_pkg.tax_table(to_r_matrix(tax_table)))
            codegen._emit("ps <- phyloseq(otu, ...)")
            return self.ps_pkg.phyloseq(*components)

    def estimate_richness(self, measures: Optional[List[str]] = None) -> pd.DataFrame:
        """Calculate alpha diversity metrics."""
        with localconverter(_converter):
            kwargs = {"measures": ro.StrVector(measures)} if measures else {}
            codegen._emit("richness <- estimate_richness(ps, ...)")
            res = self.ps_pkg.estimate_richness(self.obj, **kwargs)
            return to_pandas(res)

    def run_ordination(self, **kwargs) -> pd.DataFrame:
        """
        Perform ordination analysis and return coordinates as a pandas DataFrame.
        
        Args:
            **kwargs: Arguments for phyloseq::ordinate (method, distance, etc.)
        """
        r_kwargs = filter_kwargs(kwargs, self._PARAMS_ORDINATION)
        
        with localconverter(_converter):
            codegen._emit("ordination <- ordinate(ps, ...)")
            ord_obj = self.ps_pkg.ordinate(self.obj, **r_kwargs)
            
            # Extract coordinates and ensure conversion to a standard R data.frame
            vectors = ord_obj.rx2("vectors")
            df_r = self.base_pkg.as_data_frame(vectors)
            
            return to_pandas(df_r)
