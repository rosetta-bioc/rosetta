import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from typing import Optional, List, Dict, Any

from .._bridge import _converter, to_r_matrix, to_r_dataframe, to_pandas
from .._deps import ensure_installed
from .._errors import RDataError

def _phyloseq_available():
    """Check if Seurat is installed in the R environment."""
    try:
        from rosetta._deps import is_installed
        return is_installed("Phyloseq")
    except Exception:
        return False

class Phyloseq:
    """Class-based wrapper for phyloseq microbiome analysis."""

    def __init__(self, otu_table: pd.DataFrame, sample_data: Optional[pd.DataFrame] = None, tax_table: Optional[pd.DataFrame] = None):
        if otu_table.empty:
            raise RDataError("otu_table must not be empty")
        if (otu_table < 0).any().any():
            raise RDataError("OTU table contains negative values")
        ensure_installed("phyloseq")
        self.ps_pkg = importr("phyloseq")
        self.methods_pkg = importr("methods")
        self.base_pkg = importr("base")
        self.obj = self._create_object(otu_table, sample_data, tax_table)

    def _create_object(self, otu_table: pd.DataFrame, sample_data, tax_table):
        """Build the R phyloseq object from validated inputs."""
        r_otu = to_r_matrix(otu_table)
        with localconverter(_converter):
            components = [self.ps_pkg.otu_table(r_otu, taxa_are_rows=True)]
            if sample_data is not None:
                components.append(self.ps_pkg.sample_data(to_r_dataframe(sample_data)))
            if tax_table is not None:
                components.append(self.ps_pkg.tax_table(to_r_matrix(tax_table)))
            return self.ps_pkg.phyloseq(*components)

    def estimate_richness(self, measures: Optional[List[str]] = None) -> pd.DataFrame:
        """Calculate alpha diversity metrics."""
        with localconverter(_converter):
            kwargs = {"measures": ro.StrVector(measures)} if measures else {}
            return to_pandas(self.ps_pkg.estimate_richness(self.obj, **kwargs))

    def run_ordination(self, method: str = "PCoA", distance: str = "bray", **kwargs) -> pd.DataFrame:
        """Perform ordination analysis and return coordinates as a pandas DataFrame."""
        with localconverter(_converter):
            ord_obj = self.ps_pkg.ordinate(self.obj, method=method, distance=distance, **kwargs)
            
            # Extract the vectors (coordinates)
            vectors = ord_obj.rx2("vectors")
            
            # Force conversion: Ensure the matrix is seen as a data.frame before pandas conversion
            # Using base::as.data.frame ensures to_pandas() receives a standard tabular structure
            df_r = self.base_pkg.as_data_frame(vectors)
            
            return to_pandas(df_r)
