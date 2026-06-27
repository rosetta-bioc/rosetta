import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from typing import Dict, Any

from .._bridge import _converter, to_r_matrix, to_pandas, to_r_df
from .._deps import ensure_installed
from .._errors import RDataError

class Seurat:
    """Class-based wrapper for Seurat single-cell analysis pipeline."""

    def __init__(self, counts: pd.DataFrame):
        if counts.empty:
            raise RDataError("counts must not be empty")
        if (counts < 0).any().any():
            raise RDataError("Count matrix contains negative values")
        ensure_installed("Seurat")
        self.seurat_pkg = importr("Seurat")
        self.sobj_pkg = importr("SeuratObject")
        self.base_pkg = importr("base")
        self.methods_pkg = importr("methods")
        
        self.obj = self._create_object(counts)

    def _create_object(self, counts: pd.DataFrame):
        r_counts = to_r_matrix(counts)
        with localconverter(_converter):
            return self.sobj_pkg.CreateSeuratObject(counts=r_counts)

    def run_sctransform(self, **kwargs):
        """Standardized normalization using SCTransform."""
        with localconverter(_converter):
            self.obj = self.seurat_pkg.SCTransform(self.obj, verbose=False, **kwargs)
        return self

    def run_standard_pipeline(self, n_variable_features=2000, n_pcs=10, resolution=0.5, **kwargs):
        """Replicates your original 'one-stop' pipeline logic."""
        dims = ro.IntVector(range(1, n_pcs + 1))
        
        with localconverter(_converter):
            self.obj = self.seurat_pkg.NormalizeData(self.obj, verbose=False)
            self.obj = self.seurat_pkg.FindVariableFeatures(self.obj, nfeatures=n_variable_features, verbose=False)
            self.obj = self.seurat_pkg.ScaleData(self.obj, verbose=False)
            self.obj = self.seurat_pkg.RunPCA(self.obj, npcs=n_pcs, verbose=False)
            self.obj = self.seurat_pkg.FindNeighbors(self.obj, dims=dims, verbose=False)
            self.obj = self.seurat_pkg.FindClusters(self.obj, resolution=resolution, verbose=False, **kwargs)
            self.obj = self.seurat_pkg.RunUMAP(self.obj, dims=dims, verbose=False)
        return self
    
    def find_markers(self, ident_1: str, ident_2: str = None, group_by: str = None, **kwargs) -> pd.DataFrame:
        """Find differentially expressed markers between groups."""
        with localconverter(_converter):
            if group_by:
                self.obj = ro.r("`Idents<-`")(self.obj, value=group_by)

            find_kwargs = {"ident.1": ident_1}
            if ident_2 is not None:
                find_kwargs["ident.2"] = ident_2
            find_kwargs.update(kwargs)

            markers = self.seurat_pkg.FindMarkers(self.obj, **find_kwargs)
        return to_pandas(markers)

    def get_results(self) -> Dict[str, Any]:
        """Extracts and formats the results."""
        with localconverter(_converter):
            meta = to_pandas(to_r_df(self.methods_pkg.slot(self.obj, "meta.data")))
            embeddings = to_pandas(to_r_df(self.sobj_pkg.Embeddings(self.obj, reduction="umap")))
            var_features = list(self.sobj_pkg.VariableFeatures(self.obj))
            
        return {
            "clusters": meta["seurat_clusters"],
            "umap": embeddings,
            "variable_features": var_features,
        }
