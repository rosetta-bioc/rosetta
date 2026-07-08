import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from rosetta.utils import filter_kwargs
from rosetta._errors import RDataError
from typing import Dict, Any

from .._bridge import _converter, to_r_matrix, to_pandas, to_r_df, BaseWrapper
from .._deps import ensure_installed
from .._errors import RDataError

def _seurat_available():
    """Check if Seurat is installed in the R environment."""
    try:
        from rosetta._deps import is_installed
        return is_installed("Seurat")
    except Exception:
        return False

class Seurat(BaseWrapper):
    """Class-based wrapper for Seurat single-cell analysis pipeline."""
    # 1. Class-level parameter allowlists
    _PARAMS_NORMALIZE = {"normalization.method", "scale.factor", "verbose"}
    _PARAMS_VARIABLE_FEATURES = {"selection.method", "nfeatures", "verbose"}
    _PARAMS_SCALE_DATA = {"features", "vars.to.regress", "verbose"}
    _PARAMS_PCA = {"npcs", "features", "verbose", "seed.use"}
    _PARAMS_FIND_NEIGHBORS = {"dims", "reduction", "k.param", "verbose"}
    _PARAMS_FIND_CLUSTERS = {"resolution", "algorithm", "verbose", "random.seed"}

    def __init__(self, counts: pd.DataFrame):
        if counts.empty:
            raise RDataError("counts must not be empty")
        if (counts < 0).any().any():
            raise RDataError("Count matrix contains negative values")
        
        ensure_installed("Seurat")
        seurat_pkg = importr("Seurat")
        self.sobj_pkg = importr("SeuratObject")
        self.methods_pkg = importr("methods")
        
        obj = self._create_object(counts)
        
        super().__init__(obj, seurat_pkg)
        self.seurat_pkg = seurat_pkg

    # 3. Business logic methods
    def run_umap(self, **kwargs):
        # Apply filter_kwargs to ensure parameter validity
        r_kwargs = filter_kwargs(kwargs, self._PARAMS_RUN_UMAP)
        
        # Execute UMAP with filtered parameters
        self.obj = self.seurat_pkg.RunUMAP(self.obj, **r_kwargs)
        return self

    def run_sctransform(self, **kwargs):
        # Apply filter_kwargs to ensure parameter validity
        r_kwargs = filter_kwargs(kwargs, self._PARAMS_SCT)
        
        # Execute SCTransform
        self.obj = self.seurat_pkg.SCTransform(self.obj, **r_kwargs)
        return self
    
    def run_normalize(self, **kwargs):
        return self._call_r("NormalizeData", self._PARAMS_NORMALIZE, **kwargs)

    def run_find_variable_features(self, **kwargs):
        r_kwargs = filter_kwargs(kwargs, self._PARAMS_VARIABLE_FEATURES)
        self.obj = self.seurat_pkg.FindVariableFeatures(self.obj, **r_kwargs)
        return self

    def run_scale_data(self, **kwargs):
        r_kwargs = filter_kwargs(kwargs, self._PARAMS_SCALE_DATA)
        self.obj = self.seurat_pkg.ScaleData(self.obj, **r_kwargs)
        return self

    def run_scale_data(self, **kwargs):
        """Scale the data before PCA."""
        r_kwargs = filter_kwargs(kwargs, self._PARAMS_SCALE_DATA)
        self.obj = self.seurat_pkg.ScaleData(self.obj, **r_kwargs)
        return self

    def run_pca(self, **kwargs):
        """Run PCA on the Seurat object."""
        return self._call_r("RunPCA", self._PARAMS_PCA, **kwargs)

    def run_find_neighbors(self, **kwargs):
        """Find neighbors for clustering."""
        if "k_param" in kwargs:
            kwargs["k.param"] = kwargs.pop("k_param")
        return self._call_r("FindNeighbors", self._PARAMS_FIND_NEIGHBORS, **kwargs)

    def run_find_clusters(self, **kwargs):
        """Run graph-based clustering."""
        r_kwargs = filter_kwargs(kwargs, self._PARAMS_FIND_CLUSTERS)
        self.obj = self.seurat_pkg.FindClusters(self.obj, **r_kwargs)
        return self

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
