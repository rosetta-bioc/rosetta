"""limma-voom differential expression wrapper."""
import pandas as pd
from rosetta._bridge import ACTIVE_BACKEND, BaseWrapper, _converter, to_r_matrix, to_r_dataframe, to_pandas, r_nrow
from rosetta.utils.kwargs import filter_kwargs
from rosetta._deps import ensure_installed
from rosetta._errors import RDataError, RFormulaError
from rosetta.stats.design import build_contrast_matrix
from rosetta.stats.decide import run_decide_tests
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

class Limma(BaseWrapper):
    """Class-based wrapper for Limma-Voom differential expression analysis."""
    _codegen_target = "fit"

    _PARAMS_VOOMLMFIT = {"block", "correlation", "weights", "sample.weights", "span", "plot"}
    _PARAMS_EBAYES = {"trend", "robust", "proportion", "winsor.tail.p"}
    _PARAMS_TOPTABLE = {"coef", "number", "adjust.method", "p.value", "lfc"}

    def __init__(self, counts: pd.DataFrame, metadata: pd.DataFrame, design: str = "~ condition"):
        ensure_installed("limma")
        ensure_installed("edgeR")
        
        self.limma_pkg = importr("limma")
        self.edger_pkg = importr("edgeR")
        stats_pkg = importr("stats")

        # 1. Fitting logic
        obj = self._fit_model(counts, metadata, design, stats_pkg)
        
        # 2. call BaseWrapper
        super().__init__(obj, self.limma_pkg)

    def _fit_model(self, counts, metadata, design, stats_pkg):
        """Encapsulated model fitting logic."""
        if (counts < 0).any().any():
            raise RDataError("Count matrix contains negative values")
        if not set(counts.columns).issubset(set(metadata.index)):
            raise RDataError("Count matrix columns must match metadata row names")

        try:
            r_design_formula = ro.Formula(design)
            r_counts = to_r_matrix(counts)
            r_metadata = to_r_dataframe(metadata)
            r_design_matrix = stats_pkg.model_matrix(r_design_formula, data=r_metadata)
        except Exception as e:
            raise RFormulaError(f"Invalid design formula '{design}': {e}")

        with localconverter(_converter):
            dge = self.edger_pkg.DGEList(counts=r_counts)
            codegen._emit("dge <- DGEList(counts=counts)")
            dge = self.edger_pkg.calcNormFactors(dge)
            codegen._emit("dge <- calcNormFactors(dge)")
            # voomLmFit
            codegen._emit("fit <- voomLmFit(dge, design)")
            return self.edger_pkg.voomLmFit(dge, r_design_matrix)

    def apply_contrasts(self, contrast: list):
        """Apply contrast matrix to fitted model."""
        with localconverter(_converter):
            design_matrix = self.obj.rx2("design") 
            design_colnames = design_matrix.colnames
            
            contrast_mat = build_contrast_matrix(design_colnames, contrast)
            codegen._emit("fit <- contrasts.fit(fit, contrast.matrix)")
            self.obj = self.limma_pkg.contrasts_fit(self.obj, contrast_mat)
        return self
    
    def run_ebayes(self, **kwargs):
        """Perform empirical Bayes moderation."""
        # use call_r to process eBayes
        return self._call_r("eBayes", self._PARAMS_EBAYES, **kwargs)

    def get_results(self, **kwargs) -> pd.DataFrame:
        """Extract results using topTable."""
        r_kwargs = filter_kwargs(kwargs, self._PARAMS_TOPTABLE)
        if "number" not in r_kwargs:
            r_kwargs["number"] = r_nrow(self.obj)
            
        with localconverter(_converter):
            codegen._emit("res <- topTable(fit, ...)")
            res = self.limma_pkg.topTable(self.obj, **r_kwargs)
            return to_pandas(res)
