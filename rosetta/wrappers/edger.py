"""edgeR differential expression wrapper."""

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from rosetta._bridge import BaseWrapper, _converter, to_r_matrix, to_r_dataframe, to_pandas, to_r_df, r_nrow
from rosetta.utils.kwargs import filter_kwargs
from rosetta._deps import ensure_installed
from rosetta._errors import RDataError, RFormulaError

class EdgeR(BaseWrapper):
    """Class-based wrapper for edgeR quasi-likelihood analysis."""

    _PARAMS_QLFIT = {"dispersion", "robust", "winsor.tail.p", "abundance.trend"}
    _PARAMS_TEST = {"contrast", "coef", "lfc"}

    def __init__(self, counts: pd.DataFrame, metadata: pd.DataFrame, design: str = "~ condition"):
        ensure_installed("edgeR")
        self.edger_pkg = importr("edgeR")
        stats_pkg = importr("stats")

        obj = self._fit_model(counts, metadata, design, stats_pkg)
        
        super().__init__(obj, self.edger_pkg)

    def _fit_model(self, counts, metadata, design, stats_pkg):
        if (counts < 0).any().any():
            raise RDataError("Count matrix contains negative values")
        if not set(counts.columns).issubset(set(metadata.index)):
            raise RDataError("Count matrix columns must match metadata row names")
        r_counts = to_r_matrix(counts)
        r_metadata = to_r_dataframe(metadata)
        
        with localconverter(_converter):
            try:
                r_design = stats_pkg.model_matrix(ro.Formula(design), data=r_metadata)
            except Exception as e:
                raise RFormulaError(f"Invalid design formula: {e}")

            dge = self.edger_pkg.DGEList(counts=r_counts)
            dge = self.edger_pkg.calcNormFactors(dge)
            dge = self.edger_pkg.estimateDisp(dge, r_design)
            return self.edger_pkg.glmQLFit(dge, r_design)

    def run_test(self, lfc: float = 0, **kwargs):
        """Perform glmTreat (if lfc > 0) or glmQLFTest."""
        r_kwargs = filter_kwargs(kwargs, self._PARAMS_TEST)
        
        with localconverter(_converter):
            if lfc > 0:
                r_kwargs["lfc"] = lfc
                res = self.edger_pkg.glmTreat(self.obj, **r_kwargs)
            else:
                res = self.edger_pkg.glmQLFTest(self.obj, **r_kwargs)
            
            return res

    def get_results(self, res_obj, **kwargs) -> pd.DataFrame:
        """Extract results using topTags."""
        with localconverter(_converter):
            top = self.edger_pkg.topTags(res_obj, n=r_nrow(self.obj.rx2("counts")), **kwargs)
            return to_pandas(to_r_df(top))