"""DESeq2 differential expression wrapper."""

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from .._bridge import BaseWrapper, _converter, to_r_matrix, to_pandas, to_r_df
from .._deps import ensure_installed
from .._errors import RDataError, RFormulaError
from rosetta.utils.kwargs import filter_kwargs
from .. import codegen

class DESeq2(BaseWrapper):
    """Class-based wrapper for DESeq2 differential expression analysis."""
    
    # 1. Parameter allowlists
    _PARAMS_DESEQ = {"fitType", "betaPrior", "parallel"}
    _PARAMS_RESULTS = {"contrast", "lfcThreshold", "alpha"}
    _PARAMS_SHRINK = {"coef", "type", "lfcThreshold", "parallel"}

    def __init__(self, counts: pd.DataFrame, metadata: pd.DataFrame, design: str):
        ensure_installed("DESeq2")
        deseq2_pkg = importr("DESeq2")
        stats_pkg = importr("stats")
        
        # Internal creation logic
        obj = self._create_dds(counts, metadata, design, stats_pkg, deseq2_pkg)
        
        # Initialize BaseWrapper
        super().__init__(obj, deseq2_pkg)
        self.deseq_pkg = deseq2_pkg

    def _create_dds(self, counts, metadata, design, stats_pkg, deseq2_pkg):
        """Encapsulated object creation logic."""
        if (counts < 0).any().any():
            raise RDataError("Count matrix contains negative values")
        if not set(counts.columns).issubset(set(metadata.index)):
            raise RDataError("Count matrix columns must match metadata row names")

        try:
            r_design = stats_pkg.as_formula(design)
        except Exception as e:
            raise RFormulaError(f"Invalid design formula '{design}': {e}")

        r_counts = to_r_matrix(counts)
        r_metadata = to_r_df(metadata)

        with localconverter(_converter):
            try:
                # Keep your original codegen logging
                codegen._block([
                    "library(DESeq2)",
                    f"dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design={design})",
                ])
                return deseq2_pkg.DESeqDataSetFromMatrix(
                    countData=r_counts, colData=r_metadata, design=r_design
                )
            except Exception as e:
                raise RDataError(f"Failed to create DESeqDataSet: {e}")

    def run_deseq(self, **kwargs):
        """Fit the DESeq2 model."""
        with localconverter(_converter):
            codegen._emit("dds <- DESeq(dds)")
            return self._call_r("DESeq", self._PARAMS_DESEQ, **kwargs)
    
    def lfc_shrink(self, **kwargs):
        """Perform LFC shrinkage with manual parameter validation."""
        r_kwargs = filter_kwargs(kwargs, self._PARAMS_SHRINK)

        if "type" in r_kwargs and r_kwargs["type"] not in {"apeglm", "ashr", "normal"}:
            raise ValueError("Invalid shrinkage type. Must be one of {'apeglm', 'ashr', 'normal'}")

        with localconverter(_converter):
            res = self.deseq_pkg.lfcShrink(dds=self.obj, **r_kwargs)
            return to_pandas(to_r_df(res))

    def get_results(self, **kwargs):
        """Extract results as a pandas DataFrame."""
        with localconverter(_converter):
            # Using standardized filtering via BaseWrapper
            res = self.deseq_pkg.results(self.obj, **kwargs)
            return to_pandas(to_r_df(res))

# --- Legacy Bridge ---
def run_deseq2(counts, metadata, design):
    """Bridge function for legacy code."""
    model = DESeq2(counts, metadata, design)
    return model.run_deseq().obj


