"""scikit-learn compatible pipeline transformers for Rosetta."""

import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from ._bridge import to_r_matrix, to_pandas, to_r_df, _converter
from ._deps import ensure_installed
from ._errors import RDataError, RFormulaError
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr


class DESeq2Transformer(BaseEstimator, TransformerMixin):
    """scikit-learn compatible Transformer wrapper for DESeq2 differential expression analysis."""

    def __init__(self, metadata: pd.DataFrame, design: str = "~ condition", alpha: float = 0.05, lfc_threshold: float = 0.0):
        self.metadata = metadata
        self.design = design
        self.alpha = alpha
        self.lfc_threshold = lfc_threshold
        
        # Fitted attributes
        self.dds_ = None
        self.significant_genes_ = None

    def fit(self, X: pd.DataFrame, y=None):
        """Fit the DESeq2 model using gene count matrix X (samples x genes)."""
        ensure_installed("DESeq2")
        deseq2_pkg = importr("DESeq2")
        stats_pkg = importr("stats")

        # Ensure X is a standard pandas DataFrame, make sure gene and sample are not lost
        if not isinstance(X, pd.DataFrame):
            raise TypeError(
                "Input X must be a pandas DataFrame with gene names as columns and sample names as rows."
            )
        # ensure that the DataFrame has valid gene names as columns
        if X.columns.empty or not any(isinstance(c, str) for c in X.columns):
            raise ValueError("Input DataFrame X must have valid gene names as column headers.")

        # Basic validations
        if (X < 0).any().any():
            raise RDataError("Count matrix contains negative values")
        
        # 1. Automatically transpose from sklearn's (samples x genes) to R's (genes x samples)
        X_r_input = X.T 

        # Prepare R data structures and design formula
        r_design = stats_pkg.as_formula(self.design)
        r_counts = to_r_matrix(X_r_input)
        r_metadata = to_r_df(self.metadata)

        with localconverter(_converter):
            try:
                dds = deseq2_pkg.DESeqDataSetFromMatrix(
                    countData=r_counts, colData=r_metadata, design=r_design
                )
                self.dds_ = deseq2_pkg.DESeq(dds)
            except Exception as e:
                raise RDataError(f"Failed to fit DESeq2 model: {e}")

        res = deseq2_pkg.results(self.dds_, alpha=self.alpha, lfcThreshold=self.lfc_threshold)
        res_df = to_pandas(to_r_df(res))
        
        self.significant_genes_ = res_df[res_df['padj'] < self.alpha].index.tolist()
        
        return self

    def transform(self, X: pd.DataFrame) -> pd.DataFrame:
        """Transform X by filtering it down to only the significant genes found during fit()."""
        if self.significant_genes_ is None:
            raise RuntimeError("DESeq2Transformer is not fitted yet. Call 'fit' before 'transform'.")
        
        # Ensure input columns contain the identified significant genes
        valid_genes = [g for g in self.significant_genes_ if g in X.columns]
        
        if not valid_genes:
            raise ValueError("None of the significant genes found by DESeq2 exist in the input matrix columns.")
            
        return X[valid_genes]
