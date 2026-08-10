"""scikit-learn compatible pipeline transformers for Rosetta."""

import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted
from ._bridge import to_r_matrix, to_pandas, to_r_df, _converter
from ._deps import ensure_installed
from ._errors import RDataError
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr


class DESeq2Transformer(BaseEstimator, TransformerMixin):
    """scikit-learn compatible Transformer wrapper for DESeq2 differential expression analysis.

    NOTE: ``metadata`` is a constructor parameter rather than part of ``X`` so that it flows
    naturally through a Pipeline.  This means ``clone()`` / ``GridSearchCV`` will pass the
    same DataFrame object to the cloned estimator, which is fine for read-only use but means
    you cannot grid-search over metadata itself.
    """

    def __init__(
        self,
        metadata: pd.DataFrame,
        design: str = "~ condition",
        alpha: float = 0.05,
        lfc_threshold: float = 0.0,
    ):
        self.metadata = metadata
        self.design = design
        self.alpha = alpha
        self.lfc_threshold = lfc_threshold

    def fit(self, X: pd.DataFrame, y=None):
        """Fit the DESeq2 model using gene count matrix X (samples x genes)."""
        ensure_installed("DESeq2")
        deseq2_pkg = importr("DESeq2")
        stats_pkg = importr("stats")

        if not isinstance(X, pd.DataFrame):
            raise TypeError(
                "Input X must be a pandas DataFrame with gene names as columns and sample names as rows."
            )
        if X.columns.empty or not any(isinstance(c, str) for c in X.columns):
            raise ValueError("Input DataFrame X must have valid gene names as column headers.")
        if (X < 0).any().any():
            raise RDataError("Count matrix contains negative values")

        # Transpose from sklearn's (samples x genes) to R's (genes x samples)
        X_r_input = X.T

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
                raise RDataError(f"Failed to fit DESeq2 model: {e}") from e

        res = deseq2_pkg.results(self.dds_, alpha=self.alpha, lfcThreshold=self.lfc_threshold)
        res_df = to_pandas(to_r_df(res))
        self.significant_genes_ = res_df[res_df["padj"] < self.alpha].index.tolist()

        return self

    def transform(self, X: pd.DataFrame) -> pd.DataFrame:
        """Filter X to the significant genes identified during fit()."""
        check_is_fitted(self)

        valid_genes = [g for g in self.significant_genes_ if g in X.columns]
        if not valid_genes:
            raise ValueError("None of the significant genes found by DESeq2 exist in the input matrix columns.")

        return X[valid_genes]


class EdgeRTransformer(BaseEstimator, TransformerMixin):
    """scikit-learn compatible transformer wrapper for edgeR differential expression analysis.

    NOTE: ``metadata`` is a constructor parameter rather than part of ``X`` so that it flows
    naturally through a Pipeline.  This means ``clone()`` / ``GridSearchCV`` will pass the
    same DataFrame object to the cloned estimator, which is fine for read-only use but means
    you cannot grid-search over metadata itself.
    """

    def __init__(
        self,
        metadata: pd.DataFrame,
        design: str = "~ condition",
        alpha: float = 0.05,
        filter_low_counts: bool = True,
    ):
        self.metadata = metadata
        self.design = design
        self.alpha = alpha
        self.filter_low_counts = filter_low_counts

    def fit(self, X: pd.DataFrame, y=None):
        """Fit the edgeR model using gene count matrix X (samples x genes)."""
        ensure_installed("edgeR")
        import rpy2.robjects as ro

        edgeR_pkg = importr("edgeR")
        stats_pkg = importr("stats")

        if not isinstance(X, pd.DataFrame):
            raise TypeError(
                "Input X must be a pandas DataFrame with gene names as columns and sample names as rows."
            )
        if X.columns.empty or not any(isinstance(c, str) for c in X.columns):
            raise ValueError("Input DataFrame X must have valid gene names as column headers.")
        if (X < 0).any().any():
            raise RDataError("Count matrix contains negative values")

        # Transpose from sklearn's (samples x genes) to R's (genes x samples)
        X_r_input = X.T

        r_counts = to_r_matrix(X_r_input)
        r_metadata = to_r_df(self.metadata)

        try:
            with localconverter(_converter):
                dge = edgeR_pkg.DGEList(counts=r_counts)
                r_design = stats_pkg.model_matrix(
                    stats_pkg.as_formula(self.design), data=r_metadata
                )

                if self.filter_low_counts:
                    keep = edgeR_pkg.filterByExpr(dge, design=r_design)
                    dge = dge.rx(keep, True)

                dge = edgeR_pkg.calcNormFactors(dge)
                dge = edgeR_pkg.estimateDisp(dge, r_design)
                fit = edgeR_pkg.glmQLFit(dge, r_design)
                qlf = edgeR_pkg.glmQLFTest(fit)

                total_genes = ro.r["nrow"](qlf)[0]
                res = edgeR_pkg.topTags(qlf, n=total_genes)
                res_df = to_pandas(to_r_df(res))
                self.dge_ = dge
        except RDataError:
            raise
        except Exception as e:
            raise RDataError(f"Failed to fit edgeR model: {e}") from e

        # edgeR reports FDR when using topTags; fall back to PValue if absent
        pval_col = "FDR" if "FDR" in res_df.columns else "PValue"
        self.significant_genes_ = res_df[res_df[pval_col] < self.alpha].index.tolist()

        return self

    def transform(self, X: pd.DataFrame) -> pd.DataFrame:
        """Filter X to the significant genes identified during fit()."""
        check_is_fitted(self)

        if not isinstance(X, pd.DataFrame):
            X = pd.DataFrame(X)

        valid_genes = [g for g in self.significant_genes_ if g in X.columns]
        if not valid_genes:
            raise ValueError("None of the significant genes found by edgeR exist in the input matrix columns.")

        return X[valid_genes]
