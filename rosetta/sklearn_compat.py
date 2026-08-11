"""rosetta.sklearn_compat — sklearn TransformerMixin wrappers for Rosetta DE methods.

Makes DESeq2, EdgeR, and Limma first-class citizens in sklearn pipelines:

    from sklearn.pipeline import Pipeline
    from sklearn.decomposition import PCA
    from sklearn.ensemble import RandomForestClassifier
    from rosetta.sklearn_compat import DESeq2Transformer

    pipe = Pipeline([
        ('de',       DESeq2Transformer(counts, metadata)),
        ('pca',      PCA(n_components=50)),
        ('classify', RandomForestClassifier()),
    ])
    pipe.fit(X, y)

Each transformer:
  - fit(): runs the underlying DE analysis, stores result in self.result_
  - transform(): returns a numpy array suitable for downstream sklearn steps
    (samples × significant-genes log2FC matrix by default)
  - get_feature_names_out(): returns gene names for ColumnTransformer support

Requires: scikit-learn (pip install rosetta-bioc[sklearn])
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Optional

from .results import RosettaDataFrame


def _check_sklearn() -> None:
    try:
        import sklearn  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "scikit-learn is required for sklearn-compatible transformers: "
            "pip install 'rosetta-bioc[sklearn]'"
        ) from exc


class _BaseRosettaTransformer:
    """Mixin providing fit/transform/get_feature_names_out for Rosetta DE methods."""

    def __init__(
        self,
        counts: pd.DataFrame,
        metadata: pd.DataFrame,
        design: str = "~ condition",
        alpha: float = 0.05,
        lfc_threshold: float = 0.0,
        select: str = "significant",
        value: str = "log2FoldChange",
    ) -> None:
        """
        Args:
            counts: Gene count matrix (genes × samples), raw integers.
            metadata: Sample metadata, index matching count columns.
            design: R formula string.
            alpha: FDR threshold for gene selection.
            lfc_threshold: Minimum |log2FC| for gene selection.
            select: Which genes to include in transform output.
                "significant" — genes with padj/FDR < alpha (default)
                "all"         — all tested genes
            value: Column to use as feature values.
                "log2FoldChange" / "logFC" — fold change (default)
                "stat"                     — test statistic
        """
        _check_sklearn()
        self.counts = counts
        self.metadata = metadata
        self.design = design
        self.alpha = alpha
        self.lfc_threshold = lfc_threshold
        self.select = select
        self.value = value

        # Set by fit()
        self.result_: Optional[RosettaDataFrame] = None
        self.feature_names_: Optional[np.ndarray] = None

    # ------------------------------------------------------------------
    # sklearn interface
    # ------------------------------------------------------------------

    def fit(self, X=None, y=None):
        """Run the DE analysis. X and y are accepted but ignored (data is in __init__)."""
        self.result_ = self._run_de()
        self.result_._rosetta_method = self._method_name
        self.feature_names_ = self._select_genes()
        return self

    def transform(self, X=None):
        """Return a (samples × genes) numpy array of DE values.

        Rows = samples (ordered as metadata index).
        Cols = selected genes (ordered as feature_names_).

        Note: The DE result is gene-level, so each sample gets the same gene-level
        fold-change vector. To get a sample × gene expression matrix instead,
        pass value='counts' or subset counts directly.
        """
        if self.result_ is None:
            raise RuntimeError("Call fit() before transform().")

        result = self.result_
        genes = self.feature_names_

        # Build sample × gene matrix: replicate gene-level values across samples
        # This is the standard pattern for feeding DE features into a classifier
        # (each sample is represented by the fold-changes of selected genes)
        col = self._value_col()
        values = result.loc[genes, col].values.astype(float)

        n_samples = len(self.metadata)
        # shape: (n_samples, n_genes) — same vector tiled per sample
        return np.tile(values, (n_samples, 1))

    def fit_transform(self, X=None, y=None):
        return self.fit(X, y).transform(X)

    def get_feature_names_out(self, input_features=None):
        """Return gene names selected after fit()."""
        if self.feature_names_ is None:
            raise RuntimeError("Call fit() before get_feature_names_out().")
        return self.feature_names_

    # sklearn estimator API (needed for clone(), set_params(), Pipeline repr)
    def get_params(self, deep: bool = True) -> dict:
        return {
            "counts": self.counts,
            "metadata": self.metadata,
            "design": self.design,
            "alpha": self.alpha,
            "lfc_threshold": self.lfc_threshold,
            "select": self.select,
            "value": self.value,
        }

    def set_params(self, **params):
        for k, v in params.items():
            setattr(self, k, v)
        return self

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _select_genes(self) -> np.ndarray:
        result = self.result_
        padj_col = self._padj_col()
        lfc_col = self._value_col()

        mask = pd.Series([True] * len(result), index=result.index)
        if self.select == "significant":
            mask &= result[padj_col] < self.alpha
        if self.lfc_threshold > 0:
            mask &= result[lfc_col].abs() >= self.lfc_threshold

        genes = result.index[mask]
        if len(genes) == 0:
            raise ValueError(
                f"No genes passed filters (alpha={self.alpha}, "
                f"lfc_threshold={self.lfc_threshold}). "
                "Try increasing alpha or setting select='all'."
            )
        return np.array(genes)

    def _padj_col(self) -> str:
        for col in ("padj", "FDR", "adj.P.Val", "p.adjust"):
            if col in self.result_.columns:
                return col
        raise KeyError("Cannot find adjusted p-value column in results.")

    def _value_col(self) -> str:
        if self.value in self.result_.columns:
            return self.value
        # fall back to first available fold-change column
        for col in ("log2FoldChange", "logFC"):
            if col in self.result_.columns:
                return col
        raise KeyError(f"Column '{self.value}' not found in results.")

    # Subclasses implement these
    def _run_de(self) -> RosettaDataFrame:
        raise NotImplementedError

    @property
    def _method_name(self) -> str:
        raise NotImplementedError


class DESeq2Transformer(_BaseRosettaTransformer):
    """sklearn-compatible transformer backed by DESeq2.

    Example::

        from sklearn.pipeline import Pipeline
        from sklearn.decomposition import PCA
        from rosetta.sklearn_compat import DESeq2Transformer

        pipe = Pipeline([
            ('de',  DESeq2Transformer(counts, metadata, alpha=0.05)),
            ('pca', PCA(n_components=20)),
        ])
        transformed = pipe.fit_transform(None)
    """

    def __init__(self, counts, metadata, design="~ condition", alpha=0.05,
                 lfc_threshold=0.0, shrinkage=None, contrast=None,
                 select="significant", value="log2FoldChange"):
        super().__init__(counts, metadata, design, alpha, lfc_threshold, select, value)
        self.shrinkage = shrinkage
        self.contrast = contrast

    def get_params(self, deep=True):
        params = super().get_params(deep)
        params.update({"shrinkage": self.shrinkage, "contrast": self.contrast})
        return params

    def _run_de(self) -> RosettaDataFrame:
        from .pipelines import diff_expr
        return diff_expr(
            self.counts, self.metadata, self.design,
            method="deseq2", alpha=self.alpha, lfc_threshold=self.lfc_threshold,
            shrinkage=self.shrinkage, contrast=self.contrast,
        )

    @property
    def _method_name(self) -> str:
        return "deseq2"


class EdgeRTransformer(_BaseRosettaTransformer):
    """sklearn-compatible transformer backed by edgeR.

    Example::

        from rosetta.sklearn_compat import EdgeRTransformer
        tr = EdgeRTransformer(counts, metadata).fit()
        print(tr.result_.report())
    """

    def _run_de(self) -> RosettaDataFrame:
        from .pipelines import diff_expr
        return diff_expr(
            self.counts, self.metadata, self.design,
            method="edger", alpha=self.alpha, lfc_threshold=self.lfc_threshold,
        )

    @property
    def _method_name(self) -> str:
        return "edger"


class LimmaTransformer(_BaseRosettaTransformer):
    """sklearn-compatible transformer backed by limma-voom.

    Example::

        from rosetta.sklearn_compat import LimmaTransformer
        tr = LimmaTransformer(counts, metadata).fit()
        print(tr.result_.report())
    """

    def _run_de(self) -> RosettaDataFrame:
        from .pipelines import diff_expr
        return diff_expr(
            self.counts, self.metadata, self.design,
            method="limma", alpha=self.alpha, lfc_threshold=self.lfc_threshold,
        )

    @property
    def _method_name(self) -> str:
        return "limma"
