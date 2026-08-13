"""QuickResult — Wrapper for Tier 1 Quick API outputs with .report() support."""

from __future__ import annotations

from typing import Any, Dict, Iterator, Optional

import pandas as pd

from .results import RosettaDataFrame


class QuickResult:
    """Wraps Quick API outputs (dict-like results) with a .report() method.

    Supports dict-like access (result['key'], result.keys(), 'key' in result)
    and dispatches .report() to method-specific reporters.

    Args:
        data: Dictionary of result data.
        method: Name of the Quick API method (e.g., 'seurat', 'phyloseq').
        metadata: Optional dictionary of additional metadata about the analysis.
    """

    def __init__(
        self,
        data: Dict[str, Any],
        method: str,
        metadata: Optional[Dict[str, Any]] = None,
    ):
        self._data = data
        self._method = method
        self._metadata = metadata or {}

    # --- Dict-like access ---

    def __getitem__(self, key: str) -> Any:
        return self._data[key]

    def __contains__(self, key: object) -> bool:
        return key in self._data

    def __iter__(self) -> Iterator[str]:
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)

    def keys(self):
        """Return the keys of the result data."""
        return self._data.keys()

    def values(self):
        """Return the values of the result data."""
        return self._data.values()

    def items(self):
        """Return the items of the result data."""
        return self._data.items()

    def get(self, key: str, default: Any = None) -> Any:
        """Get a value by key with optional default."""
        return self._data.get(key, default)

    # --- Properties ---

    @property
    def data(self) -> Dict[str, Any]:
        """The underlying result data dictionary."""
        return self._data

    @property
    def method(self) -> str:
        """The Quick API method name."""
        return self._method

    @property
    def metadata(self) -> Dict[str, Any]:
        """Additional metadata about the analysis."""
        return self._metadata

    # --- Report ---

    def report(self, **kwargs) -> str:
        """Print a human-readable summary of the results.

        Dispatches to a method-specific reporter based on self._method.
        Returns the report string.
        """
        reporters = {
            "seurat": _report_seurat,
            "phyloseq": _report_phyloseq,
            "deseq2": _report_deseq2,
            "edger": _report_edger,
        }
        reporter = reporters.get(self._method, _report_generic)
        text = reporter(self._data, self._metadata, **kwargs)
        print(text)
        return text

    def __repr__(self) -> str:
        keys = list(self._data.keys())
        return f"QuickResult(method={self._method!r}, keys={keys})"


# --- Method-specific reporters ---


def _report_seurat(data: Dict[str, Any], metadata: Dict[str, Any], **kwargs) -> str:
    """Report for Seurat quick results."""
    lines = [
        "Seurat Quick Analysis Summary",
        "─" * 30,
    ]

    clusters = data.get("clusters")
    if clusters is not None:
        if isinstance(clusters, pd.Series):
            n_cells = len(clusters)
            unique_clusters = clusters.nunique()
            lines.append(f"Total cells:             {n_cells:,}")
            lines.append(f"Clusters found:          {unique_clusters}")
            lines.append("Cluster sizes:")
            cluster_counts = clusters.value_counts().sort_index()
            for cluster_id, count in cluster_counts.items():
                lines.append(f"  Cluster {cluster_id}: {count:,} cells")
        else:
            lines.append(f"Clusters: {clusters}")

    umap = data.get("umap")
    if umap is not None and isinstance(umap, pd.DataFrame):
        lines.append(f"UMAP dimensions:         {umap.shape[1]}")

    var_features = data.get("variable_features")
    if var_features is not None:
        lines.append(f"Variable features:       {len(var_features):,}")

    # Include any metadata
    if metadata:
        for key, val in metadata.items():
            lines.append(f"{key}: {val}")

    return "\n".join(lines)


def _report_phyloseq(data: Dict[str, Any], metadata: Dict[str, Any], **kwargs) -> str:
    """Report for phyloseq quick results."""
    lines = [
        "Phyloseq Diversity Summary",
        "─" * 30,
    ]

    diversity = data.get("diversity")
    if diversity is not None and isinstance(diversity, pd.DataFrame):
        lines.append(f"Samples:                 {len(diversity):,}")
        lines.append("Diversity metrics:")
        for col in diversity.columns:
            values = diversity[col].dropna()
            if len(values) > 0:
                lines.append(
                    f"  {col}: mean={values.mean():.3f}, "
                    f"sd={values.std():.3f}, "
                    f"range=[{values.min():.3f}, {values.max():.3f}]"
                )
    elif isinstance(data, dict):
        # If the data dict itself contains DataFrames directly
        for key, val in data.items():
            if isinstance(val, pd.DataFrame):
                lines.append(f"Samples:                 {len(val):,}")
                lines.append("Diversity metrics:")
                for col in val.columns:
                    values = val[col].dropna()
                    if len(values) > 0:
                        lines.append(
                            f"  {col}: mean={values.mean():.3f}, "
                            f"sd={values.std():.3f}, "
                            f"range=[{values.min():.3f}, {values.max():.3f}]"
                        )
                break

    # Include any metadata
    if metadata:
        for key, val in metadata.items():
            lines.append(f"{key}: {val}")

    return "\n".join(lines)


def _report_deseq2(data: Dict[str, Any], metadata: Dict[str, Any], **kwargs) -> str:
    """Report for DESeq2 quick results — delegates to RosettaDataFrame.report."""
    results_df = data.get("results")
    if results_df is not None and isinstance(results_df, (pd.DataFrame, RosettaDataFrame)):
        if isinstance(results_df, RosettaDataFrame):
            return results_df.report(**kwargs)
        # Wrap in RosettaDataFrame for reporting
        rdf = RosettaDataFrame(results_df)
        return rdf.report(**kwargs)
    return "DESeq2 QuickResult: no results DataFrame found."


def _report_edger(data: Dict[str, Any], metadata: Dict[str, Any], **kwargs) -> str:
    """Report for edgeR quick results — delegates to RosettaDataFrame.report."""
    results_df = data.get("results")
    if results_df is not None and isinstance(results_df, (pd.DataFrame, RosettaDataFrame)):
        if isinstance(results_df, RosettaDataFrame):
            return results_df.report(**kwargs)
        rdf = RosettaDataFrame(results_df)
        return rdf.report(**kwargs)
    return "edgeR QuickResult: no results DataFrame found."


def _report_generic(data: Dict[str, Any], metadata: Dict[str, Any], **kwargs) -> str:
    """Fallback reporter for unknown methods."""
    lines = [
        "QuickResult Summary",
        "─" * 30,
        f"Keys: {list(data.keys())}",
    ]
    for key, val in data.items():
        if isinstance(val, pd.DataFrame):
            lines.append(f"  {key}: DataFrame ({val.shape[0]} rows × {val.shape[1]} cols)")
        elif isinstance(val, pd.Series):
            lines.append(f"  {key}: Series ({len(val)} items)")
        elif isinstance(val, list):
            lines.append(f"  {key}: list ({len(val)} items)")
        else:
            lines.append(f"  {key}: {type(val).__name__}")
    return "\n".join(lines)
