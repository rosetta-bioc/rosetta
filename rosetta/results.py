"""RosettaDataFrame — DataFrame subclass with .report() for result summaries."""

import pandas as pd


class RosettaDataFrame(pd.DataFrame):
    """A pandas DataFrame with a .report() method for summarizing statistical results."""

    _metadata = ["_rosetta_method"]

    @property
    def _constructor(self):
        return RosettaDataFrame

    def report(self, alpha: float = 0.05) -> str:
        """Print a human-readable summary of the results.

        Detects result type by column names and formats accordingly.
        """
        text = _build_report(self, alpha)
        print(text)
        return text


def _build_report(df: pd.DataFrame, alpha: float = 0.05) -> str:
    """Build report string based on detected result type."""
    if "padj" in df.columns and "log2FoldChange" in df.columns:
        return _report_deseq2(df, alpha)
    elif "FDR" in df.columns and "logFC" in df.columns:
        return _report_edger(df, alpha)
    elif "adj.P.Val" in df.columns and "logFC" in df.columns:
        return _report_limma(df, alpha)
    elif "p.adjust" in df.columns and "GeneRatio" in df.columns:
        return _report_enrichment(df, alpha)
    else:
        return f"Rosetta Results: {len(df)} rows × {len(df.columns)} columns"


def _report_deseq2(df: pd.DataFrame, alpha: float) -> str:
    n = len(df)
    sig = df["padj"].dropna() < alpha
    n_sig = sig.sum()
    lfc = df.loc[sig[sig].index, "log2FoldChange"] if n_sig > 0 else pd.Series(dtype=float)
    up = (lfc > 0).sum()
    down = (lfc < 0).sum()

    lines = [
        "DESeq2 Results Summary",
        "─" * 30,
        f"Total genes tested:      {n:,}",
        f"Significant (padj<{alpha}): {n_sig:,} ({100*n_sig/n:.1f}%)" if n > 0 else "Significant: 0",
        f"  ↑ Upregulated:         {up:,}",
        f"  ↓ Downregulated:       {down:,}",
    ]
    if n_sig > 0:
        lines.append(f"LFC range:               [{lfc.min():.2f}, {lfc.max():.2f}]")
    return "\n".join(lines)


def _report_edger(df: pd.DataFrame, alpha: float) -> str:
    n = len(df)
    sig = df["FDR"].dropna() < alpha
    n_sig = sig.sum()
    lfc = df.loc[sig[sig].index, "logFC"] if n_sig > 0 else pd.Series(dtype=float)
    up = (lfc > 0).sum()
    down = (lfc < 0).sum()

    lines = [
        "edgeR Results Summary",
        "─" * 30,
        f"Total genes tested:      {n:,}",
        f"Significant (FDR<{alpha}):  {n_sig:,} ({100*n_sig/n:.1f}%)" if n > 0 else "Significant: 0",
        f"  ↑ Upregulated:         {up:,}",
        f"  ↓ Downregulated:       {down:,}",
    ]
    if n_sig > 0:
        lines.append(f"logFC range:             [{lfc.min():.2f}, {lfc.max():.2f}]")
    return "\n".join(lines)


def _report_limma(df: pd.DataFrame, alpha: float) -> str:
    n = len(df)
    sig = df["adj.P.Val"].dropna() < alpha
    n_sig = sig.sum()
    lfc = df.loc[sig[sig].index, "logFC"] if n_sig > 0 else pd.Series(dtype=float)
    up = (lfc > 0).sum()
    down = (lfc < 0).sum()

    lines = [
        "limma Results Summary",
        "─" * 30,
        f"Total genes tested:      {n:,}",
        f"Significant (adj.P<{alpha}): {n_sig:,} ({100*n_sig/n:.1f}%)" if n > 0 else "Significant: 0",
        f"  ↑ Upregulated:         {up:,}",
        f"  ↓ Downregulated:       {down:,}",
    ]
    if n_sig > 0:
        lines.append(f"logFC range:             [{lfc.min():.2f}, {lfc.max():.2f}]")
    return "\n".join(lines)


def _report_enrichment(df: pd.DataFrame, alpha: float) -> str:
    n = len(df)
    sig = df["p.adjust"].dropna() < alpha
    n_sig = sig.sum()

    lines = [
        "Enrichment Results Summary",
        "─" * 30,
        f"Total terms tested:      {n:,}",
        f"Significant (p.adj<{alpha}): {n_sig:,}",
    ]
    if n_sig > 0 and "Description" in df.columns:
        top = df.loc[sig[sig].index].nsmallest(5, "p.adjust")
        lines.append("Top enriched terms:")
        for _, row in top.iterrows():
            desc = row["Description"][:50]
            lines.append(f"  • {desc} (p={row['p.adjust']:.2e})")
    return "\n".join(lines)
