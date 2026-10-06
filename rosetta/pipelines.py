"""rosetta.pipelines — Complete analysis workflows in one call.

These are the "I just want results" functions. Each one runs the full
statistical pipeline and returns a RosettaDataFrame with .report().
"""
from typing import List, Optional

import pandas as pd

from .results import RosettaDataFrame


def diff_expr(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str = "~ condition",
    method: str = "deseq2",
    alpha: float = 0.05,
    lfc_threshold: float = 0.0,
    shrinkage: Optional[str] = None,
    contrast: Optional[list] = None,
) -> RosettaDataFrame:
    """Run differential expression — full pipeline, one call.

    Args:
        counts: Gene count matrix (genes × samples), raw integers.
        metadata: Sample metadata with index matching count columns.
        design: R formula string (e.g. "~ condition", "~ batch + treatment").
        method: One of "deseq2", "edger", "limma".
        alpha: FDR significance threshold.
        lfc_threshold: Minimum absolute log2 fold change.
        shrinkage: For DESeq2: "apeglm", "ashr", or "normal". None = no shrinkage.
        contrast: For DESeq2: [factor, numerator, denominator].

    Returns:
        RosettaDataFrame with .report() method.

    Example:
        >>> results = rb.pipelines.diff_expr(counts, meta, method="deseq2")
        >>> results.report()
        >>> sig_genes = results[results["padj"] < 0.05]
    """
    if method == "deseq2":
        from . import _bridge
        from .wrappers.deseq2 import DESeq2

        model = DESeq2(counts, metadata, design)
        model.run_deseq()

        if shrinkage:
            # Need coefficient name for shrinkage
            with _bridge.localconverter(_bridge._converter):
                coefs = list(model.deseq_pkg.resultsNames(model.r_obj))
            # Use last coefficient (typically the treatment effect)
            coef = coefs[-1] if coefs else None
            if coef:
                result = model.lfc_shrink(coef=coef, type=shrinkage)
                result._rosetta_method = "deseq2"
                return result

        result_kwargs = {"alpha": alpha, "lfcThreshold": lfc_threshold}
        if contrast:
            result_kwargs["contrast"] = _bridge.ro.StrVector(contrast)
        result = model.get_results(**result_kwargs)
        result._rosetta_method = "deseq2"
        return result

    elif method == "edger":
        from .wrappers.edger import EdgeR
        model = EdgeR(counts, metadata, design)
        test_result = model.run_test(lfc=lfc_threshold)
        result = model.get_results(test_result)
        result._rosetta_method = "edger"
        return result

    elif method == "limma":
        from .wrappers.limma import Limma
        model = Limma(counts, metadata, design)
        model.run_ebayes()
        result = model.get_results()
        result._rosetta_method = "limma"
        return result

    else:
        raise ValueError(f"Unknown method '{method}'. Use 'deseq2', 'edger', or 'limma'.")


def enrichment(
    gene_list: list[str],
    method: str = "go",
    organism: str = "hsa",
    org_db: str = "org.Hs.eg.db",
    ont: str = "BP",
    **kwargs,
) -> RosettaDataFrame:
    """Run pathway enrichment — full pipeline, one call.

    Args:
        gene_list: List of gene IDs (Entrez by default).
        method: One of "go", "kegg", "reactome".
        organism: KEGG organism code (default "hsa" for human).
        org_db: OrgDb for GO (default "org.Hs.eg.db").
        ont: GO ontology — "BP", "MF", or "CC".

    Returns:
        RosettaDataFrame with .report() method.

    Example:
        >>> results = rb.pipelines.enrichment(sig_genes, method="kegg")
        >>> results.report()
    """
    if method == "go":
        from .wrappers.clusterprofiler import enrich_go
        return enrich_go(gene_list, organism=org_db, ont=ont, **kwargs)
    elif method == "kegg":
        from .wrappers.clusterprofiler import enrich_kegg
        return enrich_kegg(gene_list, organism=organism, **kwargs)
    elif method == "reactome":
        from .wrappers.clusterprofiler import enrich_pathway
        return enrich_pathway(gene_list, **kwargs)
    else:
        raise ValueError(f"Unknown method '{method}'. Use 'go', 'kegg', or 'reactome'.")


def compare(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str = "~ condition",
    methods: Optional[List[str]] = None,
    alpha: float = 0.05,
) -> RosettaDataFrame:
    """Run multiple DE methods and return a comparison summary.

    This is the "which genes do all methods agree on?" function.

    Args:
        counts: Gene count matrix (genes × samples).
        metadata: Sample metadata.
        design: R formula string.
        methods: List of methods to compare. Default: ["deseq2", "edger", "limma"].
        alpha: FDR significance threshold.

    Returns:
        RosettaDataFrame with columns for each method's significance call
        and an 'n_methods' column showing agreement count.

    Example:
        >>> consensus = rb.pipelines.compare(counts, meta)
        >>> robust_genes = consensus[consensus["n_methods"] == 3]
    """
    if methods is None:
        methods = ["deseq2", "edger", "limma"]

    results = {}
    for method in methods:
        try:
            res = diff_expr(counts, metadata, design, method=method, alpha=alpha)
            # Extract significance column
            if "padj" in res.columns:
                results[method] = res["padj"] < alpha
            elif "FDR" in res.columns:
                results[method] = res["FDR"] < alpha
            elif "adj.P.Val" in res.columns:
                results[method] = res["adj.P.Val"] < alpha
        except Exception as e:
            print(f"  ⚠ {method} failed: {e}")
            continue

    if not results:
        raise RuntimeError("All methods failed")

    comparison = pd.DataFrame(results)
    comparison["n_methods"] = comparison.sum(axis=1)
    comparison = comparison.sort_values("n_methods", ascending=False)
    return RosettaDataFrame(comparison)
