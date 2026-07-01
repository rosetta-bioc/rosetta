"""rosetta — Seamless Python wrappers for R bioinformatics packages."""

import importlib.metadata as _meta

# Detect namespace collision with unrelated 'rosetta' PyPI package
try:
    _dist = _meta.distribution("rosetta")
    if "bioinformatics" not in (_dist.metadata.get("Summary", "") + _dist.metadata.get("Keywords", "")).lower():
        import warnings
        warnings.warn(
            "Both 'rosetta' and 'rosetta-bioc' are installed. "
            "Run 'pip uninstall rosetta' to avoid import conflicts.",
            ImportWarning,
            stacklevel=2,
        )
except _meta.PackageNotFoundError:
    pass

from importlib.metadata import version as _get_version
try:
    __version__ = _get_version("rosetta-bioc")
except Exception:
    __version__ = "0.2.3.dev0"

from ._errors import RDataError, RFormulaError, RPackageMissing
from .results import RosettaDataFrame
from .wrappers.deseq2 import deseq2
from .wrappers.edger import edger
from .wrappers.limma import limma_voom
from .wrappers.clusterprofiler import ORA, GSEA
from .wrappers.phyloseq import Phyloseq
from .wrappers.seurat import Seurat
from . import pipelines
from . import codegen

# Top-level convenience aliases for ORA methods
enrich_go = ORA.enrich_go
enrich_kegg = ORA.enrich_kegg
enrich_pathway = ORA.enrich_pathway
enrich_custom = ORA.enrich_custom

# Alias for backward compatibility
enrichment = ORA

# --- Tier 1: Quick API (Added for Week 5) ---

def quick_deseq2(counts, metadata, design="~ condition", alpha=0.05, **kwargs):
    """
    Tier 1 API: Quick DESeq2 differential expression.
    Fits the model and extracts results in one call.

    Args:
        counts: Gene count matrix (genes x samples).
        metadata: Sample metadata DataFrame.
        design: R formula string for the experimental design.
        alpha: Significance cutoff (FDR).
        **kwargs: Additional arguments passed to get_results().

    Returns:
        RosettaDataFrame with baseMean, log2FoldChange, lfcSE, stat, pvalue, padj.
    """
    from .wrappers.deseq2 import run_deseq2, get_results
    dds = run_deseq2(counts, metadata, design)
    return get_results(dds, alpha=alpha, **kwargs)


def quick_edger(counts, metadata, design="~ condition", **kwargs):
    """
    Tier 1 API: Quick edgeR quasi-likelihood differential expression.
    Runs the full edgeR QL pipeline in one call.

    Args:
        counts: Gene count matrix (genes x samples).
        metadata: Sample metadata DataFrame.
        design: R formula string for the experimental design.
        **kwargs: Additional arguments passed to edger().

    Returns:
        RosettaDataFrame with logFC, logCPM, F, PValue, FDR.
    """
    from .wrappers.edger import edger as _edger
    return _edger(counts, metadata, design, **kwargs)


def quick_seurat(counts, **kwargs):
    """
    Tier 1 API: Quick Seurat analysis.
    Executes standard pipeline and returns results dictionary.
    """
    return Seurat(counts).run_standard_pipeline(**kwargs).get_results()

def quick_phyloseq(otu_table, sample_data=None, measures=["Shannon"], **kwargs):
    """
    Tier 1 API: Quick phyloseq analysis.
    Computes alpha diversity metrics.
    """
    ps = Phyloseq(otu_table, sample_data=sample_data, **kwargs)
    return ps.estimate_richness(measures=measures)

# Lowercase convenience aliases (backward compatibility with pre-class API)
seurat = quick_seurat
phyloseq = quick_phyloseq
phyloseq_richness = quick_phyloseq

# --- Exports ---

__all__ = [
    # Metadata
    "__version__",
    # Tier 3 (Functional/Legacy)
    "deseq2", "edger", "limma_voom",
    "ORA", "GSEA", "enrichment",
    "enrich_go", "enrich_kegg", "enrich_pathway", "enrich_custom",
    # Tier 2 (Class-based)
    "Seurat", "Phyloseq",
    # Tier 1 (Quick API)
    "quick_deseq2", "quick_edger", "quick_seurat", "quick_phyloseq",
    # Backward-compat aliases
    "phyloseq", "phyloseq_richness", "seurat",
    # Utilities
    "pipelines", "codegen",
    "RosettaDataFrame",
    # Errors
    "RDataError", "RFormulaError", "RPackageMissing",
]
