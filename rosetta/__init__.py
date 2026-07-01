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

def quick_seurat(counts, **kwargs):
    """
    Tier 1 API: Quick Seurat analysis.
    Executes standard pipeline and returns results dictionary.
    """
    return Seurat(counts).run_standard_pipeline(**kwargs).get_results()

def quick_phyloseq(otu_table, sample_data=None, measures=["Shannon"], **kwargs):
    """
    Tier 1 API: Quick Phyloseq analysis.
    Computes alpha diversity metrics.
    """
    ps = Phyloseq(otu_table, sample_data=sample_data, **kwargs)
    return ps.estimate_richness(measures=measures)

# --- Exports ---

__all__ = [
    # Tier 3 (Functional/Legacy)
    "deseq2", "edger", "limma_voom",
    "ORA", "GSEA", "enrichment",
    "enrich_go", "enrich_kegg", "enrich_pathway", "enrich_custom",
    # Tier 2 (Class-based)
    "Seurat", "Phyloseq",
    # Tier 1 (Quick API)
    "quick_seurat", "quick_phyloseq",
    # Utilities
    "pipelines", "codegen",
    "RosettaDataFrame",
    # Errors
    "RDataError", "RFormulaError", "RPackageMissing",
]
