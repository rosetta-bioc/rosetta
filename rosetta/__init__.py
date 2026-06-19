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

# Top-level convenience aliases for ORA methods
enrich_go = ORA.enrich_go
enrich_kegg = ORA.enrich_kegg
enrich_pathway = ORA.enrich_pathway
enrich_custom = ORA.enrich_custom
from .wrappers.phyloseq import phyloseq, phyloseq_richness
from .wrappers.seurat import seurat
from . import pipelines
from . import codegen

# Alias for backward compatibility
enrichment = ORA

__all__ = [
    "deseq2", "edger", "limma_voom",
    "ORA", "GSEA", "enrichment",
    "enrich_go", "enrich_kegg", "enrich_pathway", "enrich_custom",
    "phyloseq", "phyloseq_richness", "seurat",
    "pipelines", "codegen",
    "RosettaDataFrame",
    "RDataError", "RFormulaError", "RPackageMissing",
]
