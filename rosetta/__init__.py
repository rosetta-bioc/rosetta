"""rosetta — Seamless Python wrappers for R bioinformatics packages."""

from ._errors import RDataError, RFormulaError, RPackageMissing
from .wrappers.deseq2 import deseq2
from .wrappers.edger import edger
from .wrappers.limma import limma_voom
from .wrappers.clusterprofiler import enrich_go, enrich_kegg, enrich_pathway, enrich_custom
from .wrappers.phyloseq import phyloseq, phyloseq_richness
from .wrappers.seurat import seurat

# Alias for backward compatibility
enrichment = enrich_go

__all__ = [
    "deseq2", "edger", "limma_voom", "enrichment", 
    "enrich_go", "enrich_kegg", "enrich_pathway", "enrich_custom",
    "phyloseq", "phyloseq_richness", "seurat", 
    "RDataError", "RFormulaError", "RPackageMissing"
]
