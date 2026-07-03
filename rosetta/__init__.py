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
from .quick_result import QuickResult
from .wrappers.deseq2 import deseq2, run_deseq2, get_results, lfc_shrink
from .wrappers.edger import edger
from .wrappers.limma import limma_voom
from .wrappers.normalize import vst, rlog, tmm_normalize
from .wrappers.clusterprofiler import ORA, GSEA
from .wrappers.phyloseq import Phyloseq
from .wrappers.seurat import Seurat
from .wrappers.vcf import VCF
from . import pipelines
from . import codegen
from . import plots

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
    Executes standard pipeline and returns QuickResult with dict-like access and .report().
    """
    data = Seurat(counts).run_standard_pipeline(**kwargs).get_results()
    metadata = {"pipeline": "standard", **kwargs}
    return QuickResult(data=data, method="seurat", metadata=metadata)

def quick_phyloseq(otu_table, sample_data=None, measures=["Shannon"], **kwargs):
    """
    Tier 1 API: Quick phyloseq analysis.
    Computes alpha diversity metrics and returns QuickResult with .report().
    """
    ps = Phyloseq(otu_table, sample_data=sample_data, **kwargs)
    diversity_df = ps.estimate_richness(measures=measures)
    data = {"diversity": diversity_df, "measures": measures}
    metadata = {"measures": measures}
    return QuickResult(data=data, method="phyloseq", metadata=metadata)

def quick_locate_variants(filepath, genome="hg38", txdb=None, region="all"):
    """
    Tier 1 API: Read VCF and annotate variant locations in one call.

    Args:
        filepath: Path to VCF/VCF.gz/BCF file.
        genome: Genome assembly name (e.g., "hg38", "hg19").
        txdb: TxDb package name. Defaults to TxDb.Hsapiens.UCSC.<genome>.knownGene.
        region: One of 'all', 'coding', 'intron', 'fiveUTR', 'threeUTR',
                'intergenic', 'spliceSite', 'promoter'.

    Returns:
        DataFrame with CHROM, POS, LOCATION, QUERYID, TXID, GENEID.

    Raises:
        RDataError: If the file does not exist.
    """
    from .wrappers.variant_annotation import read_vcf, locate_variants

    if txdb is None:
        txdb = f"TxDb.Hsapiens.UCSC.{genome}.knownGene"
    vcf = read_vcf(filepath, genome)
    return locate_variants(vcf, txdb=txdb, region=region)


def quick_predict_coding(filepath, genome="hg38", txdb=None, bsgenome=None):
    """
    Tier 1 API: Read VCF and predict coding consequences in one call.

    Args:
        filepath: Path to VCF/VCF.gz/BCF file.
        genome: Genome assembly name (e.g., "hg38", "hg19").
        txdb: TxDb package name. Defaults to TxDb.Hsapiens.UCSC.<genome>.knownGene.
        bsgenome: BSgenome package name. Defaults to BSgenome.Hsapiens.UCSC.<genome>.

    Returns:
        DataFrame with CONSEQUENCE, REFCODON, VARCODON, REFAA, VARAA, GENEID.

    Raises:
        RDataError: If the file does not exist.
    """
    from .wrappers.variant_annotation import read_vcf, predict_coding

    if txdb is None:
        txdb = f"TxDb.Hsapiens.UCSC.{genome}.knownGene"
    if bsgenome is None:
        bsgenome = f"BSgenome.Hsapiens.UCSC.{genome}"
    vcf = read_vcf(filepath, genome)
    return predict_coding(vcf, txdb=txdb, bsgenome=bsgenome)


# Lowercase convenience aliases (backward compatibility with pre-class API)
seurat = quick_seurat
phyloseq = quick_phyloseq
phyloseq_richness = quick_phyloseq

# --- Exports ---

__all__ = [
    # Metadata
    "__version__",
    # Tier 3 (Functional/Legacy)
    "deseq2", "run_deseq2", "get_results", "lfc_shrink",
    "edger", "limma_voom",
    "vst", "rlog", "tmm_normalize",
    "ORA", "GSEA", "enrichment",
    "enrich_go", "enrich_kegg", "enrich_pathway", "enrich_custom",
    # Tier 2 (Class-based)
    "Seurat", "Phyloseq", "VCF",
    # Tier 1 (Quick API)
    "quick_deseq2", "quick_edger", "quick_seurat", "quick_phyloseq",
    "quick_locate_variants", "quick_predict_coding",
    # Backward-compat aliases
    "phyloseq", "phyloseq_richness", "seurat",
    # Utilities
    "pipelines", "codegen", "plots",
    "RosettaDataFrame",
    "QuickResult",
    # Errors
    "RDataError", "RFormulaError", "RPackageMissing",
]
