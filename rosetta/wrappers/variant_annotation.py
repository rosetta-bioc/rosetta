"""VariantAnnotation — VCF reading, variant location, and coding prediction.

Tier 3 functional API: maps 1:1 to R/Bioconductor VariantAnnotation functions.
All functions accept Python types and return pandas DataFrames (or R objects
for intermediate pipeline steps).
"""

import os
import pandas as pd

from .._bridge import ACTIVE_BACKEND, _converter, to_pandas, to_r_df
from .._deps import ensure_installed
from .._errors import RDataError
from .. import codegen

# Conditionally import rpy2 components based on the active backend
if ACTIVE_BACKEND == "rpy2":
    import rpy2.robjects as ro
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr
else:
    ro = None
    localconverter = None
    importr = None


# --- Tier 3: Functional API ---


def read_vcf(filepath: str, genome: str = "", param: dict = None) -> "ro.RS4":
    """
    Read a VCF/BCF file into an R CollapsedVCF object.

    Args:
        filepath: Path to VCF/VCF.gz/BCF file.
        genome: Genome assembly name (e.g., "hg38", "hg19").
        param: Dict with keys 'info', 'geno', 'samples', 'which' to filter.

    Returns:
        R CollapsedVCF object (use vcf_to_dataframe() to convert).

    Raises:
        RDataError: If the file does not exist.
    """
    filepath = str(filepath)
    if not os.path.isfile(filepath):
        raise RDataError(f"VCF file not found: {filepath}")

    ensure_installed("VariantAnnotation")
    va = importr("VariantAnnotation")

    with localconverter(_converter):
        if param:
            r_param = _build_scan_param(param)
            codegen._emit(f'vcf <- readVcf("{filepath}", "{genome}", param=param)')
            return va.readVcf(filepath, genome, param=r_param)
        else:
            codegen._emit(f'vcf <- readVcf("{filepath}", "{genome}")')
            return va.readVcf(filepath, genome)


def vcf_to_dataframe(vcf) -> pd.DataFrame:
    """
    Convert a VCF object to a pandas DataFrame.

    Merges rowRanges (CHROM, POS, REF, ALT, QUAL, FILTER) with INFO columns.

    Args:
        vcf: R CollapsedVCF object (from read_vcf()).

    Returns:
        DataFrame with fixed fields and INFO columns.
    """
    ensure_installed("VariantAnnotation")

    with localconverter(_converter):
        # Use R-side code to produce a clean data.frame with character columns
        # that rpy2 can convert without hitting complex Bioconductor types.
        codegen._emit("df <- as.data.frame(rowRanges(vcf))")
        codegen._emit("info_df <- as.data.frame(info(vcf))")

        r_code = """
        function(vcf) {
            rr <- SummarizedExperiment::rowRanges(vcf)
            df <- data.frame(
                CHROM = as.character(GenomicRanges::seqnames(rr)),
                POS = BiocGenerics::start(rr),
                REF = as.character(VariantAnnotation::ref(vcf)),
                QUAL = VariantAnnotation::qual(vcf),
                FILTER = VariantAnnotation::filt(vcf),
                stringsAsFactors = FALSE
            )
            # INFO columns
            info_df <- as.data.frame(VariantAnnotation::info(vcf))
            # Convert list-columns to character for safe transfer
            for (col in names(info_df)) {
                if (is.list(info_df[[col]])) {
                    info_df[[col]] <- vapply(info_df[[col]],
                        function(x) paste(x, collapse=","), character(1))
                }
            }
            if (ncol(info_df) > 0) {
                df <- cbind(df, info_df)
            }
            df
        }
        """
        extract_fn = ro.r(r_code)
        r_df = extract_fn(vcf)
        return to_pandas(r_df)


def scan_vcf_header(filepath: str) -> dict:
    """
    Read VCF header metadata without loading variant data.

    Args:
        filepath: Path to VCF/VCF.gz/BCF file.

    Returns:
        Dict with keys 'info', 'geno', 'samples'.

    Raises:
        RDataError: If the file does not exist.
    """
    filepath = str(filepath)
    if not os.path.isfile(filepath):
        raise RDataError(f"VCF file not found: {filepath}")

    ensure_installed("VariantAnnotation")
    va = importr("VariantAnnotation")
    base = importr("base")

    with localconverter(_converter):
        codegen._emit(f'hdr <- scanVcfHeader("{filepath}")')
        hdr = va.scanVcfHeader(filepath)

        codegen._emit("info(hdr)")
        r_info = ro.r("VariantAnnotation::info")
        info_df = to_pandas(to_r_df(r_info(hdr)))

        codegen._emit("geno(hdr)")
        geno_df = to_pandas(to_r_df(base.as_data_frame(va.geno(hdr))))

        codegen._emit("samples(hdr)")
        r_samples = ro.r("VariantAnnotation::samples")
        samples = list(r_samples(hdr))

        return {
            "info": info_df,
            "geno": geno_df,
            "samples": samples,
        }


def locate_variants(vcf, txdb=None, region: str = "all") -> pd.DataFrame:
    """
    Annotate variant genomic locations (coding, intronic, UTR, etc.).

    Args:
        vcf: R VCF object (from read_vcf()).
        txdb: TxDb package name (str) or R TxDb object.
        region: One of 'all', 'coding', 'intron', 'fiveUTR', 'threeUTR',
                'intergenic', 'spliceSite', 'promoter'.

    Returns:
        DataFrame with CHROM, POS, LOCATION, QUERYID, TXID, GENEID.
    """
    ensure_installed("VariantAnnotation")
    va = importr("VariantAnnotation")
    base = importr("base")

    txdb_obj = _resolve_txdb(txdb)
    region_obj = _make_region(region)

    with localconverter(_converter):
        codegen._emit(f"loc <- locateVariants(vcf, txdb, {region}())")
        loc = va.locateVariants(vcf, txdb_obj, region_obj)
        return to_pandas(to_r_df(base.as_data_frame(loc)))


def predict_coding(vcf, txdb=None, bsgenome=None) -> pd.DataFrame:
    """
    Predict amino acid coding changes for variants in coding regions.

    Args:
        vcf: R VCF object (from read_vcf()).
        txdb: TxDb package name (str) or R TxDb object.
        bsgenome: BSgenome package name (str) or R BSgenome object.

    Returns:
        DataFrame with CONSEQUENCE, REFCODON, VARCODON, REFAA, VARAA, GENEID.
    """
    ensure_installed("VariantAnnotation")
    va = importr("VariantAnnotation")
    base = importr("base")

    txdb_obj = _resolve_txdb(txdb)
    bsg_obj = _resolve_bsgenome(bsgenome)

    with localconverter(_converter):
        codegen._emit("coding <- predictCoding(vcf, txdb, bsgenome)")
        coding = va.predictCoding(vcf, txdb_obj, bsg_obj)
        return to_pandas(to_r_df(base.as_data_frame(coding)))


def write_vcf(vcf, filepath: str, index: bool = False) -> str:
    """
    Write a VCF object to file.

    Args:
        vcf: R VCF object (from read_vcf()).
        filepath: Output file path.
        index: Whether to create a tabix index.

    Returns:
        The output filepath string.
    """
    ensure_installed("VariantAnnotation")
    va = importr("VariantAnnotation")

    with localconverter(_converter):
        codegen._emit(
            f'writeVcf(vcf, "{filepath}", index={str(index).upper()})'
        )
        va.writeVcf(vcf, filepath, index=index)

    return filepath


# --- Private Helpers ---


def _build_scan_param(param: dict):
    """Convert Python param dict to R ScanVcfParam object.

    Args:
        param: Dict with optional keys 'info', 'geno', 'samples', 'which'.
               'which' should map chromosome names to (start, end) tuples.

    Returns:
        R ScanVcfParam object.
    """
    ensure_installed("VariantAnnotation")
    va = importr("VariantAnnotation")

    kwargs = {}

    if "info" in param:
        kwargs["info"] = (
            ro.StrVector(param["info"]) if param["info"] else ro.NA_Character
        )
    if "geno" in param:
        kwargs["geno"] = (
            ro.StrVector(param["geno"]) if param["geno"] else ro.NA_Character
        )
    if "samples" in param:
        kwargs["samples"] = ro.StrVector(param["samples"])
    if "which" in param:
        ensure_installed("GenomicRanges")
        ensure_installed("IRanges")
        gr_pkg = importr("GenomicRanges")
        ir_pkg = importr("IRanges")

        regions = param["which"]
        seqnames = []
        starts = []
        ends = []
        for chrom, (s, e) in regions.items():
            seqnames.append(chrom)
            starts.append(int(s))
            ends.append(int(e))

        kwargs["which"] = gr_pkg.GRanges(
            seqnames=ro.StrVector(seqnames),
            ranges=ir_pkg.IRanges(
                start=ro.IntVector(starts), end=ro.IntVector(ends)
            ),
        )

    return va.ScanVcfParam(**kwargs)


def _resolve_txdb(txdb):
    """Resolve a TxDb from package name string or pass through R object.

    For string inputs, the R package name IS the object name — load
    the package and retrieve the object of the same name.
    """
    if txdb is None:
        raise RDataError("txdb argument is required (package name or R object)")
    if isinstance(txdb, str):
        ensure_installed(txdb)
        codegen._emit(f"library({txdb})")
        return ro.r(f"library({txdb}); get(\"{txdb}\")")
    return txdb


def _resolve_bsgenome(bsgenome):
    """Resolve a BSgenome from package name string or pass through R object.

    For string inputs, the R package name IS the object name — load
    the package and retrieve the object of the same name.
    """
    if bsgenome is None:
        raise RDataError("bsgenome argument is required (package name or R object)")
    if isinstance(bsgenome, str):
        ensure_installed(bsgenome)
        codegen._emit(f"library({bsgenome})")
        return ro.r(f"library({bsgenome}); get(\"{bsgenome}\")")
    return bsgenome


def _make_region(region: str):
    """Convert region string to an R VariantType object.

    Args:
        region: One of 'all', 'coding', 'intron', 'fiveUTR', 'threeUTR',
                'intergenic', 'spliceSite', 'promoter'.

    Returns:
        Instantiated R VariantType object (e.g., AllVariants()).

    Raises:
        ValueError: If region is not a recognized type.
    """
    _REGION_MAP = {
        "all": "AllVariants",
        "coding": "CodingVariants",
        "intron": "IntronVariants",
        "fiveUTR": "FiveUTRVariants",
        "threeUTR": "ThreeUTRVariants",
        "intergenic": "IntergenicVariants",
        "spliceSite": "SpliceSiteVariants",
        "promoter": "PromoterVariants",
    }

    if region not in _REGION_MAP:
        raise ValueError(
            f"Unknown region '{region}'. Must be one of: {list(_REGION_MAP.keys())}"
        )

    # Import VariantAnnotation and call the constructor
    ensure_installed("VariantAnnotation")
    va = importr("VariantAnnotation")
    constructor = getattr(va, _REGION_MAP[region])
    return constructor()
