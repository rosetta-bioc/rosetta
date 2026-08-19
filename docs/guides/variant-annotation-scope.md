# Scope: `rosetta.variant_annotation` Module

## Vision

```python
import rosetta as rb

# Tier 1 — one-liner variant annotation
results = rb.quick_variant_annotation("sample.vcf.gz", genome="hg38")
results.report()
```
```
Variant Annotation Summary
──────────────────────────────
Total variants:           1,247
  Coding:                   312 (25.0%)
  Intronic:                 589 (47.2%)
  Intergenic:               234 (18.8%)
  UTR (5'/3'):              112 (9.0%)
Consequences:
  Synonymous:               187
  Missense:                  98
  Nonsense:                  15
  Frameshift:                12
```

## API Design (Three Tiers)

### Tier 1 — Quick API

```python
# Variant location annotation
locations = rb.quick_locate_variants("sample.vcf.gz", genome="hg38")
# → RosettaDataFrame with CHROM, POS, REF, ALT, LOCATION, GENEID, TXID

# Coding consequence prediction
coding = rb.quick_predict_coding("sample.vcf.gz", genome="hg38")
# → RosettaDataFrame with CHROM, POS, REF, ALT, CONSEQUENCE, REFCODON, VARCODON, REFAA, VARAA

# Read VCF to DataFrame (no annotation, just parse)
variants = rb.read_vcf("sample.vcf.gz")
# → RosettaDataFrame with CHROM, POS, ID, REF, ALT, QUAL, FILTER, + INFO columns
```

### Tier 2 — Class-based (stateful, chainable)

```python
vcf = rb.VCF("sample.vcf.gz", genome="hg38")

# Chain operations
results = (vcf
    .filter(QUAL > 30, FILTER == "PASS")
    .locate_variants()
    .predict_coding()
    .to_dataframe()
)

# Access sub-results
vcf.info()        # → DataFrame of INFO fields
vcf.geno()        # → dict of genotype DataFrames
vcf.header()      # → dict of header metadata
```

### Tier 3 — Functional (full control, maps 1:1 to R functions)

```python
from rosetta.wrappers.variant_annotation import (
    read_vcf,
    scan_vcf_header,
    locate_variants,
    predict_coding,
    write_vcf,
)

# Read with params
vcf = read_vcf("sample.vcf.gz", genome="hg38",
               param={"info": ["AF", "DP"], "geno": ["GT", "AD"],
                      "which": {"chr22": (16000000, 17000000)}})

# Annotate
locations = locate_variants(vcf, txdb="TxDb.Hsapiens.UCSC.hg38.knownGene",
                           region="all")

# Predict coding changes
coding = predict_coding(vcf, txdb="TxDb.Hsapiens.UCSC.hg38.knownGene",
                       bsgenome="BSgenome.Hsapiens.UCSC.hg38")

# Write back
write_vcf(vcf, "output.vcf", index=True)
```

## Module Structure

```
rosetta/wrappers/variant_annotation.py    # Tier 3 functional API
rosetta/wrappers/__init__.py              # register exports
rosetta/__init__.py                       # add Tier 1 quick_* and Tier 2 VCF class
tests/test_variant_annotation.py          # unit tests
```

## Implementation: `rosetta/wrappers/variant_annotation.py`

```python
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
    base = importr("base")

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
```

## Tier 1 Quick API (in `__init__.py`)

```python
def quick_locate_variants(filepath, genome="hg38", txdb=None, region="all", **kwargs):
    """Tier 1 API: Read VCF + annotate locations in one call."""
    from .wrappers.variant_annotation import read_vcf, locate_variants
    if txdb is None:
        txdb = f"TxDb.Hsapiens.UCSC.{genome}.knownGene"
    vcf = read_vcf(filepath, genome)
    return locate_variants(vcf, txdb=txdb, region=region, **kwargs)

def quick_predict_coding(filepath, genome="hg38", txdb=None, bsgenome=None, **kwargs):
    """Tier 1 API: Read VCF + predict coding consequences in one call."""
    from .wrappers.variant_annotation import read_vcf, predict_coding
    if txdb is None:
        txdb = f"TxDb.Hsapiens.UCSC.{genome}.knownGene"
    if bsgenome is None:
        bsgenome = f"BSgenome.Hsapiens.UCSC.{genome}"
    vcf = read_vcf(filepath, genome)
    return predict_coding(vcf, txdb=txdb, bsgenome=bsgenome, **kwargs)

def read_vcf(filepath, genome="", **kwargs):
    """Tier 1 API: Read VCF to DataFrame."""
    from .wrappers.variant_annotation import read_vcf as _read_vcf, vcf_to_dataframe
    vcf = _read_vcf(filepath, genome, **kwargs)
    return vcf_to_dataframe(vcf)
```

## Tier 2 Class API

```python
class VCF:
    """Stateful VCF analysis object with chainable methods."""

    def __init__(self, filepath: str, genome: str = "hg38", param: dict = None):
        from .wrappers.variant_annotation import read_vcf
        self._vcf = read_vcf(filepath, genome, param=param)
        self._genome = genome
        self._locations = None
        self._coding = None

    def filter(self, qual_min=None, filter_pass=True, regions=None):
        """Filter variants (subset the R object)."""
        # ... subset logic ...
        return self

    def locate_variants(self, txdb=None, region="all"):
        """Annotate variant locations."""
        from .wrappers.variant_annotation import locate_variants
        if txdb is None:
            txdb = f"TxDb.Hsapiens.UCSC.{self._genome}.knownGene"
        self._locations = locate_variants(self._vcf, txdb=txdb, region=region)
        return self

    def predict_coding(self, txdb=None, bsgenome=None):
        """Predict coding consequences."""
        from .wrappers.variant_annotation import predict_coding
        if txdb is None:
            txdb = f"TxDb.Hsapiens.UCSC.{self._genome}.knownGene"
        if bsgenome is None:
            bsgenome = f"BSgenome.Hsapiens.UCSC.{self._genome}"
        self._coding = predict_coding(self._vcf, txdb=txdb, bsgenome=bsgenome)
        return self

    def to_dataframe(self):
        """Return the most recent result as a DataFrame."""
        if self._coding is not None:
            return self._coding
        if self._locations is not None:
            return self._locations
        from .wrappers.variant_annotation import vcf_to_dataframe
        return vcf_to_dataframe(self._vcf)

    def report(self):
        """Print a summary report."""
        # ... location/consequence summary ...
        pass
```

## Dependencies

```python
# In _deps.py, add:
# "VariantAnnotation" → BiocManager::install("VariantAnnotation")
# TxDb/BSgenome packages are resolved dynamically by genome name
```

## Why This Fits

1. **Same pattern** as DESeq2/edgeR wrappers — read data, call R, return DataFrame
2. **Pandas in, pandas out** — VCF becomes a DataFrame; annotations are DataFrames
3. **Genome-aware defaults** — `genome="hg38"` auto-resolves TxDb and BSgenome
4. **`.report()`** summarizes variant consequences (like DESeq2 summary)
5. **No R code for users** — all the complexity is hidden behind Python calls
6. **You maintain VariantAnnotation** — perfect dogfooding loop

## Effort Estimate

- Tier 3 functional API: ~2 hours (straightforward rpy2 wrapping)
- Tier 1 quick API: ~1 hour (compose Tier 3 functions)
- Tier 2 class API: ~2 hours (state management, chaining)
- Tests: ~2 hours
- Documentation: ~1 hour

**Total: ~8 hours to a working module with full test coverage.**
