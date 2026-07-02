"""VCF — Tier 2 class-based API for VariantAnnotation.

Stateful VCF analysis object with chainable methods.
Wraps the Tier 3 functional API in rosetta.wrappers.variant_annotation.
"""

from __future__ import annotations

import os

import pandas as pd

from .._errors import RDataError


class VCF:
    """Stateful VCF analysis object with chainable methods.

    Example:
        >>> vcf = rb.VCF("sample.vcf.gz", genome="hg38")
        >>> results = vcf.locate_variants().predict_coding().to_dataframe()
        >>> vcf.report()
    """

    def __init__(self, filepath: str, genome: str = "hg38", param: dict | None = None):
        """Read a VCF file and store the R VCF object.

        Args:
            filepath: Path to VCF/VCF.gz/BCF file.
            genome: Genome assembly name (e.g., "hg38", "hg19").
            param: Dict with keys 'info', 'geno', 'samples', 'which' to filter.

        Raises:
            RDataError: If the file does not exist.
        """
        filepath = str(filepath)
        if not os.path.isfile(filepath):
            raise RDataError(f"VCF file not found: {filepath}")

        from .variant_annotation import read_vcf

        self._vcf = read_vcf(filepath, genome=genome, param=param)
        self._filepath = filepath
        self._genome = genome
        self._locations: pd.DataFrame | None = None
        self._coding: pd.DataFrame | None = None

    def locate_variants(self, txdb=None, region: str = "all") -> "VCF":
        """Annotate variant genomic locations (coding, intronic, UTR, etc.).

        Args:
            txdb: TxDb package name (str) or R TxDb object.
                  Defaults to TxDb.Hsapiens.UCSC.<genome>.knownGene.
            region: One of 'all', 'coding', 'intron', 'fiveUTR', 'threeUTR',
                    'intergenic', 'spliceSite', 'promoter'.

        Returns:
            self (for method chaining).
        """
        from .variant_annotation import locate_variants as _locate_variants

        if txdb is None:
            txdb = f"TxDb.Hsapiens.UCSC.{self._genome}.knownGene"

        self._locations = _locate_variants(self._vcf, txdb=txdb, region=region)
        return self

    def predict_coding(self, txdb=None, bsgenome=None) -> "VCF":
        """Predict amino acid coding changes for variants in coding regions.

        Args:
            txdb: TxDb package name (str) or R TxDb object.
                  Defaults to TxDb.Hsapiens.UCSC.<genome>.knownGene.
            bsgenome: BSgenome package name (str) or R BSgenome object.
                      Defaults to BSgenome.Hsapiens.UCSC.<genome>.

        Returns:
            self (for method chaining).
        """
        from .variant_annotation import predict_coding as _predict_coding

        if txdb is None:
            txdb = f"TxDb.Hsapiens.UCSC.{self._genome}.knownGene"
        if bsgenome is None:
            bsgenome = f"BSgenome.Hsapiens.UCSC.{self._genome}"

        self._coding = _predict_coding(self._vcf, txdb=txdb, bsgenome=bsgenome)
        return self

    def to_dataframe(self) -> pd.DataFrame:
        """Return the most recent result as a DataFrame.

        Priority: coding > locations > raw VCF fixed fields + INFO.
        """
        if self._coding is not None:
            return self._coding
        if self._locations is not None:
            return self._locations
        from .variant_annotation import vcf_to_dataframe

        return vcf_to_dataframe(self._vcf)

    def header(self) -> dict:
        """Return VCF header metadata (info, geno, samples).

        Returns:
            Dict with keys 'info', 'geno', 'samples'.
        """
        from .variant_annotation import scan_vcf_header

        return scan_vcf_header(self._filepath)

    def report(self) -> None:
        """Print a summary report of variant annotations."""
        lines = []
        lines.append("Variant Annotation Summary")
        lines.append("─" * 30)

        if self._locations is not None:
            total = len(self._locations)
            lines.append(f"Total annotated locations:  {total:,}")
            if "LOCATION" in self._locations.columns:
                counts = self._locations["LOCATION"].value_counts()
                for loc, count in counts.items():
                    pct = 100 * count / total if total > 0 else 0
                    lines.append(f"  {loc:<22} {count:>6} ({pct:.1f}%)")

        if self._coding is not None:
            total = len(self._coding)
            lines.append(f"Coding predictions:        {total:,}")
            if "CONSEQUENCE" in self._coding.columns:
                counts = self._coding["CONSEQUENCE"].value_counts()
                for cons, count in counts.items():
                    lines.append(f"  {cons:<22} {count:>6}")

        if self._locations is None and self._coding is None:
            # No annotations yet — report raw variant count
            from .variant_annotation import vcf_to_dataframe

            df = vcf_to_dataframe(self._vcf)
            lines.append(f"Total variants:            {len(df):,}")
            lines.append("(No annotations computed yet — call .locate_variants() or .predict_coding())")

        print("\n".join(lines))
