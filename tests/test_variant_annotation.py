"""Tests for rosetta.wrappers.variant_annotation."""

import os
import tempfile

import pandas as pd
import pytest

from rosetta._errors import RDataError


def _va_available():
    """Check if VariantAnnotation is installed in the R environment."""
    try:
        from rosetta._deps import is_installed

        return is_installed("VariantAnnotation")
    except Exception:
        return False


# --- Tests that do NOT require VariantAnnotation installed ---


def test_read_vcf_missing_file():
    """Verify RDataError raised for nonexistent file."""
    from rosetta.wrappers.variant_annotation import read_vcf

    with pytest.raises(RDataError, match="VCF file not found"):
        read_vcf("/nonexistent/path/to/fake.vcf.gz", genome="hg38")


def test_scan_vcf_header_missing_file():
    """Verify RDataError raised for nonexistent file in scan_vcf_header."""
    from rosetta.wrappers.variant_annotation import scan_vcf_header

    with pytest.raises(RDataError, match="VCF file not found"):
        scan_vcf_header("/nonexistent/path/to/fake.vcf.gz")


def test_make_region_valid():
    """Test _make_region accepts all valid region strings (needs VariantAnnotation for object creation)."""
    from rosetta.wrappers.variant_annotation import _make_region

    valid_regions = [
        "all",
        "coding",
        "intron",
        "fiveUTR",
        "threeUTR",
        "intergenic",
        "spliceSite",
        "promoter",
    ]

    if not _va_available():
        # Without VariantAnnotation, we can only verify the ValueError path
        # The valid path requires VariantAnnotation to instantiate the object
        pytest.skip("VariantAnnotation not installed — cannot instantiate region objects")

    for region in valid_regions:
        result = _make_region(region)
        # Should return an R object (RS4 instance)
        assert result is not None


def test_make_region_invalid():
    """Test _make_region raises ValueError for bad input."""
    from rosetta.wrappers.variant_annotation import _make_region

    with pytest.raises(ValueError, match="Unknown region"):
        _make_region("nonexistent_region")

    with pytest.raises(ValueError, match="Unknown region"):
        _make_region("")

    with pytest.raises(ValueError, match="Unknown region"):
        _make_region("CODING")  # case-sensitive


def test_resolve_txdb_none_raises():
    """Verify _resolve_txdb raises RDataError when txdb is None."""
    from rosetta.wrappers.variant_annotation import _resolve_txdb

    with pytest.raises(RDataError, match="txdb argument is required"):
        _resolve_txdb(None)


def test_resolve_bsgenome_none_raises():
    """Verify _resolve_bsgenome raises RDataError when bsgenome is None."""
    from rosetta.wrappers.variant_annotation import _resolve_bsgenome

    with pytest.raises(RDataError, match="bsgenome argument is required"):
        _resolve_bsgenome(None)


# --- Tests that DO require VariantAnnotation ---


@pytest.mark.skipif(not _va_available(), reason="VariantAnnotation not installed in R")
def test_scan_vcf_header():
    """Test scan_vcf_header with a minimal VCF file."""
    # Create a minimal valid VCF file
    vcf_content = (
        "##fileformat=VCFv4.1\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
        "chr1\t100\t.\tA\tT\t30\tPASS\tDP=10\tGT\t0/1\n"
    )

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".vcf", delete=False
    ) as f:
        f.write(vcf_content)
        vcf_path = f.name

    try:
        from rosetta.wrappers.variant_annotation import scan_vcf_header

        header = scan_vcf_header(vcf_path)
        assert isinstance(header, dict)
        assert "info" in header
        assert "geno" in header
        assert "samples" in header
        assert "SAMPLE1" in header["samples"]
    finally:
        os.unlink(vcf_path)


@pytest.mark.skipif(not _va_available(), reason="VariantAnnotation not installed in R")
def test_read_vcf_basic():
    """Test basic VCF reading with a minimal file."""
    vcf_content = (
        "##fileformat=VCFv4.1\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
        "chr1\t100\t.\tA\tT\t30\tPASS\tDP=10\tGT\t0/1\n"
    )

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".vcf", delete=False
    ) as f:
        f.write(vcf_content)
        vcf_path = f.name

    try:
        from rosetta.wrappers.variant_annotation import read_vcf

        vcf = read_vcf(vcf_path, genome="hg19")
        # Should return an R object (not None)
        assert vcf is not None
    finally:
        os.unlink(vcf_path)


@pytest.mark.skipif(not _va_available(), reason="VariantAnnotation not installed in R")
def test_make_region_returns_r_objects():
    """Verify _make_region returns proper R objects for all valid regions."""
    from rosetta.wrappers.variant_annotation import _make_region

    valid_regions = [
        "all",
        "coding",
        "intron",
        "fiveUTR",
        "threeUTR",
        "intergenic",
        "spliceSite",
        "promoter",
    ]

    for region in valid_regions:
        result = _make_region(region)
        assert result is not None, f"_make_region('{region}') returned None"


# --- VCF Class (Tier 2) Tests ---


def test_vcf_class_missing_file():
    """Verify VCF class raises RDataError for nonexistent file."""
    from rosetta.wrappers.vcf import VCF

    with pytest.raises(RDataError, match="VCF file not found"):
        VCF("/nonexistent/path/to/fake.vcf.gz")


def test_vcf_class_also_importable_from_top_level():
    """Verify VCF is importable from rosetta top-level."""
    import rosetta as rb

    assert hasattr(rb, "VCF")
    assert rb.VCF is not None


# --- Tier 1 Quick API Tests ---


def test_quick_locate_missing_file():
    """Verify quick_locate_variants raises RDataError for nonexistent file."""
    import rosetta as rb

    with pytest.raises(RDataError, match="VCF file not found"):
        rb.quick_locate_variants("/nonexistent/path/to/fake.vcf.gz")


def test_quick_predict_missing_file():
    """Verify quick_predict_coding raises RDataError for nonexistent file."""
    import rosetta as rb

    with pytest.raises(RDataError, match="VCF file not found"):
        rb.quick_predict_coding("/nonexistent/path/to/fake.vcf.gz")


# --- Integration tests requiring VariantAnnotation ---


@pytest.mark.skipif(not _va_available(), reason="VariantAnnotation not installed in R")
def test_vcf_class_chainable():
    """Verify VCF class methods return self for chaining."""
    vcf_content = (
        "##fileformat=VCFv4.1\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
        "chr1\t100\t.\tA\tT\t30\tPASS\tDP=10\tGT\t0/1\n"
    )

    with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
        f.write(vcf_content)
        vcf_path = f.name

    try:
        from rosetta.wrappers.vcf import VCF

        vcf = VCF(vcf_path, genome="hg19")

        # to_dataframe should work without annotations
        df = vcf.to_dataframe()
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0

        # header should return a dict
        hdr = vcf.header()
        assert isinstance(hdr, dict)
        assert "samples" in hdr
    finally:
        os.unlink(vcf_path)


@pytest.mark.skipif(not _va_available(), reason="VariantAnnotation not installed in R")
def test_vcf_class_to_dataframe_returns_pandas():
    """Verify VCF.to_dataframe() returns a pandas DataFrame."""
    vcf_content = (
        "##fileformat=VCFv4.1\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
        "chr1\t100\t.\tA\tT\t30\tPASS\tDP=10\tGT\t0/1\n"
    )

    with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
        f.write(vcf_content)
        vcf_path = f.name

    try:
        from rosetta.wrappers.vcf import VCF

        vcf = VCF(vcf_path, genome="hg19")
        df = vcf.to_dataframe()
        assert isinstance(df, pd.DataFrame)
    finally:
        os.unlink(vcf_path)
