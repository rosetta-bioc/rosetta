"""Tests for VCF class-based wrapper and variant_annotation functional API."""
from __future__ import annotations

import os
import pandas as pd
import pytest
from unittest.mock import MagicMock, patch, call


# =========================================================================
# Helpers
# =========================================================================

def _va_patches(va_pkg=None, ro_mock=None, extra=None):
    """Common patches for variant_annotation module."""
    va_pkg = va_pkg or MagicMock(name="VariantAnnotation")
    ro_mock = ro_mock or MagicMock(name="ro")
    patches = [
        patch("rosetta.wrappers.variant_annotation.ensure_installed"),
        patch("rosetta.wrappers.variant_annotation.importr", return_value=va_pkg),
        patch("rosetta.wrappers.variant_annotation.localconverter"),
        patch("rosetta.wrappers.variant_annotation._converter", MagicMock()),
        patch("rosetta.wrappers.variant_annotation.ro", ro_mock),
        patch("rosetta.wrappers.variant_annotation.to_pandas", return_value=pd.DataFrame({"x": [1]})),
        patch("rosetta.wrappers.variant_annotation.to_r_df", return_value=MagicMock()),
        patch("rosetta.wrappers.variant_annotation.codegen", MagicMock()),
    ]
    if extra:
        patches.extend(extra)
    return patches


# =========================================================================
# variant_annotation — read_vcf (lines 54-56)
# =========================================================================

def test_read_vcf_file_not_found():
    from rosetta.wrappers.variant_annotation import read_vcf
    from rosetta._errors import RDataError
    with pytest.raises(RDataError, match="not found"):
        read_vcf("/nonexistent/path.vcf.gz")


def test_read_vcf_no_param(tmp_path):
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text("##fileformat=VCFv4.2\n")

    va_pkg = MagicMock()
    va_pkg.readVcf.return_value = MagicMock()

    with patch("rosetta.wrappers.variant_annotation.ensure_installed"), \
         patch("rosetta.wrappers.variant_annotation.importr", return_value=va_pkg), \
         patch("rosetta.wrappers.variant_annotation.localconverter"), \
         patch("rosetta.wrappers.variant_annotation._converter", MagicMock()), \
         patch("rosetta.wrappers.variant_annotation.codegen", MagicMock()):
        from rosetta.wrappers.variant_annotation import read_vcf
        result = read_vcf(str(vcf_file), genome="hg38")

    va_pkg.readVcf.assert_called_once_with(str(vcf_file), "hg38")


def test_read_vcf_with_param(tmp_path):
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text("##fileformat=VCFv4.2\n")

    va_pkg = MagicMock()
    va_pkg.readVcf.return_value = MagicMock()
    mock_ro = MagicMock()

    with patch("rosetta.wrappers.variant_annotation.ensure_installed"), \
         patch("rosetta.wrappers.variant_annotation.importr", return_value=va_pkg), \
         patch("rosetta.wrappers.variant_annotation.localconverter"), \
         patch("rosetta.wrappers.variant_annotation._converter", MagicMock()), \
         patch("rosetta.wrappers.variant_annotation.codegen", MagicMock()), \
         patch("rosetta.wrappers.variant_annotation.ro", mock_ro), \
         patch("rosetta.wrappers.variant_annotation._build_scan_param", return_value=MagicMock()) as mock_bsp:
        from rosetta.wrappers.variant_annotation import read_vcf
        read_vcf(str(vcf_file), genome="hg38", param={"info": ["AF"]})

    mock_bsp.assert_called_once_with({"info": ["AF"]})
    # called with param kwarg
    assert va_pkg.readVcf.call_count == 1


# =========================================================================
# variant_annotation — locate_variants / predict_coding (lines 194-230)
# =========================================================================

def test_locate_variants_calls_r():
    from rosetta.wrappers.variant_annotation import locate_variants
    va_pkg = MagicMock()
    base_pkg = MagicMock()
    txdb_obj = MagicMock()
    region_obj = MagicMock()
    va_pkg.locateVariants.return_value = MagicMock()

    with patch("rosetta.wrappers.variant_annotation.ensure_installed"), \
         patch("rosetta.wrappers.variant_annotation.importr", side_effect=[va_pkg, base_pkg]), \
         patch("rosetta.wrappers.variant_annotation.localconverter"), \
         patch("rosetta.wrappers.variant_annotation._converter", MagicMock()), \
         patch("rosetta.wrappers.variant_annotation.codegen", MagicMock()), \
         patch("rosetta.wrappers.variant_annotation.to_pandas", return_value=pd.DataFrame()), \
         patch("rosetta.wrappers.variant_annotation.to_r_df", return_value=MagicMock()), \
         patch("rosetta.wrappers.variant_annotation._resolve_txdb", return_value=txdb_obj), \
         patch("rosetta.wrappers.variant_annotation._make_region", return_value=region_obj):
        locate_variants(MagicMock(), txdb="TxDb.Hsapiens.UCSC.hg38.knownGene", region="coding")

    va_pkg.locateVariants.assert_called_once()


def test_predict_coding_calls_r():
    from rosetta.wrappers.variant_annotation import predict_coding
    va_pkg = MagicMock()
    base_pkg = MagicMock()
    txdb_obj = MagicMock()
    bsg_obj = MagicMock()
    va_pkg.predictCoding.return_value = MagicMock()

    with patch("rosetta.wrappers.variant_annotation.ensure_installed"), \
         patch("rosetta.wrappers.variant_annotation.importr", side_effect=[va_pkg, base_pkg]), \
         patch("rosetta.wrappers.variant_annotation.localconverter"), \
         patch("rosetta.wrappers.variant_annotation._converter", MagicMock()), \
         patch("rosetta.wrappers.variant_annotation.codegen", MagicMock()), \
         patch("rosetta.wrappers.variant_annotation.to_pandas", return_value=pd.DataFrame()), \
         patch("rosetta.wrappers.variant_annotation.to_r_df", return_value=MagicMock()), \
         patch("rosetta.wrappers.variant_annotation._resolve_txdb", return_value=txdb_obj), \
         patch("rosetta.wrappers.variant_annotation._resolve_bsgenome", return_value=bsg_obj):
        predict_coding(MagicMock(), txdb="TxDb.pkg", bsgenome="BSgenome.pkg")

    va_pkg.predictCoding.assert_called_once()


# =========================================================================
# variant_annotation — scan_vcf_header (lines 169-179)
# =========================================================================

def test_scan_vcf_header_file_not_found():
    from rosetta.wrappers.variant_annotation import scan_vcf_header
    from rosetta._errors import RDataError
    with pytest.raises(RDataError, match="not found"):
        scan_vcf_header("/no/such/file.vcf")


def test_scan_vcf_header_returns_dict(tmp_path):
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text("##fileformat=VCFv4.2\n")

    va_pkg = MagicMock()
    base_pkg = MagicMock()
    mock_ro = MagicMock()
    mock_ro.r.return_value = MagicMock(return_value=MagicMock())
    info_df = pd.DataFrame({"Number": ["1"], "Type": ["Float"], "Description": ["x"]})
    geno_df = pd.DataFrame({"Number": ["1"], "Type": ["Integer"]})

    with patch("rosetta.wrappers.variant_annotation.ensure_installed"), \
         patch("rosetta.wrappers.variant_annotation.importr", side_effect=[va_pkg, base_pkg]), \
         patch("rosetta.wrappers.variant_annotation.localconverter"), \
         patch("rosetta.wrappers.variant_annotation._converter", MagicMock()), \
         patch("rosetta.wrappers.variant_annotation.codegen", MagicMock()), \
         patch("rosetta.wrappers.variant_annotation.ro", mock_ro), \
         patch("rosetta.wrappers.variant_annotation.to_r_df", return_value=MagicMock()), \
         patch("rosetta.wrappers.variant_annotation.to_pandas", side_effect=[info_df, geno_df]):
        from rosetta.wrappers.variant_annotation import scan_vcf_header
        result = scan_vcf_header(str(vcf_file))

    assert "info" in result
    assert "geno" in result
    assert "samples" in result


# =========================================================================
# variant_annotation — _build_scan_param (lines 244-275)
# =========================================================================

def test_build_scan_param_info_geno_samples():
    from rosetta.wrappers.variant_annotation import _build_scan_param
    va_pkg = MagicMock()
    mock_ro = MagicMock()
    mock_ro.StrVector.side_effect = lambda x: x
    mock_ro.NA_Character = "NA"

    with patch("rosetta.wrappers.variant_annotation.ensure_installed"), \
         patch("rosetta.wrappers.variant_annotation.importr", return_value=va_pkg), \
         patch("rosetta.wrappers.variant_annotation.ro", mock_ro):
        _build_scan_param({"info": ["AF", "DP"], "geno": ["GT"], "samples": ["S1"]})

    va_pkg.ScanVcfParam.assert_called_once()
    kwargs = va_pkg.ScanVcfParam.call_args[1]
    assert "info" in kwargs
    assert "geno" in kwargs
    assert "samples" in kwargs


def test_build_scan_param_empty_info_uses_na():
    from rosetta.wrappers.variant_annotation import _build_scan_param
    va_pkg = MagicMock()
    mock_ro = MagicMock()
    mock_ro.NA_Character = "NA_CHAR"

    with patch("rosetta.wrappers.variant_annotation.ensure_installed"), \
         patch("rosetta.wrappers.variant_annotation.importr", return_value=va_pkg), \
         patch("rosetta.wrappers.variant_annotation.ro", mock_ro):
        _build_scan_param({"info": [], "geno": []})

    kwargs = va_pkg.ScanVcfParam.call_args[1]
    assert kwargs["info"] == "NA_CHAR"
    assert kwargs["geno"] == "NA_CHAR"


def test_build_scan_param_which():
    from rosetta.wrappers.variant_annotation import _build_scan_param
    va_pkg = MagicMock()
    gr_pkg = MagicMock()
    ir_pkg = MagicMock()
    mock_ro = MagicMock()
    mock_ro.StrVector.side_effect = lambda x: x
    mock_ro.IntVector.side_effect = lambda x: x

    with patch("rosetta.wrappers.variant_annotation.ensure_installed"), \
         patch("rosetta.wrappers.variant_annotation.importr", side_effect=[va_pkg, gr_pkg, ir_pkg]), \
         patch("rosetta.wrappers.variant_annotation.ro", mock_ro):
        _build_scan_param({"which": {"chr1": (1000, 5000), "chr2": (2000, 8000)}})

    gr_pkg.GRanges.assert_called_once()
    ir_pkg.IRanges.assert_called_once()


# =========================================================================
# variant_annotation — _resolve_txdb / _resolve_bsgenome (lines 281-311)
# =========================================================================

def test_resolve_txdb_none_raises():
    from rosetta.wrappers.variant_annotation import _resolve_txdb
    from rosetta._errors import RDataError
    with pytest.raises(RDataError, match="txdb argument is required"):
        _resolve_txdb(None)


def test_resolve_txdb_string_loads_package():
    from rosetta.wrappers.variant_annotation import _resolve_txdb
    mock_ro = MagicMock()
    mock_ro.r.return_value = MagicMock()

    with patch("rosetta.wrappers.variant_annotation.ensure_installed"), \
         patch("rosetta.wrappers.variant_annotation.codegen", MagicMock()), \
         patch("rosetta.wrappers.variant_annotation.ro", mock_ro):
        _resolve_txdb("TxDb.Hsapiens.UCSC.hg38.knownGene")

    mock_ro.r.assert_called_once()


def test_resolve_txdb_object_passthrough():
    from rosetta.wrappers.variant_annotation import _resolve_txdb
    obj = MagicMock()
    result = _resolve_txdb(obj)
    assert result is obj


def test_resolve_bsgenome_none_raises():
    from rosetta.wrappers.variant_annotation import _resolve_bsgenome
    from rosetta._errors import RDataError
    with pytest.raises(RDataError, match="bsgenome argument is required"):
        _resolve_bsgenome(None)


def test_resolve_bsgenome_string_loads_package():
    from rosetta.wrappers.variant_annotation import _resolve_bsgenome
    mock_ro = MagicMock()
    mock_ro.r.return_value = MagicMock()

    with patch("rosetta.wrappers.variant_annotation.ensure_installed"), \
         patch("rosetta.wrappers.variant_annotation.codegen", MagicMock()), \
         patch("rosetta.wrappers.variant_annotation.ro", mock_ro):
        _resolve_bsgenome("BSgenome.Hsapiens.UCSC.hg38")

    mock_ro.r.assert_called_once()


def test_resolve_bsgenome_object_passthrough():
    from rosetta.wrappers.variant_annotation import _resolve_bsgenome
    obj = MagicMock()
    assert _resolve_bsgenome(obj) is obj


# =========================================================================
# variant_annotation — _make_region (lines ~315-340)
# =========================================================================

def test_make_region_valid_keys():
    from rosetta.wrappers.variant_annotation import _make_region
    va_pkg = MagicMock()
    for region in ("all", "coding", "intron", "fiveUTR", "threeUTR",
                   "intergenic", "spliceSite", "promoter"):
        with patch("rosetta.wrappers.variant_annotation.importr", return_value=va_pkg), \
             patch("rosetta.wrappers.variant_annotation.ensure_installed"):
            _make_region(region)


def test_make_region_invalid_raises():
    from rosetta.wrappers.variant_annotation import _make_region
    with pytest.raises(ValueError, match="Unknown region"):
        _make_region("badregion")


# =========================================================================
# VCF class (vcf.py lines 60-144)
# =========================================================================

def _make_vcf(tmp_path, genome="hg38"):
    vcf_file = tmp_path / "test.vcf.gz"
    vcf_file.write_text("##fileformat=VCFv4.2\n")
    vcf_r = MagicMock(name="vcf_r")

    with patch("rosetta.wrappers.vcf.VCF.__init__", lambda self, *a, **k: None):
        from rosetta.wrappers.vcf import VCF
        obj = VCF.__new__(VCF)
        obj._vcf = vcf_r
        obj._filepath = str(vcf_file)
        obj._genome = genome
        obj._locations = None
        obj._coding = None
    return obj, vcf_r, str(vcf_file)


def test_vcf_init_file_not_found():
    from rosetta.wrappers.vcf import VCF
    from rosetta._errors import RDataError
    with pytest.raises(RDataError, match="not found"):
        VCF("/no/such/file.vcf.gz")


def test_vcf_init_calls_read_vcf(tmp_path):
    vcf_file = tmp_path / "test.vcf.gz"
    vcf_file.write_text("##fileformat=VCFv4.2\n")
    mock_vcf = MagicMock()

    # VCF.__init__ does `from .variant_annotation import read_vcf` — patch at source
    with patch("rosetta.wrappers.variant_annotation.read_vcf", return_value=mock_vcf) as mock_read:
        from rosetta.wrappers.vcf import VCF
        obj = VCF(str(vcf_file), genome="hg19", param={"info": ["AF"]})

    mock_read.assert_called_once_with(str(vcf_file), genome="hg19", param={"info": ["AF"]})
    assert obj._vcf is mock_vcf
    assert obj._genome == "hg19"


def test_vcf_locate_variants_default_txdb(tmp_path):
    obj, vcf_r, path = _make_vcf(tmp_path)
    loc_df = pd.DataFrame({"LOCATION": ["coding", "intron"]})

    # method does `from .variant_annotation import locate_variants as _locate_variants`
    with patch("rosetta.wrappers.variant_annotation.locate_variants", return_value=loc_df) as mock_lv:
        result = obj.locate_variants()

    mock_lv.assert_called_once_with(vcf_r, txdb="TxDb.Hsapiens.UCSC.hg38.knownGene", region="all")
    assert result is obj
    assert obj._locations is loc_df


def test_vcf_locate_variants_custom_txdb(tmp_path):
    obj, vcf_r, _ = _make_vcf(tmp_path)
    custom_txdb = MagicMock()

    with patch("rosetta.wrappers.variant_annotation.locate_variants", return_value=pd.DataFrame()) as mock_lv:
        obj.locate_variants(txdb=custom_txdb, region="coding")

    mock_lv.assert_called_once_with(vcf_r, txdb=custom_txdb, region="coding")


def test_vcf_predict_coding_default_packages(tmp_path):
    obj, vcf_r, _ = _make_vcf(tmp_path)
    coding_df = pd.DataFrame({"CONSEQUENCE": ["nonsynonymous"]})

    with patch("rosetta.wrappers.variant_annotation.predict_coding", return_value=coding_df) as mock_pc:
        result = obj.predict_coding()

    mock_pc.assert_called_once_with(
        vcf_r,
        txdb="TxDb.Hsapiens.UCSC.hg38.knownGene",
        bsgenome="BSgenome.Hsapiens.UCSC.hg38",
    )
    assert result is obj
    assert obj._coding is coding_df


def test_vcf_predict_coding_custom_packages(tmp_path):
    obj, vcf_r, _ = _make_vcf(tmp_path, genome="hg19")
    txdb = MagicMock()
    bsg = MagicMock()

    with patch("rosetta.wrappers.variant_annotation.predict_coding", return_value=pd.DataFrame()) as mock_pc:
        obj.predict_coding(txdb=txdb, bsgenome=bsg)

    mock_pc.assert_called_once_with(vcf_r, txdb=txdb, bsgenome=bsg)


def test_vcf_to_dataframe_priority_coding(tmp_path):
    obj, _, _ = _make_vcf(tmp_path)
    obj._coding = pd.DataFrame({"CONSEQUENCE": ["nonsynonymous"]})
    obj._locations = pd.DataFrame({"LOCATION": ["coding"]})

    result = obj.to_dataframe()
    assert "CONSEQUENCE" in result.columns


def test_vcf_to_dataframe_priority_locations(tmp_path):
    obj, vcf_r, _ = _make_vcf(tmp_path)
    obj._locations = pd.DataFrame({"LOCATION": ["intron"]})

    result = obj.to_dataframe()
    assert "LOCATION" in result.columns


def test_vcf_to_dataframe_raw_fallback(tmp_path):
    obj, vcf_r, _ = _make_vcf(tmp_path)
    raw_df = pd.DataFrame({"CHROM": ["chr1"], "POS": [100]})

    # method does `from .variant_annotation import vcf_to_dataframe`
    with patch("rosetta.wrappers.variant_annotation.vcf_to_dataframe", return_value=raw_df) as mock_vtd:
        result = obj.to_dataframe()

    mock_vtd.assert_called_once_with(vcf_r)
    assert "CHROM" in result.columns


def test_vcf_header(tmp_path):
    obj, _, path = _make_vcf(tmp_path)
    header_dict = {"info": pd.DataFrame(), "geno": pd.DataFrame(), "samples": ["S1"]}

    # method does `from .variant_annotation import scan_vcf_header`
    with patch("rosetta.wrappers.variant_annotation.scan_vcf_header", return_value=header_dict) as mock_sh:
        result = obj.header()

    mock_sh.assert_called_once_with(path)
    assert result["samples"] == ["S1"]


def test_vcf_report_no_annotations(tmp_path, capsys):
    obj, vcf_r, _ = _make_vcf(tmp_path)
    raw_df = pd.DataFrame({"CHROM": ["chr1"] * 5})

    with patch("rosetta.wrappers.variant_annotation.vcf_to_dataframe", return_value=raw_df):
        obj.report()

    out = capsys.readouterr().out
    assert "5" in out
    assert "No annotations" in out


def test_vcf_report_with_locations(tmp_path, capsys):
    obj, _, _ = _make_vcf(tmp_path)
    obj._locations = pd.DataFrame({"LOCATION": ["coding", "coding", "intron"]})

    obj.report()
    out = capsys.readouterr().out
    assert "3" in out
    assert "coding" in out


def test_vcf_report_with_coding(tmp_path, capsys):
    obj, _, _ = _make_vcf(tmp_path)
    obj._coding = pd.DataFrame({"CONSEQUENCE": ["nonsynonymous", "synonymous"]})

    obj.report()
    out = capsys.readouterr().out
    assert "2" in out
    assert "nonsynonymous" in out
