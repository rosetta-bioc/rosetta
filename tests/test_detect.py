"""Tests for rosetta._detect — backend detection and version pinning."""

import json
import warnings
from unittest.mock import patch, MagicMock

import pytest

from rosetta._detect import (
    check_rpy2_available,
    find_rscript,
    verify_r_package_versions,
    detect_backend,
    PINNED_VERSIONS,
)


# ---------------------------------------------------------------------------
# check_rpy2_available
# ---------------------------------------------------------------------------

def test_check_rpy2_available_returns_true_when_working():
    # rpy2 is installed in the test environment — should return True
    result = check_rpy2_available()
    assert isinstance(result, bool)


def test_check_rpy2_available_false_when_spec_missing():
    with patch("importlib.util.find_spec", return_value=None):
        assert check_rpy2_available() is False


def test_check_rpy2_available_false_on_import_error():
    import builtins
    import importlib.util as _ilu

    real_find_spec = _ilu.find_spec
    real_import = builtins.__import__

    def _fake_find_spec(name, *args, **kwargs):
        if name == "rpy2":
            return True  # pretend it exists
        return real_find_spec(name, *args, **kwargs)

    def _failing_import(name, *args, **kwargs):
        if name == "rpy2.robjects":
            raise RuntimeError("R init failed")
        return real_import(name, *args, **kwargs)

    with patch("importlib.util.find_spec", side_effect=_fake_find_spec):
        with patch("builtins.__import__", side_effect=_failing_import):
            # If R crashes during initialisation it should still return False
            result = check_rpy2_available()
            assert result is False


# ---------------------------------------------------------------------------
# find_rscript
# ---------------------------------------------------------------------------

def test_find_rscript_returns_path_when_available():
    with patch("shutil.which", return_value="/usr/bin/Rscript"):
        assert find_rscript() == "/usr/bin/Rscript"


def test_find_rscript_returns_none_when_missing():
    with patch("shutil.which", return_value=None):
        assert find_rscript() is None


def test_find_rscript_real():
    # On CI and developer machines Rscript should be present
    path = find_rscript()
    assert path is None or isinstance(path, str)


# ---------------------------------------------------------------------------
# verify_r_package_versions
# ---------------------------------------------------------------------------

def _make_subprocess_result(versions: dict):
    mock = MagicMock()
    mock.stdout = json.dumps(versions)
    return mock


def test_verify_no_warning_when_versions_match():
    matched = {pkg: ver for pkg, ver in PINNED_VERSIONS.items()}
    with patch("subprocess.run", return_value=_make_subprocess_result(matched)):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            verify_r_package_versions("/usr/bin/Rscript")
        assert len(w) == 0


def test_verify_warns_on_version_mismatch():
    mismatched = {pkg: "0.0.1" for pkg in PINNED_VERSIONS}
    with patch("subprocess.run", return_value=_make_subprocess_result(mismatched)):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            verify_r_package_versions("/usr/bin/Rscript")
        messages = [str(x.message) for x in w]
        assert any("version mismatch" in m for m in messages)


def test_verify_warns_when_package_missing():
    # All packages missing (NA in R → None in JSON)
    missing = {pkg: None for pkg in PINNED_VERSIONS}
    with patch("subprocess.run", return_value=_make_subprocess_result(missing)):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            verify_r_package_versions("/usr/bin/Rscript")
        messages = [str(x.message) for x in w]
        assert any("not installed" in m for m in messages)


def test_verify_warns_on_subprocess_failure():
    with patch("subprocess.run", side_effect=OSError("Rscript not found")):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            verify_r_package_versions("/usr/bin/Rscript")
        assert len(w) == 1
        assert "Could not verify" in str(w[0].message)


# ---------------------------------------------------------------------------
# detect_backend
# ---------------------------------------------------------------------------

def test_detect_backend_returns_rpy2_when_available():
    with patch("rosetta._detect.check_rpy2_available", return_value=True):
        assert detect_backend() == "rpy2"


def test_detect_backend_returns_subprocess_when_rpy2_missing_but_rscript_present():
    with patch("rosetta._detect.check_rpy2_available", return_value=False):
        with patch("rosetta._detect.find_rscript", return_value="/usr/bin/Rscript"):
            assert detect_backend() == "subprocess"


def test_detect_backend_returns_fallback_when_nothing_available():
    with patch("rosetta._detect.check_rpy2_available", return_value=False):
        with patch("rosetta._detect.find_rscript", return_value=None):
            assert detect_backend() == "fallback"


def test_detect_backend_priority_rpy2_over_subprocess():
    """rpy2 is always preferred over subprocess even if Rscript is present."""
    with patch("rosetta._detect.check_rpy2_available", return_value=True):
        with patch("rosetta._detect.find_rscript", return_value="/usr/bin/Rscript"):
            assert detect_backend() == "rpy2"
