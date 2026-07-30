"""Tests for rosetta._deps."""

from unittest.mock import patch

import pytest

from rosetta._deps import is_installed, ensure_installed
from rosetta._errors import RPackageMissing

if not is_installed("stats"):
    pytest.skip("Base R 'stats' package not available", allow_module_level=True)

def test_is_installed_base():
    """base R package 'stats' should always be installed."""
    assert is_installed("stats") is True


def test_is_installed_missing():
    assert is_installed("nonexistent_fake_pkg_xyz") is False


def test_ensure_installed_raises():
    with pytest.raises(RPackageMissing, match="nonexistent_fake_pkg_xyz"):
        ensure_installed("nonexistent_fake_pkg_xyz")


def test_ensure_installed_passes():
    ensure_installed("stats")  # should not raise
