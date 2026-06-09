"""Rosetta exception classes for R error translation."""


class RosettaError(Exception):
    """Base exception for rosetta."""


class RPackageMissing(RosettaError):
    """Required R package is not installed."""

    def __init__(self, package: str):
        self.package = package
        super().__init__(f"R package '{package}' is not installed. Install with: R -e 'BiocManager::install(\"{package}\")'")


class RFormulaError(RosettaError):
    """Invalid R design formula."""


class RDataError(RosettaError):
    """Incompatible input data for R function."""
class RosettaSecurityError(RosettaError):
    """Exception raised for security-related issues."""

