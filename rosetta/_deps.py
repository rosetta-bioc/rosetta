"""R package detection and installation via BiocManager."""

from ._bridge import ACTIVE_BACKEND, _converter, _get_base
from ._errors import RPackageMissing

# Conditionally import rpy2 components based on the active backend
if ACTIVE_BACKEND == "rpy2":
    from rpy2.robjects.conversion import localconverter
else:
    localconverter = None


def is_installed(package: str) -> bool:
    """Check if an R package is installed."""
    with localconverter(_converter):
        result = _get_base().requireNamespace(package, quietly=True)
        return bool(result[0])


def ensure_installed(package: str) -> None:
    """Ensure an R package is installed, raising RPackageMissing if not."""
    if not is_installed(package):
        raise RPackageMissing(package)
