"""Environment detection, backend selection, and R package version pinning."""

import subprocess
import shutil
import json
import warnings
from typing import Optional

# Pin key R package versions to ensure scientific reproducibility in Week 9 validation tests
PINNED_VERSIONS = {
    "DESeq2": "1.40.0",
    "edgeR": "3.40.0",
    "limma": "3.54.0"
}

def check_rpy2_available() -> bool:
    """Check if rpy2 is installed and can successfully initialize the R runtime."""
    try:
        import importlib.util
        if importlib.util.find_spec("rpy2") is None:
            return False
        
        import rpy2.robjects as ro
        ro.r("1 + 1")
        return True
    except Exception:
        return False

def find_rscript() -> Optional[str]:
    """Locate the Rscript executable in the system path, returning None if not found."""
    rscript_path = shutil.which("Rscript")
    return rscript_path

def verify_r_package_versions(rscript_path: str):
    """
    Verify R package versions against PINNED_VERSIONS (Version Pinning).
    Crucial for reproducibility in Week 9 validation tests.
    """
    check_script = """
    pkgs <- c("DESeq2", "edgeR", "limma")
    versions <- sapply(pkgs, function(p) {
        if (requireNamespace(p, quietly = TRUE)) {
            as.character(packageVersion(p))
        } else {
            NA
        }
    })
    cat(jsonlite::toJSON(versions, auto_unbox = TRUE))
    """
    try:
        result = subprocess.run(
            [rscript_path, "-e", check_script],
            capture_output=True,
            text=True,
            check=True
        )
        installed_versions = json.loads(result.stdout.strip())
        for pkg, expected in PINNED_VERSIONS.items():
            actual = installed_versions.get(pkg)
            if actual is None:
                warnings.warn(f"Required R package '{pkg}' is not installed in the R environment.")
            elif actual != expected:
                warnings.warn(
                    f"R package version mismatch for '{pkg}': "
                    f"expected {expected}, found {actual}. Results may vary."
                )
    except Exception as e:
        warnings.warn(f"Could not verify R package versions: {e}")

def detect_backend() -> str:
    """
    Auto-detection logic for Week 8:
    Returns 'rpy2' if available and stable, otherwise falls back to 'subprocess'.
    """
    if check_rpy2_available():
        return "rpy2"
    
    # Check if Rscript is present; if not, fallback gracefully instead of crashing
    rscript_path = find_rscript()
    if rscript_path is None:
        # if there is no R in Python environment, return sth to install the module safely
        return "fallback"
        
    return "subprocess"