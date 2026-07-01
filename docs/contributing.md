# Contributing

## Development setup

```bash
git clone https://github.com/rosetta-bioc/rosetta.git
cd rosetta
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Requires R 4.0+ with BiocManager installed. See [R Packages](setup/r-packages.md) for which R packages are needed for testing.

---

## Running tests

```bash
# Full test suite
pytest

# Specific module
pytest tests/wrappers/test_deseq2.py

# With coverage
pytest --cov=rosetta
```

Tests that require R packages will be skipped if those packages aren't installed.

---

## Code style

We use [ruff](https://docs.astral.sh/ruff/) for linting and formatting:

```bash
ruff check .
ruff format .
```

Key conventions:

- Type hints on public API functions
- Docstrings for all public functions and classes (Google style)
- No `print()` statements in library code (use `codegen._emit()` for R code output)

---

## Adding a new wrapper

Follow the pattern in `SPEC.md`:

1. **Create** `rosetta/wrappers/<package>.py`
2. **Accept** pandas DataFrames as input
3. **Return** pandas DataFrames (or `RosettaDataFrame`) as output
4. **Call** `ensure_installed("<R_PACKAGE>")` before any R calls
5. **Validate** inputs in Python before crossing the R boundary
6. **Expose** R parameters as Python keyword arguments
7. **Add** the public function to `rosetta/__init__.py` and `__all__`
8. **Write tests** in `tests/wrappers/test_<package>.py`

### Template

```python
"""<package> wrapper."""

import pandas as pd
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from .._bridge import _converter, to_r_matrix, to_pandas, to_r_df
from .._deps import ensure_installed
from .._errors import RDataError


def my_function(counts: pd.DataFrame, metadata: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """Brief description.

    Args:
        counts: Gene count matrix (genes x samples).
        metadata: Sample metadata DataFrame.
        **kwargs: Passed to R function.

    Returns:
        DataFrame with results.
    """
    if (counts < 0).any().any():
        raise RDataError("Count matrix contains negative values")

    ensure_installed("PackageName")
    pkg = importr("PackageName")

    r_counts = to_r_matrix(counts)

    with localconverter(_converter):
        result = pkg.some_function(r_counts, **kwargs)

    return to_pandas(to_r_df(result))
```

---

## Debugging best practices

### Understand before you fix

- Read error messages completely — especially R stack traces surfaced through rpy2
- Reproduce the bug with a minimal example
- Identify the root cause, not just symptoms

### Systematic investigation

- Start with logs and R stack traces
- Use `rpy2`'s logging to see raw R calls
- Test hypotheses with minimal changes

### Fix the root cause

```python
# ❌ BAD: Hiding the error
try:
    r_result = run_deseq(counts)
except Exception:
    pass

# ✅ GOOD: Handle appropriately
try:
    r_result = run_deseq(counts)
except RRuntimeError as e:
    raise RDataError(f"DESeq2 failed: {e}") from e
```

---

## PR guidelines

1. **One feature/fix per PR** — keep changes focused
2. **Write tests** for new functionality
3. **Run the full test suite** before submitting
4. **Update docstrings** if you change function signatures
5. **Don't commit** commented-out code, `print()` statements, or debugging artifacts

### PR description template

```markdown
## Summary
Brief description of the change.

## Changes
- What was added/modified/removed

## Testing
- How was this tested?
- Any manual testing steps?

## Checklist
- [ ] Tests pass (`pytest`)
- [ ] Linting passes (`ruff check .`)
- [ ] Docstrings updated
- [ ] CHANGELOG.md updated (for user-facing changes)
```

---

## Red flags in code review

Watch out for:

- Catching and ignoring exceptions
- Band-aid fixes without understanding the cause
- "I'm not sure why this works, but it does"
- Implementing something that already exists in the codebase
- Skipping tests

---

## Project structure

```
rosetta/
├── __init__.py          # Public API, exports, Quick API functions
├── _bridge.py           # rpy2 type conversion utilities
├── _deps.py             # R package availability checking
├── _errors.py           # Custom exception classes
├── codegen.py           # R code logging/printing
├── results.py           # RosettaDataFrame with .report()
├── quick_result.py      # QuickResult for Tier 1 dict-like results
├── pipelines.py         # Pipeline utilities
├── wrappers/
│   ├── deseq2.py        # DESeq2 wrapper
│   ├── edger.py         # edgeR wrapper
│   ├── limma.py         # limma-voom wrapper
│   ├── clusterprofiler.py  # ORA + GSEA wrappers
│   ├── seurat.py        # Seurat class wrapper
│   └── phyloseq.py      # phyloseq class wrapper
├── stats/               # Statistical utilities (design matrices, TREAT)
└── plots/               # Visualization (volcano, MA, PCA)
```

---

## Getting help

- **Issues:** [GitHub Issues](https://github.com/rosetta-bioc/rosetta/issues)
- **Good first issues:** labeled in the issue tracker
- **Discussions:** open an issue for design questions before implementing large changes
