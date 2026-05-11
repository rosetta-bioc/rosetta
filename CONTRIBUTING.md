# Contributing to rosetta

## Getting Started

```bash
git clone https://github.com/rosetta-bioc/rosetta.git
cd rosetta
pip install -e ".[dev]"
```

Requires R 4.0+ with BiocManager installed.

## Adding a New Wrapper

1. Create `rosetta/wrappers/<package>.py` following the pattern in `SPEC.md`
2. Accept pandas DataFrames, return pandas DataFrames
3. Call `ensure_installed("<R_PACKAGE>")` before any R calls
4. Expose R parameters as Python keyword arguments
5. Add the public function to `rosetta/__init__.py`

## Debugging Best Practices

### Understand Before You Fix

- Read error messages completely — especially R stack traces surfaced through rpy2
- Reproduce the bug consistently with a minimal example
- Identify the root cause, not just symptoms

### Systematic Investigation

- Start with logs and R stack traces
- Use `rpy2`'s logging to see raw R calls
- Test hypotheses with minimal changes — don't randomly change code hoping it works

### Fix the Root Cause

```python
# ❌ BAD: Hiding the error
try:
    r_result = run_deseq(counts)
except Exception:
    pass  # It'll probably be fine

# ✅ GOOD: Handle appropriately
try:
    r_result = run_deseq(counts)
except RRuntimeError as e:
    raise RDataError(f"DESeq2 failed: {e}") from e
```

```python
# ❌ BAD: Working around the symptom
if result["padj"].isna().all():
    result["padj"] = result["pvalue"]  # Why are they all NaN?

# ✅ GOOD: Fix the root cause — check input data
def _validate_counts(counts: pd.DataFrame):
    if (counts < 0).any().any():
        raise RDataError("Count matrix contains negative values")
```

### Don't Re-Implement Existing Solutions

Before implementing anything, check if it already exists:

```bash
git log --all --grep="feature-name"
git grep "pattern"
```

Warning signs: "this feels familiar", multiple files with similar names, duplicate functions.

### Code Archaeology for Legacy R Packages

When wrapping an unfamiliar R package:

1. **Inventory** — read the R package vignette and identify the core pipeline
2. **Map dependencies** — understand which R functions call which
3. **Trace the happy path** — run the R example in an R console first
4. **Then wrap** — only wrap what you understand

### Clean Up

- Remove temporary debugging code before committing
- Don't leave `print()` or `console.log` statements
- Don't commit commented-out code
- One fix per commit

### Validate

- Test the specific bug scenario
- Run the full test suite: `pytest`
- Check for unintended side effects in type conversion

## Bug Report Template

```markdown
## Bug Description
Brief description

## Steps to Reproduce
1. Input data characteristics (shape, types)
2. Function called with what arguments
3. Expected vs. actual result

## Root Cause
Why this happened

## Solution
What was changed and why
```

## Red Flags

Watch out for these in PRs:
- Catching and ignoring exceptions
- Band-aid fixes without understanding the cause
- "I'm not sure why this works, but it does"
- Skipping tests
- Implementing something that already exists in the codebase
