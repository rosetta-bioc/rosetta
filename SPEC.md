# rosetta — Technical Specification

## Overview

`rosetta` (`rosetta-bioc`) is a Python library that wraps R/Bioconductor bioinformatics packages (DESeq2, edgeR, limma, Seurat, phyloseq, clusterProfiler, etc.) via `rpy2`, providing a pandas-native API, reproducible codegen, and rich `.report()` summaries.

## Architecture

```
rosetta/
├── plots/               # Visualization wrappers and plotting utilities
├── stats/               # Statistical helper functions
├── utils/               # Internal utility functions
├── wrappers/            # R package wrappers
│   ├── __init__.py
│   ├── deseq2.py        # DESeq2 wrapper
│   ├── edger.py         # EdgeR wrapper
│   ├── limma.py         # limma-voom wrapper
│   ├── seurat.py        # Seurat wrapper
│   ├── clusterprofiler.py # ORA and GSEA wrapper
│   ├── phyloseq.py      # Microbiome diversity and ordination wrapper
│   ├── normalize.py     # Normalization wrapper
│   ├── vcf.py           # VCF parsing wrapper
│   └── variant_annotation.py
├── __init__.py          # Public API (Tier 1 quick functions, Tier 2 classes, aliases)
├── __main__.py          # CLI entry point
├── _bridge.py           # rpy2 session management and R↔Python type conversion
├── _deps.py             # R/Bioconductor package detection and installation
├── _detect.py           # Environment and package detection utilities
├── _errors.py           # R error translation to Python exceptions
├── codegen.py           # R code generation and tracking mechanism
├── example.py           # Built-in usage examples
├── pipelines.py         # Multi-step analysis pipelines
├── quick_result.py      # QuickResult container with .report() formatting
├── results.py           # RosettaDataFrame subclass for pandas outputs with .report()
├── sklearn_compat.py    # scikit-learn compatibility layer

```

## Core Components

### 1. Bridge Layer (`_bridge.py`)

Manages a single `rpy2` R session and handles all type conversion between Python and R.

* Key conversions:
* `pandas.DataFrame` ↔ `R data.frame`
* `numpy.ndarray` ↔ `R matrix`
* Python `dict` ↔ `R named list`
* `None` ↔ `R NULL`



### 2. Three-Tier API Design

* **Tier 1 — Quick API (`quick_*()`)**: One-liner functions (e.g., `quick_deseq2`, `quick_edger`, `quick_seurat`, `quick_phyloseq`) optimized for rapid notebook analysis. Returns `RosettaDataFrame` or `QuickResult` with built-in `.report()` methods.
* **Tier 2 — Class-based API (`Seurat()`, `DESeq2()`, `Phyloseq()`)**: Stateful, chainable wrappers providing step-by-step pipeline control.
* **Tier 3 — Escape Hatch / Functional API**: Direct access to underlying R objects (`model.r_obj`) and execution of custom R snippets (`model.run_r_script()`), along with granular functional operations.

### 3. Dependency Manager (`_deps.py` & `_detect.py`)

* Checks if required R/Bioconductor packages are installed prior to execution.
* If missing, guides or triggers installation via `BiocManager::install()`.

### 4. Error Translation (`_errors.py`)

Catches `rpy2` runtime exceptions and maps them to descriptive Python errors:

* `RPackageMissing` — R package not installed
* `RFormulaError` — invalid R design formula
* `RDataError` — incompatible input data (e.g., negative counts)

### 5. Transparency & Codegen (`codegen.py`)

* Tracks and logs executed R commands behind the scenes.
* Allows users to inspect or export equivalent native R scripts (`rb.codegen.enable()`, `rb.codegen.last()`).

## Design Decisions

* **Wrap, don't reimplement** — statistical correctness relies entirely on original R/Bioconductor packages.
* **Pandas-native** — inputs are pandas DataFrames/Series, outputs are `RosettaDataFrame` subclasses supporting `.report()`.
* **Lazy R initialization** — R session starts on first wrapper call, not upon import.
* **Reproducibility first** — seamless translation between Python interface and R execution backend.

## Dependencies

* Python 3.9+
* `rpy2 >= 3.5`
* `pandas >= 1.5`
* `numpy >= 1.23`
* R 4.0+ with BiocManager