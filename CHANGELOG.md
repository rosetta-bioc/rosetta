# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.2] - 2026-08-11

### Added
- `rosetta.sklearn_compat`: `DESeq2Transformer`, `EdgeRTransformer`, `LimmaTransformer` — sklearn `TransformerMixin` wrappers enabling use inside `sklearn.pipeline.Pipeline`
- `RosettaDataFrame.volcano()` — volcano plot helper (matplotlib, optional dep)
- `RosettaDataFrame.ma_plot()` — MA plot helper with significant-gene highlighting
- `RosettaDataFrame._rosetta_method` metadata tracking; `.report()` now includes method header
- `subprocess` fallback backend with auto-detection and Bioconductor version pinning
- Granular `**kwargs` passthrough to underlying R functions
- R escape-hatch mode for environments without rpy2

### Fixed
- `results.py`: `from __future__ import annotations` so `str | None` union syntax works on Python 3.9
- `test_detect.py`: fixed `test_check_rpy2_available_false_on_import_error` patch target for Python 3.12 compatibility (was patching `rpy2.robjects.r` — now patches `builtins.__import__`)
- `_converter` and `_bridge.py` rpy2 import guards for pure-Python fallback
- Bool/float/str R vector coercion in `filter`
- `lfc_shrink()` pipeline integration

## [0.3.0] - 2026-07-01

### Added
- Tier 1 Quick API: `quick_deseq2`, `quick_edger`, `quick_seurat`, `quick_phyloseq`
- `QuickResult` class with `.report()` for all Quick API functions
- Plotting module: `volcano()`, `ma_plot()`, `pca()` with auto-detection of DESeq2/edgeR/limma columns
- Normalization wrappers: `vst()`, `rlog()`, `tmm_normalize()`
- VariantAnnotation module: `read_vcf`, `locate_variants`, `predict_coding`, `VCF` class
- `__version__` attribute via importlib.metadata
- CHANGELOG.md
- CI test workflow (Python 3.9/3.12, R 4.4, cached Bioconductor packages)
- mkdocs documentation site with 14 pages (API reference, guides, setup)
- Distribution content (Biostars tutorial, SO answers, blog post, Twitter thread)
- Backward-compat aliases: `rb.phyloseq`, `rb.seurat`, `rb.phyloseq_richness`

### Fixed
- `_phyloseq_available()` docstring and case-sensitive package name
- `test_init.py` updated for Three-Tier `__all__` exports
- README test count (was "170+", now accurate)

## [0.2.2] - 2026-06-26

### Added
- Seurat class-based wrapper with `run_standard_pipeline()`, `find_markers()`
- Phyloseq class-based wrapper with `estimate_richness()`, `run_ordination()`
- zizmor GitHub Actions security scanning

### Fixed
- `__all__` exports updated for new class names
- `find_markers()` R parameter convention (ident.1/ident.2)
- `Idents<-` setter through rpy2

## [0.2.1] - 2026-06-18

### Added
- ORA and GSEA class wrappers for clusterProfiler
- `_run_r_enrichment()` facade pattern
- `enrichGO`, `enrichKEGG`, `gseGO`, `gseKEGG` integration

## [0.2.0] - 2026-06-09

### Added
- `rosetta/stats/` module: `treat.py`, `design.py`, `decide.py`
- edgeR `glmTreat()` quasi-likelihood support
- Test reorganization into `tests/wrappers/` and `tests/stats/`

## [0.1.0] - 2026-06-03

### Added
- Modular DESeq2 API: `run_deseq2`, `get_results`, `lfc_shrink`
- `preview_design()`, `get_results_names()`
- Docker environment with Bioconductor dependencies
- Initial test suite
