# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Tier 1 Quick API: `quick_seurat`, `quick_phyloseq` (PR #21)
- Backward-compat aliases: `rb.phyloseq`, `rb.seurat`, `rb.phyloseq_richness`
- Three-Tier API organization in `__init__.py`

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
