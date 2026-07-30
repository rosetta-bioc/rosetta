"""Tests for rosetta.__init__ module imports."""

import pytest


def test_main_imports():
    """Test that all main functions can be imported from rosetta."""
    import rosetta as rb
    
    # Tier 3 & Wrappers
    assert hasattr(rb, 'EdgeR')
    assert hasattr(rb, 'Limma')
    assert hasattr(rb, 'ClusterProfiler')
    assert hasattr(rb, 'vst')
    assert hasattr(rb, 'rlog')
    assert hasattr(rb, 'tmm_normalize')
    assert hasattr(rb, 'enrich_go')
    assert hasattr(rb, 'enrich_kegg')

    # Tier 2 (Class-based)
    assert hasattr(rb, 'Seurat')
    assert hasattr(rb, 'Phyloseq')
    assert hasattr(rb, 'VCF')

    # Tier 1 (Quick API)
    assert hasattr(rb, 'quick_seurat')
    assert hasattr(rb, 'quick_phyloseq')
    assert hasattr(rb, 'quick_locate_variants')
    assert hasattr(rb, 'quick_predict_coding')
    assert hasattr(rb, 'quick_deseq2')
    assert hasattr(rb, 'quick_edger')

    # Backward-compat aliases
    assert hasattr(rb, 'phyloseq')
    assert hasattr(rb, 'phyloseq_richness')
    assert hasattr(rb, 'seurat')
    
    # Exception classes
    assert hasattr(rb, 'RDataError')
    assert hasattr(rb, 'RFormulaError') 
    assert hasattr(rb, 'RPackageMissing')


def test_all_attribute():
    """Test that __all__ contains expected items."""
    import rosetta
    
    expected = [
        # Metadata
        "__version__",
        # Wrappers & Classes
        "DESeq2", "run_deseq2",
        "EdgeR",
        "Limma",
        "vst", "rlog", "tmm_normalize",
        "ClusterProfiler",
        "enrich_go", "enrich_kegg",
        # Tier 2
        "Seurat", "Phyloseq", "VCF",
        # Tier 1
        "quick_seurat", "quick_phyloseq", "quick_deseq2", "quick_edger",
        "quick_locate_variants", "quick_predict_coding",
        # Backward-compat aliases
        "phyloseq", "phyloseq_richness", "seurat",
        # Utilities
        "pipelines", "codegen", "plots",
        "RosettaDataFrame",
        "QuickResult",
        # Errors
        "RDataError", "RFormulaError", "RPackageMissing",
    ]
    
    assert hasattr(rosetta, '__all__')
    assert set(rosetta.__all__) == set(expected)


def test_function_callability():
    """Test that imported functions are callable."""
    import rosetta as rb
    
    assert callable(rb.run_deseq2)
    assert callable(rb.enrich_go)
    assert callable(rb.enrich_kegg)
    assert callable(rb.quick_seurat)
    assert callable(rb.quick_phyloseq)
    assert callable(rb.quick_locate_variants)
    assert callable(rb.quick_predict_coding)
    assert callable(rb.phyloseq)
    assert callable(rb.phyloseq_richness)
    assert callable(rb.seurat)
    assert callable(rb.VCF)


def test_exception_inheritance():
    """Test that imported exceptions have correct inheritance."""
    import rosetta as rb
    from rosetta._errors import RosettaError
    
    # Should be able to instantiate exceptions
    assert issubclass(rb.RDataError, RosettaError)
    assert issubclass(rb.RFormulaError, RosettaError)
    assert issubclass(rb.RPackageMissing, RosettaError)