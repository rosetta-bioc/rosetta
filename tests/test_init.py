"""Tests for rosetta.__init__ module imports."""

import pytest


def test_main_imports():
    """Test that all main functions can be imported from rosetta."""
    import rosetta as rb
    
    # Tier 3 (Functional/Legacy)
    assert hasattr(rb, 'deseq2')
    assert hasattr(rb, 'edger') 
    assert hasattr(rb, 'limma_voom')
    assert hasattr(rb, 'enrich_go')
    assert hasattr(rb, 'enrich_kegg')
    assert hasattr(rb, 'enrich_pathway')
    assert hasattr(rb, 'enrich_custom')

    # Tier 2 (Class-based)
    assert hasattr(rb, 'Seurat')
    assert hasattr(rb, 'Phyloseq')

    # Tier 1 (Quick API)
    assert hasattr(rb, 'quick_seurat')
    assert hasattr(rb, 'quick_phyloseq')

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
        # Tier 3
        "deseq2", "edger", "limma_voom",
        "ORA", "GSEA", "enrichment",
        "enrich_go", "enrich_kegg", "enrich_pathway", "enrich_custom",
        # Tier 2
        "Seurat", "Phyloseq",
        # Tier 1
        "quick_seurat", "quick_phyloseq", "quick_deseq2", "quick_edger",
        # Backward-compat aliases
        "phyloseq", "phyloseq_richness", "seurat",
        # Utilities
        "pipelines", "codegen",
        "RosettaDataFrame",
        # Errors
        "RDataError", "RFormulaError", "RPackageMissing",
    ]
    
    assert hasattr(rosetta, '__all__')
    assert set(rosetta.__all__) == set(expected)


def test_function_callability():
    """Test that imported functions are callable."""
    import rosetta as rb
    
    assert callable(rb.deseq2)
    assert callable(rb.edger)
    assert callable(rb.limma_voom)
    assert callable(rb.enrich_go)
    assert callable(rb.enrich_kegg)
    assert callable(rb.enrich_pathway)
    assert callable(rb.enrich_custom)
    assert callable(rb.quick_seurat)
    assert callable(rb.quick_phyloseq)
    assert callable(rb.phyloseq)
    assert callable(rb.phyloseq_richness)
    assert callable(rb.seurat)


def test_exception_inheritance():
    """Test that imported exceptions have correct inheritance."""
    import rosetta as rb
    from rosetta._errors import RosettaError
    
    # Should be able to instantiate exceptions
    assert issubclass(rb.RDataError, RosettaError)
    assert issubclass(rb.RFormulaError, RosettaError)
    assert issubclass(rb.RPackageMissing, RosettaError)
