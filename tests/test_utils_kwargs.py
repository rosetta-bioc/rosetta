import pytest
import logging
from rosetta.utils.kwargs import filter_kwargs

def test_filter_kwargs_basic(caplog):
    """Test that allowed parameters are kept and disallowed are filtered with warning."""
    allowed = {"n_neighbors", "metric"}
    input_kwargs = {
        "n_neighbors": 15, 
        "metric": "cosine", 
        "invalid_param": "should_be_removed"
    }
    
    with caplog.at_level(logging.WARNING):
        result = filter_kwargs(input_kwargs, allowed)
        
    # Verify filtering results
    assert "n_neighbors" in result
    assert "metric" in result
    assert "invalid_param" not in result
    assert result["n_neighbors"] == 15
    
    # Verify warning message content
    assert "Parameter 'invalid_param' is not supported" in caplog.text

def test_filter_kwargs_type_conversion():
    """Test that Python types are correctly handled for R types."""
    allowed = {"verbose", "assay"}
    input_kwargs = {
        "verbose": True, 
        "assay": None
    }
    
    result = filter_kwargs(input_kwargs, allowed)
    
    assert result["verbose"] is True
    # None is preserved as-is (or converted to R NULL if rpy2 is available);
    # either way the key must be present in the filtered result.
    assert "assay" in result

def test_empty_kwargs():
    """Test handling of empty input dictionary."""
    assert filter_kwargs({}, {"param1"}) == {}