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


# --- list → R vector coercion (requires rpy2) ---

def test_filter_kwargs_bool_list(monkeypatch):
    """bool list → ro.BoolVector."""
    import types
    fake_ro = types.SimpleNamespace(
        BoolVector=lambda v: ("BoolVector", v),
        IntVector=lambda v: ("IntVector", v),
        FloatVector=lambda v: ("FloatVector", v),
        StrVector=lambda v: ("StrVector", v),
    )
    fake_ro.r = lambda _: None
    monkeypatch.setattr("rosetta.utils.kwargs.ro", fake_ro)

    result = filter_kwargs({"flags": [True, False, True]}, {"flags"})
    assert result["flags"] == ("BoolVector", [True, False, True])


def test_filter_kwargs_int_list(monkeypatch):
    """int list → ro.IntVector."""
    import types
    fake_ro = types.SimpleNamespace(
        BoolVector=lambda v: ("BoolVector", v),
        IntVector=lambda v: ("IntVector", v),
        FloatVector=lambda v: ("FloatVector", v),
        StrVector=lambda v: ("StrVector", v),
    )
    fake_ro.r = lambda _: None
    monkeypatch.setattr("rosetta.utils.kwargs.ro", fake_ro)

    result = filter_kwargs({"dims": [1, 2, 3]}, {"dims"})
    assert result["dims"] == ("IntVector", [1, 2, 3])


def test_filter_kwargs_float_list(monkeypatch):
    """float list → ro.FloatVector."""
    import types
    fake_ro = types.SimpleNamespace(
        BoolVector=lambda v: ("BoolVector", v),
        IntVector=lambda v: ("IntVector", v),
        FloatVector=lambda v: ("FloatVector", v),
        StrVector=lambda v: ("StrVector", v),
    )
    fake_ro.r = lambda _: None
    monkeypatch.setattr("rosetta.utils.kwargs.ro", fake_ro)

    result = filter_kwargs({"weights": [0.1, 0.9]}, {"weights"})
    assert result["weights"] == ("FloatVector", [0.1, 0.9])


def test_filter_kwargs_str_list(monkeypatch):
    """str list → ro.StrVector."""
    import types
    fake_ro = types.SimpleNamespace(
        BoolVector=lambda v: ("BoolVector", v),
        IntVector=lambda v: ("IntVector", v),
        FloatVector=lambda v: ("FloatVector", v),
        StrVector=lambda v: ("StrVector", v),
    )
    fake_ro.r = lambda _: None
    monkeypatch.setattr("rosetta.utils.kwargs.ro", fake_ro)

    result = filter_kwargs({"genes": ["BRCA1", "TP53"]}, {"genes"})
    assert result["genes"] == ("StrVector", ["BRCA1", "TP53"])


def test_filter_kwargs_mixed_list_falls_back_to_strVector(monkeypatch):
    """mixed-type list → ro.StrVector of str(x)."""
    import types
    fake_ro = types.SimpleNamespace(
        BoolVector=lambda v: ("BoolVector", v),
        IntVector=lambda v: ("IntVector", v),
        FloatVector=lambda v: ("FloatVector", v),
        StrVector=lambda v: ("StrVector", v),
    )
    fake_ro.r = lambda _: None
    monkeypatch.setattr("rosetta.utils.kwargs.ro", fake_ro)

    result = filter_kwargs({"mixed": [1, "two", 3.0]}, {"mixed"})
    assert result["mixed"] == ("StrVector", ["1", "two", "3.0"])


def test_filter_kwargs_none_without_rpy2(monkeypatch):
    """None with ro=None preserved as Python None."""
    monkeypatch.setattr("rosetta.utils.kwargs.ro", None)
    result = filter_kwargs({"assay": None}, {"assay"})
    assert result["assay"] is None


def test_filter_kwargs_scalar_passthrough():
    """Non-list, non-bool, non-None scalars pass through unchanged."""
    result = filter_kwargs({"resolution": 0.5, "n_top": 100}, {"resolution", "n_top"})
    assert result["resolution"] == 0.5
    assert result["n_top"] == 100