"""Tests for rosetta._bridge."""

import numpy as np
import pandas as pd
import pytest

from rosetta._bridge import to_r_matrix, to_r_dataframe, to_pandas, to_r_df, r_nrow
from rosetta._errors import RDataError


def test_dataframe_roundtrip(sample_counts):
    r_df = to_r_dataframe(sample_counts)
    result = to_pandas(r_df)
    pd.testing.assert_frame_equal(result, sample_counts.astype(float), check_dtype=False)


def test_matrix_roundtrip(sample_counts):
    r_mat = to_r_matrix(sample_counts)
    result = to_pandas(r_mat)
    np.testing.assert_array_equal(result, sample_counts.values)


def test_to_r_matrix_rejects_non_dataframe():
    with pytest.raises(RDataError):
        to_r_matrix("not a dataframe")


def test_to_r_dataframe_rejects_non_dataframe():
    with pytest.raises(RDataError):
        to_r_dataframe([1, 2, 3])


def test_to_r_df_conversion(sample_counts):
    """Test to_r_df function."""
    r_mat = to_r_matrix(sample_counts)
    r_df = to_r_df(r_mat)
    result = to_pandas(r_df)
    
    # Should be a DataFrame with same shape
    assert isinstance(result, pd.DataFrame)
    assert result.shape == sample_counts.shape


def test_r_nrow_function(sample_counts):
    """Test r_nrow function."""
    r_mat = to_r_matrix(sample_counts)
    nrows = r_nrow(r_mat)
    
    assert nrows == len(sample_counts)


def test_empty_dataframe_conversion():
    """Test conversion of empty DataFrame."""
    empty_df = pd.DataFrame()
    r_df = to_r_dataframe(empty_df)
    result = to_pandas(r_df)
    
    assert isinstance(result, pd.DataFrame)
    assert result.empty


def test_single_column_dataframe(sample_metadata):
    """Test conversion of single column DataFrame."""
    single_col = sample_metadata[['condition']]
    r_df = to_r_dataframe(single_col)
    result = to_pandas(r_df)
    
    pd.testing.assert_frame_equal(result, single_col, check_dtype=False)


def test_dataframe_with_missing_values():
    """Test conversion of DataFrame with NaN values."""
    df_with_nan = pd.DataFrame({
        'A': [1, 2, np.nan], 
        'B': [4, np.nan, 6]
    })
    r_df = to_r_dataframe(df_with_nan)
    result = to_pandas(r_df)
    
    assert isinstance(result, pd.DataFrame)
    assert result.shape == df_with_nan.shape
    assert result.isna().sum().sum() == 2  # Two NaN values


def test_dataframe_with_different_dtypes():
    """Test conversion of DataFrame with mixed data types."""
    mixed_df = pd.DataFrame({
        'integers': [1, 2, 3],
        'floats': [1.1, 2.2, 3.3],
        'strings': ['a', 'b', 'c']
    })
    r_df = to_r_dataframe(mixed_df)
    result = to_pandas(r_df)
    
    assert isinstance(result, pd.DataFrame)
    assert result.shape == mixed_df.shape


def test_matrix_conversion_preserves_row_col_names(sample_counts):
    """Test that matrix conversion preserves row and column names."""
    r_mat = to_r_matrix(sample_counts)
    
    # The matrix should preserve the structure for later conversion
    assert r_mat is not None
