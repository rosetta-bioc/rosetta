"""R session management and bidirectional type conversion."""

import os
import contextlib
import numpy as np
import pandas as pd

# Suppress rpy2 ABI mode warning during import
with open(os.devnull, "w") as _devnull, contextlib.redirect_stderr(_devnull):
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import Converter, localconverter
    from rpy2.robjects.packages import importr

_converter = Converter("rosetta")
_converter += numpy2ri.converter
_converter += pandas2ri.converter
_converter += ro.default_converter

_base = None


def _get_base():
    """Lazily import R base package."""
    global _base
    if _base is None:
        _base = importr("base")
    return _base


def to_r_matrix(df: pd.DataFrame):
    """Convert pandas DataFrame to R matrix."""
    from ._errors import RDataError
    if not isinstance(df, pd.DataFrame):
        raise RDataError("Expected pandas DataFrame")
    with localconverter(_converter):
        return _get_base().as_matrix(ro.conversion.get_conversion().py2rpy(df))


def to_r_dataframe(df: pd.DataFrame):
    """Convert pandas DataFrame to R data.frame."""
    from ._errors import RDataError
    if not isinstance(df, pd.DataFrame):
        raise RDataError("Expected pandas DataFrame")
    with localconverter(_converter):
        return ro.conversion.get_conversion().py2rpy(df)


def to_pandas(r_obj) -> "pd.DataFrame":
    """Convert R data.frame/matrix to pandas DataFrame (with .report() method)."""
    from .results import RosettaDataFrame
    with localconverter(_converter):
        df = ro.conversion.get_conversion().rpy2py(r_obj)
    if isinstance(df, pd.DataFrame):
        return RosettaDataFrame(df)
    return df


def to_r_df(r_obj):
    """Convert an R object to R data.frame via base::as.data.frame."""
    with localconverter(_converter):
        return _get_base().as_data_frame(r_obj)


def r_nrow(r_obj):
    """Get nrow of an R object via base::nrow."""
    with localconverter(_converter):
        result = _get_base().nrow(r_obj)
        return int(result[0])  # Convert R vector to Python int
