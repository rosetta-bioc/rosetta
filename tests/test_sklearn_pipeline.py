"""Tests for scikit-learn compatible pipeline integration using airway dataset."""

import pytest
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier

from rosetta.sklearn_pipeline import DESeq2Transformer
from rosetta._deps import is_installed
from rosetta._bridge import to_pandas, to_r_df


@pytest.mark.skipif(not is_installed("airway") or not is_installed("DESeq2"), reason="Required R packages (airway, DESeq2) not installed")
def test_airway_sklearn_pipeline():
    """Test end-to-end integration of DESeq2Transformer with StandardScaler, PCA, and RandomForest."""
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr

    # Load airway dataset from R
    airway_pkg = importr("airway")
    ro.r("data(airway)")
    
    # Extract assay (counts) and colData (metadata)
    counts_r = ro.r("assay(airway)")
    coldata_r = ro.r("colData(airway)")
    
    # Standard Python DataFrame conversion using Rosetta's bridge
    counts_df = to_pandas(to_r_df(counts_r))
    metadata = to_pandas(to_r_df(coldata_r))
    metadata["condition"] = metadata["dex"]

    # Standard sklearn format: X is (samples x genes)
    X = counts_df.T  
    y = metadata["dex"].values

    # Build scikit-learn Pipeline
    pipe = Pipeline([
        ('deseq2', DESeq2Transformer(metadata=metadata, design="~ dex", alpha=0.05)),
        ('scale', StandardScaler()),
        ('pca', PCA(n_components=2)),
        ('classify', RandomForestClassifier(random_state=42))
    ])

    # Fit and predict just like a normal sklearn user
    pipe.fit(X, y)
    predictions = pipe.predict(X)
    
    # Assertions
    assert len(predictions) == len(y)
    assert "deseq2" in pipe.named_steps
    assert "pca" in pipe.named_steps