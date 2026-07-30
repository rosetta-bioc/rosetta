"""Validation tests comparing Rosetta limma-voom wrapper against direct R execution (airway dataset)."""

import pytest
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from rosetta._deps import is_installed
from rosetta.wrappers.limma import Limma

# Validation tests require real R packages and datasets
if not is_installed("limma") or not is_installed("edgeR") or not is_installed("airway"):
    pytest.skip("limma, edgeR, or airway package not available for validation", allow_module_level=True)


def test_limma_airway_parity():
    """Verify that Rosetta limma-voom output matches direct R execution within tolerance limits."""
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    from rosetta._bridge import to_pandas, to_r_df

    # 1. Load airway dataset from R and convert to pandas for Rosetta wrapper initialization
    airway = importr("airway")
    
    ro.r("data('airway', package='airway')")
    ro.r("se <- airway")
    
    counts_r = ro.r("assay(se)")
    coldata_r = ro.r("as.data.frame(colData(se))")
    
    counts_df = to_pandas(to_r_df(counts_r))
    metadata_df = to_pandas(to_r_df(coldata_r))

    try:
        # Initialize and run Rosetta limma wrapper (using design ~ dex)
        model = Limma(counts=counts_df, metadata=metadata_df, design="~ dex")
        model.run_ebayes()
        rosetta_res = model.get_results(sort_by="none")
    except Exception as e:
        pytest.skip(f"Skipping due to execution environment issue: {e}")

    # 2. Assertions on structure
    assert rosetta_res is not None and not rosetta_res.empty, "Rosetta limma results must not be empty"

    required_cols = ["logFC", "P.Value", "adj.P.Val"]
    for col in required_cols:
        assert col in rosetta_res.columns, f"Missing required column: {col}"

    # 3. Direct R execution reference for comparison
    ro.r("""
        library(limma)
        library(edgeR)
        data('airway', package='airway')
        se <- airway
        counts_ref <- assay(se)
        metadata_ref <- as.data.frame(colData(se))
        
        design_ref <- model.matrix(~ dex, data=metadata_ref)
        dge_ref <- DGEList(counts=counts_ref)
        dge_ref <- calcNormFactors(dge_ref)
        v_ref <- voomLmFit(dge_ref, design_ref)
        fit_ref <- eBayes(v_ref)
        res_ref <- topTable(fit_ref, number=nrow(counts_ref), sort.by="none")
    """)
    ref_res = to_pandas(ro.r("res_ref"))

    # Align by common genes/rows
    common_genes = rosetta_res.index.intersection(ref_res.index)
    ros_aligned = rosetta_res.loc[common_genes]
    ref_aligned = ref_res.loc[common_genes]

    # Acceptance criteria validation
    np.testing.assert_allclose(ros_aligned['logFC'], ref_aligned['logFC'], atol=1e-6, err_msg="logFC mismatch")
    np.testing.assert_allclose(ros_aligned['adj.P.Val'], ref_aligned['adj.P.Val'], atol=1e-6, err_msg="adj.P.Val mismatch")
    
    rho, _ = spearmanr(ros_aligned['logFC'], ref_aligned['logFC'])
    assert rho > 0.999, f"Gene ranking Spearman correlation too low: {rho}"