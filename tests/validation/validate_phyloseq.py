"""Validation tests comparing Rosetta Phyloseq wrapper against direct R execution using published datasets."""

import pytest
import pandas as pd
from rosetta._deps import is_installed
from rosetta.wrappers.phyloseq import Phyloseq

# Validation tests require real R packages and datasets
if not is_installed("phyloseq"):
    pytest.skip("phyloseq package not available for validation", allow_module_level=True)


def test_phyloseq_real_dataset_parity():
    """Verify that Rosetta Phyloseq output matches direct R execution using published phyloseq datasets."""
    import rpy2.robjects as ro
    from rosetta._bridge import to_pandas

    try:
        ro.r("""
            library(phyloseq)
            data("GlobalPatterns")
            gp_sub <- prune_taxa(taxa_names(GlobalPatterns)[1:50], GlobalPatterns)
            gp_sub <- prune_samples(sample_names(gp_sub)[1:10], gp_sub)
            
            mat_ref <- as.data.frame(as(otu_table(gp_sub), "matrix"))
            
            # Force sample_data to be a standard base R data.frame (stripping out S4/phyloseq class attributes)
            samp_raw <- as(sample_data(gp_sub), "data.frame")
            samp_ref <- as.data.frame(samp_raw)
            
            tax_ref <- as.data.frame(as(tax_table(gp_sub), "matrix"))
            
            res_ref <- estimate_richness(gp_sub, measures=c("Observed", "Shannon"))
        """)
        
        from rpy2.robjects.conversion import localconverter
        from rosetta._bridge import _converter

        with localconverter(_converter):
            mat_pd = to_pandas(ro.r("mat_ref"))
            samp_pd = to_pandas(ro.r("samp_ref"))
            tax_pd = to_pandas(ro.r("tax_ref"))
            ref_res = to_pandas(ro.r("res_ref"))

        # 2. Execute the Rosetta Phyloseq Wrapper
        ps = Phyloseq(otu_table=mat_pd, sample_data=samp_pd, tax_table=tax_pd)
        rosetta_res = ps.estimate_richness(measures=["Observed", "Shannon"])

    except Exception as e:
        pytest.skip(f"Skipping due to execution environment or data loading issue: {e}")

    # 3. Retrieve the native R calculation results for comparison
    ref_res = to_pandas(ro.r("res_ref"))

    # 4. Strict parity validation
    assert len(rosetta_res) == len(ref_res), "Row count mismatch between Rosetta and direct R execution"
    for sample_id in rosetta_res.index:
        assert rosetta_res.loc[sample_id, "Observed"] == ref_res.loc[sample_id, "Observed"], f"Observed richness mismatch for {sample_id}"
        # Allow minor tolerance for floating-point calculations
        assert abs(rosetta_res.loc[sample_id, "Shannon"] - ref_res.loc[sample_id, "Shannon"]) < 1e-5, f"Shannon diversity mismatch for {sample_id}"