"""
DESeq2 example — airway dataset (Himes et al. 2014)

Demonstrates the full rosetta DESeq2 pipeline:
  pandas counts matrix → R DESeq2 → pandas results

Dataset: 8 bronchial epithelial cell samples
         4 untreated + 4 dexamethasone-treated
         ~33,000 genes

Usage:
    python examples/statistics/deseq2_airway.py

Requires: R 4.0+, DESeq2, airway (Bioconductor)
"""

import time
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

import rosetta as rb


def load_airway():
    """Load the airway dataset from Bioconductor via rpy2."""
    from rosetta._deps import ensure_installed
    ensure_installed("airway")

    ro.r('library(airway)')
    ro.r('data(airway)')
    ro.r('airway <- as(airway, "DESeqDataSet")')

    # Extract count matrix
    counts_r = ro.r('as.data.frame(counts(airway))')
    metadata_r = ro.r('as.data.frame(colData(airway)[, "dex", drop=FALSE])')

    from rpy2.robjects import pandas2ri
    with (ro.default_converter + pandas2ri.converter).context():
        counts = ro.conversion.rpy2py(counts_r)
        metadata = ro.conversion.rpy2py(metadata_r)

    metadata.columns = ["condition"]
    metadata["condition"] = metadata["condition"].astype(str)

    return counts, metadata


def main():
    print("rosetta DESeq2 example — airway dataset")
    print("=" * 45)

    print("\n[1/4] Loading airway dataset...")
    counts, metadata = load_airway()
    print(f"      {counts.shape[0]:,} genes × {counts.shape[1]} samples")
    print(f"      Conditions: {metadata['condition'].value_counts().to_dict()}")

    print("\n[2/4] Fitting DESeq2 model...")
    t0 = time.time()
    dds = rb.run_deseq2(counts, metadata, design="~ condition")
    print(f"      Done in {time.time() - t0:.1f}s")

    print("\n[3/4] Extracting results (treated vs untreated)...")
    results = rb.get_results(
        dds,
        contrast=["condition", "trt", "untrt"],
        alpha=0.05,
    )
    print(f"      {len(results):,} genes tested")

    sig = results.dropna(subset=["padj"])
    sig = sig[sig["padj"] < 0.05]
    up = sig[sig["log2FoldChange"] > 0]
    down = sig[sig["log2FoldChange"] < 0]
    print(f"      {len(sig)} significant DEGs (padj < 0.05)")
    print(f"      {len(up)} up-regulated, {len(down)} down-regulated")

    print("\n[4/4] Top 10 DEGs by adjusted p-value:")
    top = results.dropna(subset=["padj"]).sort_values("padj").head(10)
    print(top[["baseMean", "log2FoldChange", "padj"]].to_string())

    print("\n✓ Complete. Results are a standard pandas DataFrame — plug into")
    print("  any Python ML pipeline, visualization library, or nodes.bio.")


if __name__ == "__main__":
    main()
