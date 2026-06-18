"""rosetta quick-start example — synthetic data, no files needed.

Run: python -m rosetta.example
"""

import numpy as np
import pandas as pd


def main():
    """Demonstrate rosetta with synthetic RNA-seq count data."""
    from .results import RosettaDataFrame

    print("🪨 rosetta — quick demo with synthetic data\n")

    # Generate fake count matrix (1000 genes × 6 samples)
    np.random.seed(42)
    n_genes, n_samples = 1000, 6
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    sample_names = [f"S{i+1}" for i in range(n_samples)]

    # Base expression + condition effect for first 100 genes
    base = np.random.negative_binomial(n=5, p=0.01, size=(n_genes, n_samples))
    base[:100, 3:] += np.random.negative_binomial(n=3, p=0.01, size=(100, 3))  # upregulated in treated

    counts = pd.DataFrame(base, index=gene_names, columns=sample_names)
    metadata = pd.DataFrame(
        {"condition": ["control"] * 3 + ["treated"] * 3},
        index=sample_names,
    )

    print(f"Count matrix: {counts.shape[0]} genes × {counts.shape[1]} samples")
    print(f"Conditions: {metadata['condition'].value_counts().to_dict()}\n")

    # Simulate DESeq2-like results (without R, for demo purposes)
    pvals = np.random.uniform(0, 1, n_genes)
    pvals[:100] = np.random.uniform(0, 0.001, 100)  # truly DE genes
    lfc = np.random.normal(0, 0.5, n_genes)
    lfc[:100] = np.random.normal(2.0, 0.8, 100)

    results = RosettaDataFrame({
        "baseMean": counts.mean(axis=1).values,
        "log2FoldChange": lfc,
        "lfcSE": np.abs(np.random.normal(0.3, 0.1, n_genes)),
        "stat": lfc / 0.3,
        "pvalue": pvals,
        "padj": np.minimum(pvals * n_genes / np.arange(1, n_genes + 1), 1.0),  # BH correction
    }, index=gene_names)

    print("─" * 40)
    results.report()
    print("─" * 40)

    # Show top genes
    sig = results[results["padj"] < 0.05].sort_values("log2FoldChange", ascending=False)
    print(f"\nTop 5 upregulated genes:")
    print(sig[["log2FoldChange", "padj"]].head().to_string())

    print("\n✓ To run with real R packages:")
    print("  results = rb.deseq2(counts, metadata, design='~ condition')")
    print("  results.report()")


if __name__ == "__main__":
    main()
