# Statistics Examples

End-to-end examples showing rosetta wrappers on real published datasets.
Each script runs a full analysis pipeline — data in, results out — and can be
run independently.

## Structure

```
examples/statistics/
├── deseq2_airway.py          # DESeq2 on the airway dataset (Himes et al. 2014)
├── edger_pasilla.py          # edgeR on the pasilla dataset
├── limma_voom_example.py     # limma-voom on a two-group RNA-seq design
└── README.md
```

## Running an Example

```bash
# From the repo root
pip install -e ".[dev]"
python examples/statistics/deseq2_airway.py
```

Requires R 4.0+ with BiocManager and the relevant Bioconductor packages installed.
Running an example for the first time will prompt you to install any missing R packages.

## Adding a New Example

1. One file per wrapper / dataset combination
2. Use only public, published datasets (airway, pasilla, GSE datasets)
3. Print a summary of results at the end — genes found, top hits, runtime
4. Keep it self-contained — no imports from other example files
