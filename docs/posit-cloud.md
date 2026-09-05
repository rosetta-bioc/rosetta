# Posit Cloud Setup

## ✅ Posit Cloud is now supported (as of July 2025)

Posit Cloud updated its base build to include the system libraries previously
missing (`libomp-dev` and friends). As of the **Python 3.12 / R 4.6.1 / Bioc 3.23**
base image, `rosetta-bioc` installs and runs without any workarounds.

---

## Tested environment

| Component | Version |
|-----------|---------|
| Python | 3.12.11 |
| R | 4.6.1 |
| Bioconductor | 3.23 |
| rosetta-bioc | 0.3.2 |
| rpy2 | 3.6.7 |
| numpy | 2.5.2 |
| pandas | 3.0.5 |

---

## Installation

Open a **Terminal** in your Posit Cloud project and run:

```bash
pip install rosetta-bioc
```

That's it. `rpy2` and all Python dependencies install from wheels — no
compilation or system package installation required.

To install DESeq2 (and other Bioconductor packages), use R's console or an
R chunk in a notebook:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
```

DESeq2 1.52.0 installs cleanly on R 4.6.1 / Bioc 3.23.

---

## Running the test suite

The tests are not bundled in the wheel, so clone the source repo first:

```bash
git clone https://github.com/rosetta-bioc/rosetta /cloud/project/rosetta-src
```

Install test dependencies:

```bash
pip install pytest hypothesis matplotlib scikit-learn
```

### Pure-Python tests (all pass, ~5 seconds)

```bash
cd /cloud/project/rosetta-src
python -m pytest tests/test_bridge.py tests/test_codegen.py tests/test_deps.py \
    tests/test_detect.py tests/test_errors.py tests/test_init.py -v --tb=short
```

Result: **45/45 passed** on the tested environment above.

### R-integration tests (memory-constrained)

The full test suite — particularly `tests/test_normalize.py`,
`tests/stats/`, and `tests/property/` — loads DESeq2 into memory alongside
Python/rpy2. These tests are killed by the **1 GB memory cap** on Posit Cloud's
free/standard tier. This is a resource limit, not a compatibility issue:
the same tests pass on machines with ≥4 GB available.

To run the R-integration tests on Posit Cloud you need a project configured
with a larger memory allocation (≥4 GB recommended).

---

## Memory summary

| Test group | Memory required | Posit Cloud free tier (1 GB) |
|------------|----------------|-------------------------------|
| Pure-Python (bridge, detect, errors, …) | < 200 MB | ✅ passes |
| R-integration (normalize, stats, wrappers) | ~1–2 GB | ❌ OOM killed |
| Property-based (hypothesis + DESeq2) | > 2 GB | ❌ OOM killed |

---

## Demo notebook

**[→ Live demo on Posit Connect Cloud](https://01a072e1-a9dc-6e68-1843-4aaebba231fe.share.connect.posit.cloud/)**

[`posit-cloud-demo.Rmd`](posit-cloud-demo.Rmd) — end-to-end DESeq2 differential
expression on the `airway` dataset, entirely in Python. Loads counts from
Bioconductor, runs DESeq2, returns a pandas DataFrame, plots a volcano.
Open in RStudio on Posit Cloud and click **Knit**. Requires a ≥4 GB memory allocation.

---

## Other options

- **Docker:** `docker pull ghcr.io/rosetta-bioc/rosetta:latest` — fully
  configured, no memory cap issues for typical workloads
- **GitHub Codespaces:** full sudo access, configurable machine size
- **Local install:** `pip install rosetta-bioc` (requires R + Bioconductor
  already installed locally)
