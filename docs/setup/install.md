# Installation

## Python package

```bash
pip install rosetta-bioc
```

### Optional dependencies

```bash
# For plotting (volcano, MA, PCA plots)
pip install rosetta-bioc[plots]

# For development
pip install rosetta-bioc[dev]
```

---

## Requirements

| Component | Minimum version |
|-----------|----------------|
| Python | 3.9+ |
| R | 4.0+ |
| rpy2 | 3.5+ |
| pandas | 1.5+ |
| numpy | 1.23+ |

---

## R setup

Rosetta calls R through rpy2 — you need a working R installation with Bioconductor packages.

### Install R

=== "macOS"

    ```bash
    # Homebrew
    brew install r

    # Or download from CRAN
    # https://cran.r-project.org/bin/macosx/
    ```

=== "Ubuntu/Debian"

    ```bash
    sudo apt-get update
    sudo apt-get install r-base r-base-dev
    ```

=== "Fedora/RHEL"

    ```bash
    sudo dnf install R R-devel
    ```

### Install Bioconductor packages

From an R console:

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "limma"))
```

Or use rosetta's install script:

```bash
Rscript install.R
```

See [R Packages](r-packages.md) for which R packages are needed for each wrapper.

---

## Docker

The easiest way to get a fully configured environment:

```bash
docker pull ghcr.io/rosetta-bioc/rosetta:latest
docker run -it ghcr.io/rosetta-bioc/rosetta python3
```

The Docker image includes:

- Python 3 with rosetta-bioc installed
- R 4.3 with all Bioconductor packages (DESeq2, edgeR, limma, Seurat, phyloseq, clusterProfiler)
- Annotation databases (org.Hs.eg.db)
- rpy2 configured and tested

### Building from source

```bash
git clone https://github.com/rosetta-bioc/rosetta.git
cd rosetta
docker build -t rosetta .
docker run -it rosetta python3
```

### Mounting local data

```bash
docker run -it -v $(pwd)/data:/data rosetta python3
```

```python
import pandas as pd
counts = pd.read_csv("/data/counts.csv", index_col=0)
```

---

## Posit Cloud

!!! warning
    Posit Cloud is **not supported**. It lacks system libraries (`libomp-dev`) required to build rpy2, and `sudo` is not available.

Use Docker instead, or try GitHub Codespaces which provides full `sudo` access.

---

## Verifying installation

```python
import rosetta as rb
print(rb.__version__)

# Quick check that R bridge works
import rpy2.robjects as ro
print(ro.r("R.version.string")[0])
```

### Checking available R packages

```python
from rosetta._deps import is_installed

print(f"DESeq2: {is_installed('DESeq2')}")
print(f"edgeR: {is_installed('edgeR')}")
print(f"limma: {is_installed('limma')}")
print(f"Seurat: {is_installed('Seurat')}")
print(f"phyloseq: {is_installed('phyloseq')}")
print(f"clusterProfiler: {is_installed('clusterProfiler')}")
```

---

## Troubleshooting

### `rpy2` won't install

Make sure R is in your PATH and `R_HOME` is set:

```bash
which R          # Should print path to R
R RHOME          # Should print R home directory
echo $R_HOME     # Should be set (or will be auto-detected)
```

On macOS with Homebrew R:

```bash
export R_HOME=$(R RHOME)
pip install rpy2
```

### `RPackageMissing` error

Install the missing R package from an R console:

```r
BiocManager::install("PackageName")
```

### Namespace collision with `rosetta` package

If you have the unrelated `rosetta` PyPI package installed:

```bash
pip uninstall rosetta
pip install rosetta-bioc
```

Rosetta will warn you on import if both packages are detected.
