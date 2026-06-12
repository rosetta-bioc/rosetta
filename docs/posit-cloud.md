# Posit Cloud Setup

## Recommended: Use Docker instead

Posit Cloud has known build issues with rpy2 (missing OpenMP library). For a guaranteed working environment:

```bash
docker pull ghcr.io/rosetta-bioc/rosetta:latest
docker run -it ghcr.io/rosetta-bioc/rosetta python3 -m rosetta
```

## If you must use Posit Cloud

In the **Terminal** tab (not R Console):

```bash
sudo apt-get install -y libomp-dev
pip install rosetta-bioc
python3 -m rosetta
```

Then install R packages:
```bash
Rscript -e "install.packages('BiocManager'); BiocManager::install(c('DESeq2','edgeR','limma'), ask=FALSE)"
```

## Notes

- The `libomp-dev` fix is needed because Posit Cloud's base image doesn't include OpenMP
- If this still fails, use the Docker container — it has everything pre-configured
- Free-tier Posit Cloud has limited RAM; large datasets may need a paid plan
