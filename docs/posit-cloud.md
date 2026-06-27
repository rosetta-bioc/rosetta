# Posit Cloud Setup

## ⚠️ Posit Cloud is not supported

Posit Cloud's environment is missing system libraries (`libomp-dev`) required to build rpy2, and `sudo` is not available to install them.

## Use Docker instead

```bash
docker pull ghcr.io/rosetta-bioc/rosetta:latest
docker run -it ghcr.io/rosetta-bioc/rosetta python3 -m rosetta
```

This gives you Python + R + Bioconductor + rosetta-bioc, fully configured, in one pull.

## Other options that work

- **Local install:** `pip install rosetta-bioc` (requires R + Bioconductor already installed)
- **GitHub Codespaces:** full sudo access, can install system deps
- **Any Linux/macOS with R:** just works
