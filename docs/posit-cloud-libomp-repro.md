# Repro: `pip install rosetta-bioc` fails on Posit Cloud

**Environment:** Posit Cloud free tier, RStudio project, Terminal tab  
**Date:** 2026-08-14  
**Package:** `rosetta-bioc` (rpy2 ≥ 3.5.15 dependency)

## Steps to Reproduce

1. Create a new RStudio project on Posit Cloud
2. Open the Terminal tab
3. Run:

```bash
pip install rosetta-bioc
```

## Root Cause

R 4.6.1 on Posit Cloud is compiled with `-fopenmp`, so `rpy2-rinterface`'s build
system picks up `-fopenmp` from `R CMD config --ldflags` and passes it to the
linker. The linker then looks for `-lomp`, which is not installed in the base image.

The fix: install `libomp-dev` (or equivalent) in the Posit Cloud base image.

## Error Output

```
Looking for R CONFIG with: /opt/R/4.6.1/lib/R/bin/R CMD config --ldflags
['-Wl,--export-dynamic -fopenmp -L/usr/local/lib -L/opt/R/4.6.1/lib/R/lib -lR ...']

clang -pthread -shared ... \
  -Wl,--export-dynamic -fopenmp \
  -Wl,-rpath,/opt/R/4.6.1/lib/R/lib

/usr/bin/ld: cannot find -lomp: No such file or directory
clang: error: linker command failed with exit code 1

error: Command '['clang', '-pthread', '-shared', '-Wl,--exclude-libs,ALL',
  'build/temp.linux-x86_64-cpython-312/build/temp.linux-x86_64-cpython-312/_rinterface_cffi_api.o',
  '-L/usr/local/lib', '-L/opt/R/4.6.1/lib/R/lib', '-L/opt/python/3.12.11/lib',
  '-lR', '-lpcre2-8', '-ldeflate', '-lzstd', '-llzma', '-lbz2', '-lz',
  '-ltirpc', '-lrt', '-ldl', '-lm', '-licuuc', '-licui18n', '-lblas',
  '-o', 'build/lib.linux-x86_64-cpython-312/_rinterface_cffi_api.abi3.so',
  '-Wl,--export-dynamic', '-fopenmp',
  '-Wl,-rpath,/opt/R/4.6.1/lib/R/lib']' returned non-zero exit status 1.

ERROR: Failed building wheel for rpy2-rinterface
ERROR: Failed to build installable wheels for some pyproject.toml based projects (rpy2-rinterface)
```

## Suggested Fix

```bash
# On the Posit Cloud base image (Ubuntu/Debian):
apt-get install -y libomp-dev
```

This installs `libomp.so` so the linker can resolve `-lomp` at build time.
This unblocks `rpy2`, `rosetta-bioc`, and any other Python/R bridge package
built against an OpenMP-enabled R installation.
