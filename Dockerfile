# Base image: Use rocker/r-ver for a stable R environment
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# 1. Install system-level dependencies required for R and Python compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gfortran libxml2-dev libssl-dev libcurl4-openssl-dev \
    libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libharfbuzz-dev libfribidi-dev libglpk-dev libgmp3-dev zlib1g-dev \
    libbz2-dev liblzma-dev libicu-dev r-base-dev \
    python3-pip python3-dev python3-setuptools python3-wheel \
    && rm -rf /var/lib/apt/lists/*

# 2. Set environment variables to ensure Python can locate R libraries
ENV R_HOME=/usr/local/lib/R
ENV LD_LIBRARY_PATH=$R_HOME/lib:$LD_LIBRARY_PATH
ENV PYTHONPATH=/app:$PYTHONPATH

# 3. Install R Bioconductor packages
# We install BiocManager first, then install the critical but smaller packages together,
# and finally install the resource-heavy packages separately.

RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

# Install essential dependencies and moderate-sized packages
RUN R -e "BiocManager::install(c('DESeq2', 'edgeR', 'limma', 'phyloseq', 'clusterProfiler', 'apeglm', 'ashr'), ask=FALSE, update=FALSE)"

# Install large/heavy packages separately to ensure build stability on memory-constrained systems
RUN R -e "BiocManager::install('org.Hs.eg.db', ask=FALSE, update=FALSE)"
RUN R -e "BiocManager::install('Seurat', ask=FALSE, update=FALSE)"
RUN R -e 'BiocManager::install("ReactomePA")'

# 4. Install Python core tools and rpy2 (pinned version to ensure compatibility)
RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir rpy2==3.5.15 && \
    pip install --no-cache-dir numpy pandas scipy

# 5. Install project dependencies
WORKDIR /app
COPY . .
RUN pip install --no-cache-dir '.[dev,plots,sklearn]'

# 6. Default command to run tests
CMD ["pytest"]