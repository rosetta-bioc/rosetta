# Glossary

Everything in this project explained from first principles. No prior biology or R knowledge assumed.

## The Problem rosetta Solves

**The two-language problem** — Bioinformaticians write their pipelines in Python but the gold-standard analysis tools (DESeq2, edgeR, limma) only exist in R. Today, switching between them means breaking your workflow, manually converting data formats, and debugging cryptic cross-language errors.

**rosetta** — A Python interface that lets you call these R tools with one line of code. You pass in a pandas DataFrame, you get a pandas DataFrame back. All the statistical rigor of the original R packages, none of the context-switching.

## Biology (What the Tools Analyze)

**Gene** — A section of DNA that contains instructions for making a protein. Humans have ~20,000 genes.

**Gene expression** — How actively a cell is reading a particular gene. Genes can be highly expressed, barely expressed, or silent. Think of it as a volume dial per gene.

**RNA-seq** — A lab technique that measures gene expression across an entire genome at once. Produces a count of how many times each gene was read in a sample.

**Count matrix** — The raw data from RNA-seq: a table where rows are genes, columns are samples, and each cell is "how many times we observed this gene in this sample."

**Differential expression analysis** — The core question: "Which genes changed between my two conditions?" (e.g., healthy vs. diseased, treated vs. untreated). This is what DESeq2, edgeR, and limma compute.

**Log-fold change (LFC)** — How much a gene's expression changed, on a log2 scale. LFC of 1 means it doubled. LFC of -1 means it halved.

**Shrinkage (lfcShrink)** — Genes with very few counts produce noisy, unreliable fold-change estimates. Shrinkage stabilizes these by pulling extreme estimates toward zero — keeping only changes supported by strong evidence.

**Contrasts** — The specific comparisons you want to make. If your experiment has conditions A, B, and C, a contrast says "compare B vs. A" or "compare C vs. the average of A and B."

**p-value** — The probability of seeing your result by random chance alone. Small p-value = likely a real biological effect.

**Adjusted p-value (padj / FDR)** — Corrected for testing thousands of genes simultaneously. Without adjustment, testing 20,000 genes guarantees false positives.

**Microbiome** — The community of microorganisms living in an environment (e.g., the human gut). Analyzed with phyloseq.

**Single-cell RNA-seq** — Measuring gene expression in individual cells rather than bulk tissue. Reveals cell-type diversity. Analyzed with Seurat.

## The R Packages (What rosetta Wraps)

**DESeq2** — The most widely-cited tool for differential expression. Rigorous statistical model, trusted by reviewers and journals.

**edgeR** — Alternative to DESeq2 with a different statistical approach. Often used as a cross-validation ("do both tools agree?").

**limma-voom** — Originally built for microarray data, extended to RNA-seq. Fast, flexible, battle-tested over 20 years.

**Seurat** — The standard for single-cell analysis. Clustering, visualization, cell-type annotation.

**clusterProfiler** — Gene set enrichment: "I found 200 differentially expressed genes — what biological pathways are they involved in?"

**phyloseq** — Microbiome community analysis: diversity, composition, comparison between conditions.

**Bioconductor** — R's curated repository for bioinformatics packages. Peer-reviewed, version-controlled, the standard distribution channel for these tools.

**apeglm / ashr** — Shrinkage algorithms that DESeq2 calls internally. Different mathematical approaches to the same stabilization goal.

## The Technology Stack

**rpy2** — The low-level Python↔R bridge. Handles raw interprocess communication but requires you to manage type conversion, R sessions, and error handling manually. Powerful but verbose.

**rosetta** — The ergonomic layer on top of rpy2. Handles conversion, validation, session management, and error translation automatically. rpy2 is to rosetta as urllib is to requests.

**pandas** — Python's standard library for tabular data. The universal format for data science in Python.

**DataFrame** — A table in memory (rows and columns). The common currency between pandas and R.

**Type conversion** — Translating between Python and R data formats. A pandas DataFrame and an R data.frame represent the same concept but are different objects in memory.

**Namespace** — Calling `DESeq2::DESeq()` instead of just `DESeq()`. Specifies which package a function belongs to, preventing conflicts. A safety practice rosetta enforces.

**Design formula** — R's notation for describing experimental structure. `~ condition` means "model expression as a function of the condition variable." Passed as a string in rosetta.

## Project Context

**GSoC (Google Summer of Code)** — A global program where Google funds contributors to build open-source software, paired with experienced mentors.

**The distribution problem** — Building great software isn't enough if nobody uses it. Adoption requires clear documentation, accessible language, and meeting users where they are.

**Subprocess fallback** — A planned resilience feature: if rpy2's installation breaks (a known issue), rosetta falls back to running R as a separate process. The user's code doesn't change.
