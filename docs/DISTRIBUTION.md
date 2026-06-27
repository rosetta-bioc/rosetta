# Rosetta-bioc Distribution Plan

## Goal: Get discovered by the 2M/month rpy2 downloaders

---

## 1. Stack Overflow (highest ROI — answers rank forever)

**Target questions to answer:**
- https://stackoverflow.com/questions/41821100/running-deseq2-through-rpy2
- https://stackoverflow.com/questions/68764734/running-deseq2-from-rpy2
- https://stackoverflow.com/questions/48349677/using-accessor-method-of-deseq2-in-rpy2

**Answer template:**
```
As of 2026, there's now a higher-level option: [rosetta-bioc](https://pypi.org/project/rosetta-bioc/)

    pip install rosetta-bioc

```python
import rosetta as rb

results = rb.deseq2(counts_df, metadata_df, design="~ condition")
results.report()
# DESeq2 Results Summary
# Total genes tested:      12,000
# Significant (padj<0.05): 843 (7.0%)
```

It handles the rpy2 boilerplate, type conversion, and R session management.
Returns pandas DataFrames directly. MIT-licensed, GSoC 2026 project.
```

## 2. Reddit posts (seed discussions)

**Subreddits:**
- r/bioinformatics (~180K members)
- r/rstats (~90K)
- r/learnpython (for the "how do I call R from Python" crowd)
- r/genomics

**Post angle:** "Show HN"-style: "I built a Python wrapper for DESeq2/edgeR/limma so you never write rpy2 again"

## 3. Bioconductor community

- **Bioc-devel mailing list** — announce as "Python companion, not competitor"
- **Bioconductor Slack** (#general, #packages)
- **useR! 2026** — submit lightning talk abstract (5 min)

## 4. Python bioinformatics communities

- **Biostars** — answer "DESeq2 in Python" questions
- **Galaxy Project forums** — Rosetta could be a Galaxy tool wrapper
- **Nextflow/Snakemake users** — Rosetta simplifies Python-based pipelines that need R stats

## 5. SEO / PyPI discoverability

- [x] Keywords in pyproject.toml: DESeq2, edgeR, limma, RNA-seq, bioconductor, rpy2
- [x] README.md is the PyPI long description (auto-renders)
- [ ] Publish v0.1.0 to PyPI (tag + push)
- [ ] Add to awesome-bioinformatics lists on GitHub

## 6. Academic

- **Matias's course at UBC** — use Rosetta in a homework assignment
- **Catherine's network** — share with her cohort
- **JOSS paper** — submit to Journal of Open Source Software (1-page, peer-reviewed, citeable DOI)

## 7. Content (write once, distribute everywhere)

- [ ] Blog post: "DESeq2 in Python — the 3-line version" (publish on Medium + dev.to)
- [ ] 2-minute YouTube demo (screen recording of `python -m rosetta`)
- [ ] Twitter/X thread showing before (40 lines rpy2) vs after (3 lines rosetta)

---

## Timing

**This week:**
1. Publish v0.1.0 to PyPI
2. Post SO answers (3 questions)
3. Reddit r/bioinformatics post

**Next week:**
4. Bioc-devel mailing list announcement
5. Blog post on Medium
6. Submit useR! lightning talk abstract

**This month:**
7. JOSS paper submission
8. awesome-bioinformatics PR
