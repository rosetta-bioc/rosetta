# JOSS Paper Revision Spec

**Source:** Reviewer feedback from Matias Salibian-Barrera (2026-07-02)
**Reference exemplar:** https://joss.theoj.org/papers/10.21105/joss.00862
**Status:** In progress — ORCIDs done (both authors); Task 1 examples DONE
(worked three-tier example added to paper.md, verified against 0.3.0 API);
Zenodo release/tag pending (see Task 2).

---

## ✅ Done
- [x] John Muirhead-Gould ORCID: `0009-0002-0470-5131`
- [x] Matias Salibian-Barrera ORCID: `0000-0003-1873-4611` (name correctly hyphenated)

## Open (Catherine)
- [ ] Catherine Chi Chung ORCID — none provided; optional, leave blank unless she supplies one.

---

## Task 1 — Add concrete tiered usage examples

**Problem (reviewer):** The paper "reads as a fairly high-level description of the
package, a bit on the abstract side." Currently only ONE code example exists
(the DESeq2 one-liner in Statement of Need). The reviewer wants examples that
show "how the different tiers differ, how they are used."

**Fix:** Add a new subsection `## Usage Examples` after `# Software Design`
(before `# Research Impact Statement`). Show the SAME analysis (DESeq2
differential expression) expressed at all three tiers so the reader sees the
ergonomics/control tradeoff directly.

### Example content to write (verify each runs against current API before committing):

```python
# Tier 1 — Quick API: one call, sensible defaults
import rosetta as rb
results = rb.quick_deseq2(counts_df, metadata_df, design="~ condition")
results.report()
```

```python
# Tier 2 — Class-based: stateful, chainable for multi-step pipelines
# NOTE: the .normalize()/.find_variable_features()/.scale()/.run_pca() chain
# below was ILLUSTRATIVE ONLY and does NOT match the real Seurat API. The
# actual chainable methods are run_standard_pipeline() / run_sctransform() /
# find_markers() / get_results(). The paper uses the verified version:
#     Seurat(sc_counts).run_standard_pipeline().find_markers(ident_1="0", ident_2="1")
from rosetta import Seurat
sc = (Seurat(counts_df)
        .normalize()
        .find_variable_features(n=2000)
        .scale()
        .run_pca())
sc.report()
```

```python
# Tier 3 — Functional: explicit control over each intermediate object
from rosetta import run_deseq2, get_results, lfc_shrink
dds = run_deseq2(counts_df, metadata_df, design="~ condition")
res = get_results(dds, contrast=("condition", "treated", "control"))
shrunk = lfc_shrink(dds, coef="condition_treated_vs_control", type="apeglm")
```

**Requirements:**
- Each snippet MUST be runnable against the current public API. Check actual
  function/method names in `rosetta/` — the Tier 2/3 names above are from the
  paper's own description and need verification against real signatures.
- Add 1-2 sentences of prose between snippets explaining WHEN a user picks each tier.
- Keep it tight — JOSS papers are short (~1000 words). Don't over-explain.
- Consider showing sample `.report()` output (a few lines) to make it concrete.

## Task 2 — Zenodo archiving setup (required by JOSS on acceptance)

Per reviewer: "upon acceptance, you will need to archive a copy of the current
version of the package with zenodo.org."

- [ ] Log into https://zenodo.org with GitHub, enable the Zenodo–GitHub webhook
      for the rosetta repo (Settings → GitHub → toggle repo ON).
- [ ] Cut a tagged GitHub release matching the current package version. The
      package is on the 0.3.x line (`rosetta.__version__ == "0.3.0"`), and the
      top-level Tier 3 export change means the release tag must be bumped
      accordingly (e.g. `v0.3.1`). Do NOT reuse `v0.1.0` — that predates the
      current API. Zenodo auto-archives the tag and mints a DOI.
- [ ] Add the resulting Zenodo DOI badge to `README.md`.
- [ ] Do the release AT acceptance time (JOSS wants the archived version to match
      the accepted paper), but wire up the integration now so it's one click.

## Task 3 — Pre-submission checklist
- [ ] `paper.bib` — confirm all `@cite` keys resolve (huber2015orchestrating,
      gautier2010rpy2, muzellec2023pydeseq2).
- [ ] `architecture.svg` present and renders.
- [ ] Build paper locally with the JOSS Docker action / `openjournals/inara` to
      confirm it compiles to PDF cleanly before submitting.
- [ ] Notify Matias when revised paper.md is pushed (he offered a closer read).

