# Rosetta GSoC Meeting Notes — May 21, 2026

**Date:** Wednesday, May 21, 2026, 4:00 PM ET  
**Attendees:** Matias Salibian-Barrera (UBC), Catherine Chi Chung, John Muirhead-Gould (Nodes Bio)  
**Transcript:** `MasterPlan/02_AreasOfFocus/Nodes_Bio/Nodes Bio, Inc./00_Strategic_Notes/2026-05-21_Rosetta_GSoC_Meeting_Transcript.json`

---

## Key Decisions

1. **Wrap, don't reimplement** — Matias strongly endorsed the wrapping approach. The Bioconductor packages (DESeq2, edgeR, limma) have 20+ years of community testing and thousands of citations. A Python reimplementation would take a decade to earn equivalent trust. Rosetta's value is giving Python users access to the *original* trusted R code.

2. **Target audience** — Not hardcore R users (they'll stay in R). The target is Python-first bioinformaticians who are *reluctantly* using R because Bioconductor has no Python equivalent.

3. **Prioritize by usage** — Matias will identify which Bioconductor packages have the highest download counts / most demand from Python users. Focus rosetta's effort on those for maximum impact.

4. **Catherine confirmed** — Python reimplementations (e.g., PyDESeq2) have known statistical discrepancies vs. the original R. This validates the wrapping approach.

## Status

- Matias has AWS authentication working; hasn't run Kiro yet
- Catherine received new access link during the meeting
- John demoed Kiro CLI + quantum computing (QPU triage of RStudio bugs on IQM Garnet)

## Next Actions

| Who | Action | Due |
|-----|--------|-----|
| **Matias** | Ask colleague (+ others) which Bioconductor packages are most in-demand for Python users | Before next meeting |
| **Catherine** | Complete dev environment setup, run all 33 pytest tests | May 26 (GSoC start) |
| **John** | Send calendar invite for next meeting | This week |
| **All** | Email asynchronously for any blockers before next sync | Ongoing |

## Next Meeting

**Wednesday, May 27, 2026 — 1:00 PM Pacific / 4:00 PM Eastern**

## Notes

- Matias's colleague who would be best suited to answer the "which packages" question is on parental leave, but Matias will expand the ask to others
- The GSoC coding period officially starts May 26 — Catherine's first week is DESeq2 deepening (lfcShrink with apeglm/ashr/normal, contrast specification)
- The `shrink` parameter is already scaffolded in `rosetta/wrappers/deseq2.py`
