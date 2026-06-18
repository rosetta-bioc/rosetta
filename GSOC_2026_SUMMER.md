# GSoC 2026 Summer — Hardening Bioconductor Wrappers

**Contributor:** Catherine Chi Chung  
**Mentors:** John Muirhead-Gould (john@nodes.bio), Prof. Matias Salibian-Barrera (matias@stat.ubc.ca)  
**Organization:** R Project for Statistical Computing  
**Project size:** Large (350 hours, ~12 weeks)  
**Dates:** May 26 – August 25, 2026

---

## Project Summary

Catherine will deepen and harden rosetta's Bioconductor wrappers — adding publication-critical parameters, building a three-tier API, implementing a subprocess fallback, and expanding test coverage. See [GSOC_PROPOSAL.md](GSOC_PROPOSAL.md) for the full technical scope and [SPEC.md](SPEC.md) for the current architecture.

## Communication

| Channel | Purpose | Cadence |
|---------|---------|---------|
| GitHub Issues/PRs | All technical work, code review | Continuous |
| Email thread (Catherine + John + Matias) | Weekly status update | Every Monday |
| Video call | Sync on blockers, design decisions | Weekly (schedule TBD in bonding period) |

All design decisions and technical discussions should happen in GitHub issues so the community can follow along. Email/calls are for coordination, not architecture.

## 12-Week Timeline

### Community Bonding (May 8 – May 25)

- [ ] Set up local dev environment (R 4.0+, rpy2, Bioconductor)
- [ ] Run existing test suite (`pytest`) — all 33 tests should pass
- [ ] Read `SPEC.md` and `CONTRIBUTING.md`
- [ ] Review Catherine's existing test repo: [usacchung/Rosetta-Test-Implementation](https://github.com/usacchung/Rosetta-Test-Implementation)
- [ ] Schedule weekly call time with John and Matias
- [ ] Open a GitHub issue: "GSoC 2026 Tracking" to log weekly progress

### Phase 1: Deepen Wrappers (Weeks 1–4, May 26 – June 22)

**Goal:** Add publication-critical parameters to all six wrappers.

| Week | Focus | Deliverable |
|------|-------|-------------|
| 1 | DESeq2 | `lfcShrink()` (apeglm/ashr/normal), `contrast`, `resultsNames()`, LFC thresholds |
| 2 | edgeR + limma | Contrast matrices, TREAT testing, quasi-likelihood, `decideTests()` |
| 3 | clusterProfiler | GSEA support, KEGG/Reactome pathways, custom gene sets |
| 4 | Seurat + phyloseq | `FindMarkers()`, SCTransform, ordination, differential abundance |

**Acceptance criteria:**
- Each new parameter has at least 2 pytest tests (happy path + edge case)
- Docstrings with usage examples for every new kwarg
- PR per wrapper, reviewed and merged

### Phase 2: Three-Tier API + Subprocess Fallback (Weeks 5–8, June 23 – July 20)

**Goal:** Restructure the API and add a fallback backend.

| Week | Focus | Deliverable |
|------|-------|-------------|
| 5 | Tier 1 (quick defaults) | Refactor existing single-call API, ensure backward compatibility |
| 6 | Tier 2 (granular kwargs) | Expose all R parameters as Python kwargs with validation |
| 7 | Tier 3 (R escape hatch) | Expose underlying rpy2 objects for advanced users |
| 8 | Subprocess fallback | `subprocess` + `Rscript` + JSON backend, auto-detection, version pinning |

**Midterm evaluation (July 14–18):**
- Phases 1–2 code merged to main
- All existing tests still pass (no regressions)
- New test count target: 60+ (up from 33)

### Phase 3: Testing + Documentation (Weeks 9–12, July 21 – August 25)

| Week | Focus | Deliverable |
|------|-------|-------------|
| 9 | Validation suite | Compare rosetta output vs direct R output for published datasets (airway, pasilla) |
| 10 | Integration tests | Real R tests (not mocked) for all wrappers, CI configuration |
| 11 | Documentation | Vignette-style tutorial notebook, API reference, update README |
| 12 | Polish + blog post | Final cleanup, GSoC summary blog post for R community |

**Final evaluation (August 25–September 1):**
- 80%+ test coverage
- All wrappers validated against direct R output
- Tutorial notebook published
- Blog post draft submitted

## How to Submit Work

1. Create a feature branch: `gsoc/week-N-description`
2. Open a draft PR early — don't wait until the code is perfect
3. Request review from `@jmg421` (John) and tag Matias on stats-heavy decisions
4. Each PR should include: code, tests, and docstring updates
5. Keep PRs focused — one wrapper or one feature per PR

## Key Resources

- [SPEC.md](SPEC.md) — Architecture and design decisions
- [CONTRIBUTING.md](CONTRIBUTING.md) — Code style, debugging practices, PR guidelines
- [GSOC_PROPOSAL.md](GSOC_PROPOSAL.md) — Full project proposal with technical details
- [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [edgeR user guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)
- [limma user guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)

## Evaluation Criteria

Per GSoC guidelines, evaluations assess:
- **Code quality** — clean, tested, documented
- **Communication** — responsive, proactive about blockers
- **Deliverables** — milestones met on time
- **Community engagement** — work done in the open, PRs reviewable

## Questions?

Open a GitHub issue or email John and Matias. No question is too small during community bonding.
