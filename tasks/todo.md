# Curation And Tissue Logic Audit

## Spec

- [x] Read project instructions in `AGENTS.md`.
- [x] Check for prior lessons in `tasks/lessons.md`.
  Result: file does not exist yet.
- [x] Inspect `tsarina` CTA curation and tissue restriction code paths.
  Focus: reproductive tissue sets, HPA thresholding, CTA evidence/tier annotations, and any documentation or test mismatches.
- [x] Inspect `hitlist` source classification code as used by `tsarina`.
  Focus: tissue category loading, healthy/reproductive/thymus labeling, cancer-specific classification, and potential annotation drift.
- [x] Cross-check implementation against tests and exposed docs.
  Focus: places where code behavior differs from README/docs wording or expected curation semantics.
- [x] Write a report with concrete findings.
  Deliverable: a markdown report with severity, evidence, and recommended fixes.

## Review

- Report written to `tasks/curation_tissue_logic_report.md`.
- Highest-risk finding: `tsarina/iedb.py` bypasses key `hitlist` curation inputs, so PMID overrides are not applied in `tsarina` scans.
- Highest-risk CTA finding: the default CTA set includes 16 genes whose shipped synthesized restriction is `SOMATIC`.
