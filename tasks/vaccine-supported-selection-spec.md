# Restricted, contributing CTA vaccine selection (2026-10-06)

## Objective

Design strict and loose ten-proteoform antigens with all MAGE family members
excluded except MAGEA4. Retain the original mortality-weighted individual
prevalence objective and whole-sequence identity collapse. Use actual MS
evidence, preserve non-CTA 8-mer subtraction and all final sequence constraints.
The configurable glob can narrow the exclusion to MAGEA when requested.

## Selection contract

- Add reproducible gene-symbol glob exclusions with explicit symbol exceptions;
  any excluded member vetoes an identical-sequence group. Exceptions affect
  candidate eligibility only, never CTA curation or the non-CTA background.
- Preserve existing ranked top-k semantics by default. Add an explicit supported
  selection mode, where k counts contributing proteoforms after subtraction,
  qualifying panel MS and minimum whole-segment length feasibility.
- Inspect eligible positive-score candidates in mortality rank order, in bounded
  batches. Record every inspected candidate's rejection/selection reason.
- In supported mode reserve one shortest ligand-preserving whole segment per
  selected proteoform before filling remaining sequence capacity by rank. A tight
  cap must not allow a high-ranked protein's extra pieces to crowd out all
  sequence from a lower-ranked selected protein.
- If fewer than k targets can contribute, write the audit and fail explicitly;
  never describe a smaller construct as a ten-protein design.

## Actual-MS contract (Hitlist #644 / Tsarina #189)

User clarification: keep the exact MS-observed peptide, infer every sample
allele with affinity below 1000 nM, and allow untyped MS with panel prediction.
Add explicit sample-affinity support mode alongside the previous presentation
tier mode. Do not require best-of-haplotype assignment or observed exact pMHC
identity. Preserve observation-to-allele assignment tables for both typed and
untyped evidence; study-wide allele pools are not individual sample genotypes.

Hitlist 1.64.7's nonbinding partition is not sufficient evidence of MS elution.
Vaccine selection requires positive MS modality: an explicit mass-spectrometry
assay method, curated MS-only supplement provenance with blank assay method, or
an explicit supplied MS-modality label when the method is absent. An explicit
non-MS method cannot be overridden. Unknown and non-MS rows are rejected with
their complete source metadata and reason. Preserve accepted/rejected input
tables and hashes even when no support remains. Require current imported
Hitlist, and do not restamp the stale full observations cache.

## Real validation and reporting

- Reuse the unchanged, hashed OncoRef/Ensembl input snapshots; make a new fresh
  scoped raw scan of every peptide needed for the expanded candidate set, with
  released Hitlist scanner/collector classification and current dedup/exclusions.
- Do not reuse the previous peptide-scoped MS file for newly inspected targets.
- Both definitions use global54_abc, RNA, HBB/HBB_FI UTRs, polyA 120, 1000-aa /
  3500-nt limits, padding step 2, beam 6, rounds 10. Record any inability to reach
  ten without relaxing scientific or length constraints.
- Publish selected/cancer/MS/HLA tables, complete funnels and dropout reasons,
  figures, sequence FASTAs and provenance. Replace old provisional examples with
  correctly labeled, newly validated results. FDA TECELRA approval motivates
  MAGEA4 eligibility; it does not approve this vaccine or other peptide-HLA pairs.
- Independently reconcile ten contributing full-sequence groups, actual MS
  modality, native/assembled ligand coordinates, non-CTA 8-mer absence, final
  junction windows, translation, length arithmetic and output hashes.

## Release

Regression-test exclusions, support backfill, tight-cap target reservation,
missing/explicitly non-MS evidence, failure audits and CLI propagation. Run
format, lint and full model-enabled tests. Version 1.33.1, feature-branch PR,
exact-head CI, merge, clean-main deploy and published artifact verification.
Link upstream #644 and downstream #189, and report remaining upstream work.
