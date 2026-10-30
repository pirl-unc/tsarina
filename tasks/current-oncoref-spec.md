# Current OncoRef selection migration

The user chose the current selection policy on 2026-10-29. The dependency-only
scope in the first PR draft is superseded by this specification.

## Contract

- Default human CTA membership comes exclusively from
  `oncoref.cta_selection.results(species="human", category="CTA")`.
  Preserve protein entity/member identities, policy/reference/source receipts,
  flags, and the distinction between Unsafe, Neither and null/unassigned.
- Expose the unmodified protein-level selection table and a gene-level join
  for existing consumers. Start that join from the current decision universe,
  then attach historical HPA annotations and Tsarina's MS overlay. Missing
  HPA data remain missing; legacy flags never define current admission.
- Keep explicitly named historical HPA access and strict/loose vaccine modes
  for reproducibility. Current membership helpers and default ranking must
  not reapply HPA confidence, vital-RNA cutoffs or a MAGE-family exclusion.
  Explicit user filters and downstream measured MS screening remain available.
- Default peptide/reference partitioning uses current selected member genes.
  All other coding genes, including Neither/unassigned, remain background.
  Do not treat Neither as never expressed or transplant older Disposable/Silent
  exclusions. Document this consumer background policy separately from OncoRef's
  independently versioned exact-protein policy. Sequence reference provenance
  remains distinct from the selection's Ensembl-115 gene-space receipt.
- Add `current` as the human vaccine definition/default; retain explicit
  `strict`, `loose`, and `both` historical comparisons. Dog frozen inputs retain
  their existing reviewed policy and cannot inherit human selection.
- Update tests to verify current-member parity, newly admitted genes absent
  from HPA, CTAG2 admission, null status preservation, background partitioning,
  and default CLI/API behavior. Preserve legacy evidence checks through the
  explicit historical accessor. Record the complete old/current membership delta.
- Fix the independently reproduced current Hitlist vocabulary compatibility
  failure without converting ambiguous typing into exact observed restriction.

## Work

- [x] Audit upstream selection schema, reference metadata and legacy boundaries.
- [x] Implement current membership/evidence and migrate default consumers.
- [x] Update CLI, vaccine provenance, compatibility tests and documentation.
- [x] Run formatting, lint, all tests, real-model checks and dependency checks.
- [ ] Rewrite PR #219 around the final scope, pass CI, merge, deploy and verify.

The shared environment remains constrained by Pirlygenes #658; use the isolated
release environment and preserve sibling work. Upstream policy thresholds are
consumed as published; this migration does not reimplement biological rules.

## Membership audit

The complete 225-added / 291-removed gene-ID delta is recorded in
[oncoref-1.23.1-membership-migration.json](oncoref-1.23.1-membership-migration.json).
Current 522 protein decisions flatten to 558 member genes. Historical HPA
annotations are absent for 92 selected members. Four selected identities have
no symbol and retain their stable gene IDs as display names.

Hitlist provenance compatibility is tracked in
[Tsarina #220](https://github.com/pirl-unc/tsarina/issues/220). Current-mode
vaccine backgrounds no longer inherit gene membership from same-symbol
alternate haplotypes; historical modes retain their explicit behavior.

Evidence-plot integration is tracked in
[Tsarina #221](https://github.com/pirl-unc/tsarina/issues/221). Default figures
select gene IDs, keep missing-HPA confidence strata, and count member loci.
The updated confidence figure was rendered and visually checked.

## Pre-merge verification

`./format.sh`, `./lint.sh`, strict MkDocs and isolated `pip check` pass.
`./test.sh --run-mhcflurry`: **761 passed, 2 optional NetMHCpan skips,
16 warnings, 86% coverage**. Real MHCflurry, Pepsickle, codon optimization
and local Ensembl checks pass. The final targeted current-policy, report and
historical-authority checks also pass. PR #219 carries merge/release evidence.
