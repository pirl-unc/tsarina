# Vaccine atlas and length-budget selection

Build a reusable self-contained website beside each scientific vaccine report;
combine strict/loose reports without rerunning models. Publish real results
through the existing MkDocs GitHub Pages deployment. Visitor copy describes
biology, inputs, evidence and limitations without development history. Support
arbitrary target count, allele panel and DNA/RNA settings.

Keep ranked/supported top-k modes. Add budget selection requiring an aa or
total-nt cap, screening all positive-scoring eligible targets and allocating
whole native ligand-bearing segments by marginal mortality-weighted expression
gain per amino acid, then incidence gain and new peptide/HLA support. Recompute
gain after each addition. Use max(p) as an overlap-conservative surrogate, not
actual patient coverage. Allow additional pieces of an existing protein to add
ligands. Preserve selected pieces during junction optimization; report the
heuristic, finite candidate scope, constraints and stop reason.

Explain actual OncoRef tissue scopes: strict testis/ovary/placenta; loose adds
cervix/endometrium/epididymis/fallopian tube/prostate/seminal vesicle/vagina.
Breast is not included. Explain thymus and protein/RNA curation conventions;
each definition has its own non-CTA background. Compare CT83 membership,
ranking and all screening steps independently of which construct selected it.

Calculate HLA carrier probabilities by summed allele frequencies within locus
under Hardy–Weinberg equilibrium; combine loci under linkage equilibrium.
Use published global-reference proxy values, never panel-normalize or silently
clip invalid sums; missing frequencies are explicit. Show sensitivity to
measured, typed-inferred and untyped-inferred restrictions. For cancer p95
marginals show union bounds [max(p), min(1,sum(p))], missing measurements and
represented burden. Mortality/incidence curves have separate denominators.
Show cumulative proteins and final-order segments versus actual aa length,
distinct MS-observed peptides and pMHC pairs colored by strongest restriction
evidence. Do not call these T-cell-validated epitopes or clinical coverage.

Map full source proteins with CTA-specific intervals, retained pieces and
positive-MS observations from cancer, healthy nonreproductive, reproductive
and unknown samples. Query full proteins, including removed regions. Preserve
tissue, disease, cell-line/ex-vivo status, HLA restriction and study references.
Do not label lung-cancer observations as healthy lung. Mark cancer-only within
the queried corpus, not intrinsically cancer-specific. Show normal-tissue
overlap with chosen pieces. After the user requested source-verified cardiac
exclusion, inspect the original Atlas donor/sample tables and publication.
Use Hitlist's donor-resolved nonmalignant heart/brain/lung adapter and blacklist
audit, with an explicit one-donor threshold for this conservative design.
Subtract 8-mers from qualifying primary nonmalignant HLA-I peptides. IEDB
healthy flags alone do not establish donor/sample verification. Keep the gate
optional and retain audit-only comparisons. Neither mode establishes safety.

Package HTML/CSS/JS/fonts and use one payload builder for fresh reports and a
standalone vaccine-report command accepting completed manifests. Verify saved
artifact hashes. Export coverage/tissue CSVs, SVG/PNG figures, sequences and
source tables. Test analytical coverage, missing/invalid data, marginal budget
selection, coordinates, deduplication, tissue classification, arbitrary configs,
single/both-definition output and packaging. Independently verify real results,
run format/lint/full tests and strict MkDocs, check browser controls/layout,
then PR/CI/merge/clean-main PyPI deployment and live Pages verification.
