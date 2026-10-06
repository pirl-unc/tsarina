# Mortality-prioritized CTA vaccine design

## Public interface

Add `tsarina vaccine` and a Python API. Default to k=10, strict OncoRef core
CTAs, the existing global54_abc HLA panel, 8-mer subtraction, 8–11-mer class-I
ligands, MHCflurry, maximum 10-residue terminal padding, and RNA output.
`--cta-definition both` produces separate strict/loose designs, never mixes
their non-CTA backgrounds. Loose means OncoRef extended reproductive tissue
scope, including prostate, not every nominated CTA. Allow explicit HLA lists,
cohort mapping JSON, padding range, junction affinity threshold (1000 nM),
linker candidates (direct then AAY), DNA/RNA, custom or named Vaxrank UTRs,
polyA length, amino-acid and total-nucleotide limits, and output directory.

## Selection and identity

Use OncoRef's public genome-scope proteoform within-sample prevalence API at threshold .95;
identical-protein genes are summed before within-sample expression ranking.
Do not confuse this with the 95th percentile of TPM across patients.
Rank by sum over distinct burden categories of mortality share times p95
prevalence. This is an additive prioritization score, not unique-patient union
coverage or an estimate of preventable deaths. Retain incidence shares/counts,
mortality shares/counts, prevalence and sample denominators separately.
Default cohort mapping uses broad histology cohorts from the OncoRef report;
where multiple cohorts represent a category use sample-count-weighted
prevalence, counting that category's mortality once. Record this empirical
mixture assumption. Do not substitute narrow subtypes or residual categories.
Missing measurements remain missing; label scores as observed partial scores,
report unavailable genes/cohorts, and never fabricate zero incidence/counts.
Preserve full upstream burden/proteoform/CTA source tables and versions.

Resolve each candidate to the longest translated Ensembl protein, grouping
byte-identical full sequences. Verify its OncoRef proteoform expression key;
unregistered identical sequences require an upstream correction rather than
adding already-computed prevalences. Keep every contributing gene/transcript
and sequence SHA256. Rank top k before filtering; explicitly report dropped
targets instead of silently changing the requested cancer-priority selection.

## Sequence and evidence funnel

Walk all translated coding proteins, including noncanonical isoforms, for
non-CTA 8-mers. Mark the UNION of all residues covered by a shared 8-mer as
forbidden. The complement is a set of contiguous CTA-specific native intervals;
also break at ambiguous amino acids. Resolve exact-symbol HSCHR alternate
haplotype annotations of curated primary CTA genes before background screening,
recording every alias; primary non-CTA loci are never admitted by this rule
(Tsarina #187). Expression IDs and prevalence are unchanged. A proteoform shared
with an independent non-CTA gene
therefore has no surviving interval. Sharing with any other CTA is allowed.

Load current Hitlist human class-I MS observations with pushdown filtering to
peptides from surviving intervals. Reuse Tsarina's existing monoallelic,
sample-genotype/deconvolved, and unrestricted-MS inference and tier cutoffs;
do not promote peptide-only MS to measured allele restriction. Keep all
qualified supported peptides and repeated occurrences, not only top-per-cell.
Retain each specific interval with at least one support hit; trim only its
terminal padding, never cut an observed target ligand or stitch across removed
native sequence. Track raw → specific → MS-supported → padded → assembled
aa, fractions, interval counts, panel alleles, evidence tiers and drop reasons.

## Construct optimization

Search segment order and independent terminal-padding choices with a bounded,
deterministic beam search. Direct joins are preferred; consider AAY only when
needed. Optimize lexicographically: junction peptide/allele predictions below
1000 nM, then strongest binding, then predicted intersegment cleavage, then
linker cost/retained context. Use a real proteasomal predictor (Pepsickle via
mhctools), explicitly distinguish predictive evidence from proven cleavage.
Never replace unavailable model predictions with fabricated scores. Enumerate
all non-native 8–11-mers across EVERY final assembly boundary, including the
start methionine and both sides/inside linkers; include windows crossing
multiple short pieces. Re-score the actual final translated construct and
report unresolved binders, all junction scores, and all boundary cleavage
probabilities. Allow failure-on-unresolved-junctions as a strict option.

Length limits include added start M, linkers and stop codon; nt limit also
includes UTRs/polyA. Drop lower-priority whole supported pieces if needed,
preserving intact ligands and logging exclusions. Produce a single construct,
or explicit no-feasible-construct failure, never silent truncation.
Use Vaxrank's sourced HBB/HBB_FI named UTRs and the same DnaChisel codon
optimization method through an optional vaccine extra; custom UTR strings are allowed.
Avoid Vaxrank's full variant-pipeline dependency because it pins OncoRef 1.8.206,
which lacks the loose CTA API (Vaxrank #579, PirlyGenes #632). The shared
environment is restored to a consistent pinned stack; use the released 1.8.207
wheel in a temporary reference location for initial scientific checks. Final
validation uses a separate vaccine environment with Hitlist 1.64.7 and
OncoRef 1.8.207; reject older imported Hitlist code in live vaccine design.
DNA versus RNA changes alphabet,
not translation; validate round-trip translation. No automatic signal peptide,
MITD or additional translated elements without their requested biological role.

## Audit artifacts and validation

Write protein/CDS/full FASTAs, JSON manifest with configuration, provenance,
hashes, complete ranking and cohort contributions, interval and ligand tables,
per-protein funnel, construct layers, junction and cleavage tables, readable
Markdown report, and SVG figures (mortality/prevalence heatmap, sequence funnel,
construct map). The vaccine extra supplies the plotting dependency.
The focused OncoRef CTA registry omits RBMY1F/J despite identical sequences;
the genome registry includes them (OncoRef #567). Use the genome registry and
genome-scope within-sample ranking consistently, including every identical
protein group before the p95 calculation. This changes the reference rank axis
from earlier CTA-scope plots and must be recorded in every report.

Ship comprehensive method documentation and a reproducible offline fixture
example that exercises the full pipeline and is clearly labeled synthetic.
Exercise live OncoRef, Ensembl, Hitlist and real prediction models locally;
integration tests opt into model loading. Test exact grouping, missing data,
8-mer overlaps, ligand preservation, multi-boundary junctions, padding/linkers,
constraints, alphabets, reports, and CLI. Run format, lint and full tests;
bump to 1.33.0, PR, CI, merge, deploy from clean main and verify PyPI artifacts.
