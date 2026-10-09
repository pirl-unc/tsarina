# Canine evidence and vaccine policy

Implement the offline migration in #208 without making Canvax or unpublished
upstream patches runtime dependencies. Human inputs, defaults and outputs stay
compatible. A dog construct is exploratory, never a validated vaccine.

## Input contract

A versioned JSON design bundle contains taxon 9615; a reference identity with
assembly, annotation release and asset SHA256s; every complete translated source
occurrence; per-occurrence reproductive-restriction admission and reasons; named
cohort sample/donor metadata; exact-sequence RNA abundance bounds; positive-MS
facets with independent peptide, host and MHC identities; per-allele model
capabilities; and explicitly observed or simulated tumor/genotype pairs.
Reference keys are the SHA256 of canonical JSON reference metadata. Protein
groups are recomputed from full sequences; all occurrences remain searchable.
No symbol-based inheritance, human CTA mapping or cross-species orthology rescue.
Rejected and unknown occurrences remain in the audit and specificity background.

An admitted occurrence requires a supplied policy/normal assessment with allowed
tissues and provenance. Verified healthy primary RNA or MS outside that scope
supplies counterevidence. This is an import of upstream reviewed evidence, not a
replacement for Canvax's expression allocation or CTA discovery. Unknown normal
assessment cannot admit a candidate. Reference/unit/identity mismatches fail.

## Ranking and coverage

The caller chooses one named untreated primary-tumor cohort and an absolute TPM
threshold from the bundle. Independent reconciled dogs are the denominator.
Repeated biopsies use all-positive lower/any-possible upper bounds; missing RNA
stays unknown. Unresolved donor/pooling/QC/preparation rows remain visible and
are not pooled into that denominator. No mortality or within-sample p95 priors.

Use target-specific expression plus qualified DLA evidence on supplied full
genotypes, retaining more than two same-locus alleles. Joint coverage lives in
the shared coverage module (#201), with observed versus simulated pairing,
unsupported and missing mass, and lower/upper bounds. No HWE or breed/global
population claim without an explicitly defined sample frame.

## Native support and assembly

Retain exact positive MS-observed peptide strings only. Human-host DLA
monoallelic evidence is distinct from endogenous dog presentation. Typed
multiallelic support requires affinity below the configured cutoff to an allele
actually in that sample; untyped support requires explicit opt-in. Never assign
an observed long peptide to a different nested peptide. Model capability is
length-specific; sequence-extrapolated selection requires explicit opt-in.

Subtract non-admitted occurrence 8-mers and verified healthy-primary MS 8-mers,
including longer normal peptides. Reuse native interval/segment objects and
bounded assembly search. Audit the entire final product, including initiating
methionine and synthetic linkers, against the entire background. Missing DLA
prediction entries fail rather than pass. Human Pepsickle and percentile gates
are not implicit dog defaults. Unassessed cleavage contributes no preference.
The final construct retains unassessed status even when evaluated junctions pass.

Positive-MS admission consumes an explicit reviewed source boolean and reason,
retaining raw method/response/qualitative fields. Host taxonomy and actual MS
sample identities are independent facets. The genotype sampling frame must
name the selected RNA cohort. Portable protein/provenance drift is rejected.

Frozen affinity tables enable reproducible offline native and junction scoring;
an injected runtime predictor is supported with explicit callback provenance.
Live model/key/capability integration remains tracked upstream. Generic exact
reverse translation is the dog default because the installed codon-table library
has no dog table; explicit alternative codon species is recorded, never called
dog optimization. UTRs require explicit sequences in canine mode. All existing
AA/nt limits, stop, polyA, padding and linker constraints apply.

## Deliverables and validation

Add species/input-bundle/cohort/exploratory CLI options and preserve portable
VaccineInputs by adding an optional canine policy field at the end. Emit raw
bundle hash, all source tables, ranking, independent-dog prevalence, exact-MS
assignments/rejections, funnel, native maps, capability and cumulative paired
coverage, final sequence/whole-product audit and plain-language HTML/Markdown.
RNA tables are emitted even if downstream support is absent.

Offline fixture tests: repeated dogs, missing expression, identical adult-heart
source veto, forbidden normal MS long-peptide overlap, duplicate DLA loci,
unsupported/missing mass, exact versus inferred restriction, no nested evidence,
mixed references/taxa/units, exploratory opt-in, incomplete prediction archive,
whole-product boundary veto, exact nt translation/cap and existing human tests.
Run format, lint, full suite with real human models, strict docs build, wheel
asset/import checks, PR CI, merge and clean-main PyPI deployment.

Dependencies: oncoref#571, hitlist#660/#661, mhctools#544, mhcflurry#490,
vaxrank#596; Tsarina #199/#191/#192/#201 remain independently tracked.
