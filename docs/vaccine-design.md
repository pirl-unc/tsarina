# Mortality-prioritized CTA vaccine design

`tsarina vaccine` selects cancer-testis antigen **proteoforms**, retains native
CTA-specific regions with panel-supported MS ligands, and assembles one
DNA/RNA antigen with complete tables, figures, and a final junction audit.

The [real-data validation](vaccine-validation.md) includes strict/loose selected
proteoforms, per-cancer incidence/mortality/p95 tables, retention figures,
sequence files and unresolved junction predictions.

```sh
python -m venv .venv-vaccine
. .venv-vaccine/bin/activate
pip install 'tsarina[vaccine]'
pyensembl install --release 112 --species human
mhcflurry downloads fetch models_class1_presentation

# Register MS data via tsarina data; build with tsarina build observations.
tsarina vaccine --top-k 10 --cta-definition both --auto-fetch \
  --panel global54_abc --vaccine-type mrna \
  --include-utrs --utr-5p HBB --utr-3p HBB_FI --poly-a-length 120 \
  --max-length-aa 686 --max-length-nt 2500 -o vaccine-out

# Custom alleles, coding DNA, fine-grained edge search, direct joins only.
tsarina vaccine -k 5 --hla 'HLA-A*02:01,HLA-A*24:02,HLA-B*07:02' \
  --vaccine-type dna --max-padding 10 --padding-step 1 --linkers '' \
  --require-clean-junctions -o dna-out

# Ten contributing proteins; MAGEA4 is the sole eligible MAGE-family member.
# Exact MS peptides qualify by affinity to any sample allele; untyped MS can
# use panel predictions. Quote the glob so your shell does not expand it.
tsarina vaccine -k 10 --selection-mode supported --cta-definition both \
  --exclude-gene-pattern 'MAGE*' --allow-gene MAGEA4 \
  --ms-support-mode sample-affinity --ms-affinity-nm 1000 --allow-untyped-ms \
  --include-utrs --poly-a-length 120 --max-length-aa 686 --max-length-nt 2500 \
  --padding-step 2 --beam-width 6 --optimization-rounds 10 -o magea4-ten-out
```

For the published RNA designs, the complete construct budget is **2500 nt**.
The HBB 5′ UTR is 50 nt, HBB_FI 3′ UTR is 268 nt, polyA is 120 nt and the
stop codon is 3 nt. The encoded protein budget is therefore
`floor((2500 - 50 - 268 - 120 - 3) / 3) = 686 aa`, including the initiating
methionine and all linkers. A full 686-aa construct uses 2499 total nt.
Custom UTR/polyA settings change this calculation; `--max-length-nt` always
constrains the full nucleotide sequence.

The vaccine extra needs OncoRef >=1.8.207 for loose CTAs. Current Vaxrank and
PirlyGenes releases pin 1.8.206, so use a compatible separate environment
until coordinated validation ([Vaxrank #579](https://github.com/openvax/vaxrank/issues/579),
[PirlyGenes #632](https://github.com/pirl-unc/pirlygenes/issues/632)). The base
Tsarina installation retains its existing reference minimum. OncoRef's
genome-scope expression calculation can need complete per-sample matrices;
`--auto-fetch` permits their download. Missing references fail explicitly.

Live vaccine design also requires **Hitlist >=1.64.7**, including the merged
correctness fixes [#637](https://github.com/pirl-unc/hitlist/pull/637) and
performance/default-CTA changes [#641](https://github.com/pirl-unc/hitlist/pull/641).
The command checks the version of the actually imported code. Both that
version and installed distribution metadata are recorded; Hitlist owns the
freshness/rebuild check for observation and mapping artifacts.

## Selection and expression units

```
prevalence(p,h) = fraction of cohort h with p in each sample's top 5% expression ranks
prevalence(p,c) = sum_h(n_h * prevalence(p,h)) / sum_h(n_h)
score(p)       = sum_c(world_mortality_pct(c) / 100 * prevalence(p,c))
```

The reference axis is OncoRef's **genome-scope collapsed transcriptome**.
Genes encoding an identical canonical protein have their TPMs summed within
each patient **before** expression ranking. Noncoding singleton genes remain
on this transcriptome axis. This is an RNA-expression proxy, not measured
protein abundance or the 95th percentile of TPM across patients.

The genome registry includes identical CTA proteins omitted from the focused
registry, such as RBMY1F/J ([OncoRef #567](https://github.com/pirl-unc/oncoref/issues/567)).
Its rank axis therefore differs from older CTA-scope plots. The loader checks
full canonical sequence identity against expression keys and refuses a
disagreement. Gene prevalence fractions cannot be added to repair one.

Top-k ranks individual additive scores, with deterministic name-based ties.
It does not optimize distinct-patient union coverage or estimate preventable
deaths. Each mapped cancer category's mortality share is counted once.
Sample-weighted cohorts are an empirical mixture, not worldwide patient
prevalence; represented histologies can be narrower than the burden category.
Missing measurements stay missing in the cohort table. Incomplete scores
are observed partial scores and are flagged for cautious comparison.

`--selection-mode ranked` preserves top-k candidate selection before downstream
filtering. `--selection-mode supported` instead inspects eligible positive-score
proteoforms in rank order and requires k targets that contribute native sequence
after specificity, panel MS support and minimum whole-segment length checks.
It reserves one shortest ligand-preserving segment per target before filling
extra pieces by rank. A shortfall writes an audit and fails; it never silently
returns fewer contributing proteins. `selection_screen.csv` records the complete
ranking with exclusion, unsupported, length-rejected, selected and uninspected
states. A bounded scan batch can inspect lower-ranked targets beyond the selected
k; those remain auditable.

`--selection-mode budget` removes the protein count cap and requires an amino
acid or total nucleotide limit. It screens every eligible positive-scoring
proteoform and greedily adds whole ligand-bearing pieces by marginal
mortality-weighted lower expression-union gain per amino acid. Ties use
incidence gain, new supported HLA alleles and new exact observed peptides.
When these gains per amino acid are equal, prefer another piece from an
already selected protein, then a longer contiguous native piece, then protein
rank. This reduces additional proteins and fragmentation without overriding
higher coverage/evidence scores or restoring excluded sequence.
The expression surrogate is `max(prevalence)` per cancer;
additional pieces from the same protein can add ligand evidence without
counting expression twice. Missing measurements add no known expression gain.
The allocation is a heuristic, with its decisions in `budget_allocation.csv`;
it does not optimize actual patient overlap or guarantee a global optimum.
`--top-k` is ignored in this mode.

This lexicographic priority lets a very small positive expression gain outrank
many additional MS peptides. Its HLA gain counts new alleles across the entire
construct, so it can miss broader presentation of PRAME when those alleles
already support another protein. The [coverage and evidence comparisons](vaccine-evidence-tradeoffs.md)
quantify these limitations and compare alternative allocations without changing
the current command's default scoring.

Repeatable `--exclude-gene-pattern` accepts gene-symbol globs;
`--allow-gene` provides exact symbol exceptions. Any excluded member vetoes an
identical-full-sequence group. These controls affect candidate eligibility,
not CTA membership, expression grouping or the non-CTA background. MAGEA4 is
the target of FDA-approved [TECELRA](https://www.fda.gov/vaccines-blood-biologics/cellular-gene-therapy-products/tecelra);
that approval does not establish this vaccine's safety or approve its other
peptide-HLA assignments.

Global cancer incidence, mortality, absolute counts, CTA p95 prevalence,
and sample denominators are separate table columns. Incidence is not an
additional multiplicative score factor. OncoRef's current curated reference
uses GLOBOCAN 2022 shares; missing absolute counts are never synthesized,
and rounded shares are preserved without renormalization. Source/scope
follow-ups: [OncoRef #542](https://github.com/pirl-unc/oncoref/issues/542),
[#543](https://github.com/pirl-unc/oncoref/issues/543).

The exact default mapping is in `DEFAULT_CANCER_COHORTS` and each manifest.
Override it with `--cancer-cohorts mapping.json`, containing disjoint lists:

```json
{"lung": ["LUAD", "LUSC"], "colorectal": ["COAD", "READ"]}
```

No residual mortality category or nested subtype is counted automatically a
second time. Requested references must be available or explicitly fetched.

## Definitions and the native sequence funnel

Strict uses OncoRef's canonical core reproductive scope: **testis, ovary and
placenta**. Loose adds **cervix, endometrium, epididymis, fallopian tube,
prostate, seminal vesicle and vagina** to the RNA numerator. Neither adds
breast or thymus. Thymus remains in the default RNA fraction denominator but
is excluded from somatic maxima and the reproductive protein flag; protein
annotations use broader reproductive tissue conventions. These are curated
RNA/protein restriction rules with exceptions, not absolute absence gates.
It does not admit every nominated or warning-tier CTA. Each definition has
its own non-CTA background and independent design under `strict/` or `loose/`.

1. Collapse byte-identical longest translated proteins, preserving all gene,
   transcript and protein IDs. Ranked mode selects top-k before downstream
   filtering; supported mode explicitly backfills contributing targets.
   For specificity, same-symbol Ensembl `HSCHR` alternate-haplotype annotations
   inherit the curated primary gene's CTA membership; they are another
   annotation of that gene, not an independent non-CTA locus. The exact
   resolutions are in `background_cta_aliases.csv`
   ([Tsarina #187](https://github.com/pirl-unc/tsarina/issues/187)). This does
   not change expression identities, sum extra prevalence, or admit ordinary
   non-CTA loci by symbol/sequence similarity.
2. Screen every translated coding non-CTA isoform, including IG/TR germline
   coding segments. Remove the union of every residue covered by a shared
   8-mer; break at ambiguous residues. Retain complementary native stretches.
   Sharing with another CTA is allowed, including an unselected CTA; sharing
   with any independent non-CTA source is disqualifying, including an identical
   full protein.
3. Query live Hitlist human class-I observations for exact peptides entirely
   within those stretches. Require positive MS modality; explicitly non-MS
   fluorescence, stability and structural assays and unknown-modality records
   are rejected with complete metadata in `rejected_ms_observations.csv`.
   Explicitly negative assay results are also rejected; an MS method alone
   does not establish detection. Separate positive observations remain eligible.
   Hitlist's nonbinding flag alone is insufficient ([#644](https://github.com/pirl-unc/hitlist/issues/644)).
   Curated MS-only supplements may have blank methods; supplied data with blank
   methods must declare `assay_modality=mass_spectrometry`.
   Default presentation mode uses monoallelic percentile <=2, sample/deconvolved
   <=1 and unrestricted peptide MS plus prediction <=0.5. Sample-affinity mode
   instead accepts every typed sample panel allele with affinity below
   `--ms-affinity-nm` (default 1000), without requiring a best-allele assignment
   or presentation-percentile cutoff. `--allow-untyped-ms` permits prediction
   against the panel when exact sample typing is unavailable; study-wide allele
   pools are not sample genotypes. `ms_assignments.csv` links observations to
   qualified alleles. The peptide must be exactly observed; a longer observed
   peptide does not support a different unobserved nested epitope. Inferred
   allele support is always distinct from measured restriction.
4. Keep stretches containing at least one qualifying ligand, preserving all
   supported ligands and repeated occurrences. Trim unsupported terminal
   padding to at most 10 aa by default; internal native context stays intact.
5. Fit whole pieces in protein rank order under construct limits, prioritizing
   broader allele support within a protein. Record all exclusions; never cut
   a supported ligand to force a piece to fit. Supported mode first reserves one
   shortest whole ligand-bearing segment per selected target.

`funnel.csv` reports raw, specific, MS-supported, maximally padded, and
assembled aa lengths/piece counts, fractions retained, support and dropout
reasons. MS-supported length means the full specific interval before end
trimming. Native/API/CSV coordinates are **zero-based, half-open**.
`assembled_ms_ligand_count` and `assembled_pmhc_count` count distinct retained
peptides and peptide-HLA pairs per proteoform. `hla_support_counts.csv` counts
retained peptides by allele and evidence tier. Shared peptides across different
proteoforms and multi-allele predictions must not be counted as independent MS
observations or patients.

### Source-verified nonmalignant tissue exclusion

By default, `--normal-ms-policy audit` displays normal-tissue observations.
Use `--normal-ms-policy exclude --normal-ms-atlas-dir PATH` to also remove
residues covered by 8-mers from donor-resolved nonmalignant **heart, brain and
lung HLA-I** observations. `PATH` contains the original HLA Ligand Atlas
2020.12 `peptides`, `donors` and `sample_hits` TSVs, as `.tsv.gz` files or the
release ZIP's `HLA_*.tsv` names. Hitlist >=1.65.1 reads and fingerprints them.
The vaccine default threshold is **one distinct donor**, configurable with
`--normal-ms-min-donors`; repeated samples from one donor are not independent
replication. All three organs fall outside both CTA definitions.

The [Atlas data](https://hla-ligand-atlas.org/data) and
[primary study](https://doi.org/10.1136/jitc-2020-002071) establish primary
autopsy tissue from donors without diagnosed malignancy. Donors could have
other diseases. Cell lines, tumor-adjacent sources and unresolved donor
records do not establish this exclusion. IEDB tissue flags alone are audit
evidence. The blacklist does not require predicted affinity to the vaccine
panel, and donor HLA genotypes are not measured peptide restrictions.

`normal_ms_summary.csv`, `normal_ms_audit.csv` and
`normal_ms_exclusions.csv` preserve exact source records, donor counts and
qualification decisions. `cta_specific_before_normal_ms.csv` and the added
funnel stage distinguish sequence lost to the normal-MS gate. Other normal
observations remain visible warnings. Neither absence of an observation nor
this exclusion establishes vaccine safety.

## Construct search and design elements

The deterministic beam search compares complete constructs using adjacent
swaps, reversal/rotations, N/C padding changes, and direct/AAY joins.
Equivalent padding choices at native boundaries share one search state.
At a length cap, linker additions can jointly trim adjacent terminal padding,
preserving every retained ligand without requiring a worse intermediate join.
`--min-padding`, `--max-padding`, `--padding-step`, `--beam-width`, and
`--optimization-rounds` control its search budget. This heuristic does not
guarantee a global optimum.

The lexicographic objective minimizes junction-window/allele affinity
predictions below 1000 nM (`--junction-affinity-nm`), then their log affinity
burden, then maximizes mean predicted boundary cleavage, then minimizes
linker length and retains more native context. The final translated product
is audited, including initiating M, linker interiors and both linker edges,
and windows crossing multiple boundaries. All scores are retained; missing
or invalid predictions fail. `--require-clean-junctions` fails after writing
the audit if binders remain.

Pepsickle's human-only in-vivo model predicts boundary cleavage through
mhctools: eight preceding residues, the cleavage residue itself, and eight
following residues where available. Its identity is in the manifest.
Cleavage and binding predictions do not establish processing, presentation,
immunogenicity or absence of cross-reactivity.

UTR names/sequences and command concepts follow Vaxrank. `--include-utrs`
defaults to HBB 5' and tandem HBB_FI 3'. Custom A/C/G/T/U strings or `none`
are accepted. `--poly-a-length` defaults to 0. DnaChisel codon optimization
uses `use_best_codon` with exact translation enforcement and `h_sapiens`
by default. No signal peptide or trafficking domain is added implicitly.
RNA/mRNA are aliases and use U; DNA describes an insert without promoter or
plasmid backbone. FASTA does not encode a chemical cap or modified nucleosides.

The aa limit includes initiation M and linkers. Total nt includes UTRs,
CDS/stop and polyA. `--max-construct-length-aa/nt` alias `--max-length-aa/nt`.

Sources: [Vaxrank mRNA library](https://github.com/openvax/vaxrank/blob/main/vaxrank/mrna_library.py),
[RefSeq NM_000518.5](https://www.ncbi.nlm.nih.gov/nuccore/NM_000518.5),
[Sahin et al.](https://doi.org/10.1038/nature23003),
[Holtkamp et al.](https://doi.org/10.1182/blood-2006-04-015024),
[MHCflurry 2.0](https://pubmed.ncbi.nlm.nih.gov/32711842/),
[Pepsickle](https://pubmed.ncbi.nlm.nih.gov/34478497/).

## Output bundle and API

Every successful construct writes a self-contained `website/` alongside the
scientific report. With `--definition both`, a combined website also appears
at the output root. View the [published Vaccine Atlas](vaccine-results/index.html)
for strict/loose budget designs and ten-protein comparisons.

| Artifact | Meaning |
| --- | --- |
| `report.md`, SVG figures | Selected proteins/cancers, retention funnel, construct map |
| `ranking.csv` | All candidates, identities, scores and selection flags |
| `cancer_cohorts.csv`, `cancer_summary.csv` | Incidence, mortality, p95, denominators, missingness and score contributions |
| `*_reference.csv`, `expression_availability.csv` | Upstream definitions and source provenance |
| `background_cta_aliases.csv` | Audited same-gene alternate-haplotype identity resolutions |
| `specific_intervals.csv`, `ligands.csv`, `funnel.csv` | Native regions, MS support, per-protein retention |
| `hla_panel.csv` | Existing Tsarina frequency/provenance audit for requested alleles |
| `layers.csv`, `excluded_segments.csv` | Source-to-construct coordinates and whole-piece exclusions |
| `junctions.csv`, `cleavage.csv`, `search_history.csv` | Final predictions and optimization history |
| `protein.fasta`, `cds.fasta`, `full.fasta` | Single antigen, CDS/stop, complete DNA/RNA |
| `manifest.json` | Configuration/results, versions, model identities, input/sequence/output hashes |
| `website/index.html` | Interactive protein/cancer/HLA results, definitions and sources |
| `website/downloads/*/cumulative-*.csv` | Protein/segment cumulative coverage bounds and MS counts |
| `protein_ms_map.csv`, `full_protein_ms_observations.csv` | Exact native peptide positions and sample context, including removed regions |

The website includes SVG/PNG figures for HLA frequencies and measured/inferred
support, protein additions, MS peptides versus final construct prefix length,
and full-protein tissue maps. Each observed peptide is counted once in overall
MS totals, with its strongest HLA evidence tier; these are not T-cell-validated
epitopes. Large website CSV downloads are losslessly gzip-compressed.

HLA carrier reach sums published CIWD Table A2 allele frequencies within each
locus, applies Hardy–Weinberg equilibrium to two copies, then assumes linkage
equilibrium across loci. Missing CIWD frequencies are explicit and omitted
from the known-allele proxy. Frequencies are never panel-normalized. This donor
reference is not census-weighted global population coverage. Cancer curves
show marginal union bounds `max(p)` to `min(1, sum(p))`; missing measurements
widen the upper bound. Mortality/incidence use separate global burden shares
without renormalizing the represented cancer subset. HLA and cancer curves
are not multiplied into a clinical coverage estimate.

Render any named collection of completed reports without rerunning models:

```sh
tsarina vaccine-report \
  --report strict=designs/strict --report loose=designs/loose \
  --output-dir vaccine-website --analysis-date 2026-10-07
python -m http.server --directory vaccine-website 8000
```

The renderer verifies every saved manifest artifact hash before reading it.
Use a fresh website output directory to avoid mixing results.
Reports can have different panels, DNA/RNA settings and target counts. Relative
assets require no CDN or scientific software to view once served over HTTP.

Regional HLA frequencies are proxy evidence, not guaranteed population
coverage. The broad default `global54_abc` can be replaced by a named panel
or an exact HLA-A/B/C list. Use a fresh output directory for each design.

```python
from tsarina.vaccine import design_vaccine
from tsarina.vaccine_construct import VaccineConfig

result = design_vaccine(
    VaccineConfig(top_k=10, definition="strict", panel="global54_abc"),
    auto_fetch=True, output_dir="vaccine-audit",
)
```

`VaccineInputs` supports versioned offline data, with prevalence already
collapsed and ranked correctly. Injected prediction callbacks are explicitly
labeled in the manifest; deterministic example scores are not evidence.

Run the fully labeled synthetic example without reference downloads:

```sh
python examples/cta_vaccine_demo.py /tmp/cta-vaccine-demo
```

It generates every output table and figure, demonstrating identical-protein
collapse, cancer ranking, shared-end subtraction, allele-specific MS evidence,
construct coordinates and DNA/RNA encoding. Its illustrative values are
explicitly marked in the report and manifest.
