# Canine cancer-antigen evidence and DLA design

`tsarina vaccine --species canine` imports a frozen dog evidence bundle, ranks
exact protein sequences by tumor RNA prevalence in independent dogs, and retains
contiguous native regions with exact positive-MS and DLA affinity support.
Outputs include source tables, peptide assignments, sequence-retention funnels,
native-protein maps, cumulative coverage/MS figures and an offline HTML report.

The canine mode does not fetch human CTAs, HPA expression, TCGA p95 or human
mortality priors. It implements the offline migration in
[#208](https://github.com/pirl-unc/tsarina/issues/208). **No real canine target
list or validated vaccine is supplied.** The runnable example uses invented
proteins, dogs, evidence and affinity values to demonstrate the contracts.

## [Offline example](canine-example/index.html)

Install `tsarina[vaccine]`. From a Tsarina checkout:

```sh
tsarina vaccine --species canine \
  --input-bundle tests/fixtures/canine-vaccine.json \
  --canine-cohort osteosarcoma --allow-exploratory-dla \
  --selection-mode supported --top-k 2 --lengths 8,9 \
  --min-padding 0 --max-padding 0 --optimization-rounds 1 \
  --vaccine-type rna --poly-a-length 6 --max-length-nt 90 \
  -o canine-out

# Copy verified saved outputs without executing any model.
tsarina vaccine-report --report canine=canine-out -o canine-website
```

Open `canine-out/index.html`. Omit `--allow-exploratory-dla` to produce the RNA
evidence report while excluding the example's extrapolated DLA support. A bundle
can declare an empty panel and fully missing genotype data for an RNA-only
report. Unknown normal assessments cannot admit targets.

The example assembles two protein groups. A third group is rejected because an
independent adult-heart locus encodes the identical complete protein. A longer
verified healthy-heart MS peptide removes shared 8-mers from the first protein;
its exact ligand remains in an eight-residue region. Three invented dogs, one
with repeated biopsies, supply the prevalence denominator. An unresolved donor
is audited separately. One genotype contains three distinct DLA-88 alleles;
unsupported alleles and 20% missing genotype mass remain in the coverage bounds.
The output is **exploratory and unassessed**, even when evaluated junctions pass.

## Frozen bundle contract

The schema identifier is `tsarina.canine.v1`. See the
[complete fixture](https://github.com/pirl-unc/tsarina/blob/main/tests/fixtures/canine-vaccine.json)
and [importer](https://github.com/pirl-unc/tsarina/blob/main/tsarina/canine_inputs.py).

| Field | Required meaning |
| --- | --- |
| `taxon` | `9615`; reference and reference-scoped records must agree |
| `reference`, `reference_key` | Assembly accession, annotation/source version and assembly/annotation asset SHA256s; key is canonical-JSON SHA256 |
| `background_complete` | Explicit declaration of a complete translated background for this reference; a target-only FASTA cannot establish specificity |
| `sources` | Source IDs with HTTPS URL, version, license and asset SHA256; every evidence row names an inventoried source |
| `occurrences` | Every complete protein occurrence, gene/transcript/protein IDs, exact full-sequence SHA256, restriction status/reason and normal assessment |
| `restriction_policy` | Name/version, `strict` or `loose`, explicit allowed normal tissues, inventoried normal-assessment source, and a separate `somatic_tpm_threshold` |
| `expression_policy` | Positive absolute tumor-expression threshold and `TPM` unit |
| `cohorts`, `samples` | Named histology/quantification and population frame; study, reconciled donor, specimen, preparation, breed, treatment, pooling, QC and unit |
| `rna_bounds` | One lower/upper abundance interval per sample and exact protein group, with allocation provenance; duplicate/group-ID/unit mismatches fail |
| `ms_hits` | Exact peptide, raw method/response/result/context, reviewed boolean `is_ms_observation` and admission reason, MHC class, observation/sample IDs, independent host/peptide-source/MHC taxa, sample kind/health/tissue, restriction kind and actual sample alleles |
| `capabilities` | Per-allele tier, reason, executable lengths, model/version, weights and MHC-sequence hashes, validation source; empirical support also requires canine validation taxon/hash |
| `panel` | Explicit canonical DLA class-I allele list; `--alleles` can override it |
| `predictions` | Frozen peptide–DLA `affinity_nm` values; missing required native or junction entries fail |
| `prediction_capability_sha256` | Canonical hash binding a nonempty affinity archive to its declared capability inventory |
| `genotype_population` | Named sampling frame for the explicitly selected RNA cohort, and unobserved `missing_mass` in `[0,1]` |
| `tumor_genotype_pairs` | Independent observed matching dogs or explicitly simulated cohort pairings; full allele lists, breed/source, positive weights |

Canonical JSON uses sorted keys, compact separators, UTF-8 and no NaN. Protein
identity preserves I/L and requires the complete uppercase sequence without a
terminal stop. Asset hashes bind declared provenance; the importer does not
download or independently reconstruct those assets. The input copy, source
tables and manifest preserve all rejected and unknown occurrences and decisions.

Reviewed sequence-group RNA bounds must come from an upstream allocation
analysis. The importer does not infer transcript abundance from gene totals,
add gene prevalence fractions or rescue candidates by human orthology. It
recomputes exact full-protein groups and independent-dog prevalence, and vetoes
contradictory healthy-primary somatic RNA. All source occurrences, including
unrelated same-symbol loci, remain in the normal-sharing search.

## Prevalence and coverage

Only QC-passing, unpooled, untreated primary bulk tumors enter the named cohort.
Cell lines, treated tumors, unresolved donor IDs and pooling/QC failures remain
in `cohort_sample_audit.csv`. Histology and quantification must agree with the
cohort; study, breed and other sample metadata remain available for stratified
analysis. Reconciliation of donor aliases belongs upstream.

For each exact protein and independent dog, the expression lower bound is one
only when every biopsy has a measured lower TPM above the threshold. The upper
bound is one if any biopsy could pass or any measurement is missing. Prevalence
is the sum of these indicators divided by **all identified admitted dogs**;
missing RNA is not zero. This is absolute-expression evidence in the supplied
cohort, not a worldwide estimate, TCGA rank prevalence or a confidence interval.

`paired_target_coverage` in the shared coverage module uses each selected
protein's RNA bounds and its own qualified DLA alleles on explicit full
genotypes. Adding a presenter to PRAME-like target A can improve A's joint reach
even if that allele was already represented by target B. Observed matching and
simulated independent-cohort pairings are labelled separately. Unknown
expression, unsupported genotypes and missing population mass widen bounds;
allele frequencies are neither invented nor renormalized. No HWE, cross-locus
independence or two-copy assumption is imposed.

Unsupported genotype mass counts dogs carrying any unassessed allele. Such a dog
can also have demonstrated reach through another allele, so this mass overlaps
the coverage bounds and cannot be added to them as a separate category.

Coverage figures describe final-order construct prefixes, not independently
optimized shorter constructs or clinical protection. MS curves count distinct
observed peptide strings, split into observed monoallelic restriction and
inferred-only restriction; assay observations and peptide–DLA pairs have separate
counts. DLA bar charts show peptide evidence, **not allele frequencies**.

## Restriction policy, exact MS and assembly

There is no default canine reproductive-tissue policy. `strict` and `loose` are
labels for the frozen bundle's explicit allowed tissues and normal thresholds.
Comparing definitions requires two independently reviewed bundles; `--cta-definition
both` cannot manufacture them. A loose bundle must be requested explicitly with
`--cta-definition loose`. The human prostate/epididymis allowances do not transfer.

The design subtracts every 8-mer in a rejected or unknown reference occurrence,
then every 8-mer of positive, verified healthy-primary dog MS outside the allowed
tissues. A `verified_normal` MS row needs a `normal_verification_source` tracing
its donor/specimen; tumor-adjacent samples and cell lines are not healthy-primary
counterevidence. Longer normal peptides also exclude shared shorter windows.
Native maps distinguish exact observed ligands from healthy-MS 8-mer overlaps.

An exact positive-MS peptide must bind below `--ms-affinity-nm` to a supplied
sample allele in the requested panel. Its context must be class-I MHC ligand
elution; bulk proteomics and class-II assays remain rejected audit records.
The reviewed source's positive-MS decision takes precedence over raw binding
qualitative labels; raw methods/responses and the admission reason remain visible.
This follows Hitlist's positive-MS contract without reinterpreting its historical
index partition or using nonbinding as evidence of MS.
Monoallelic assignments remain separate
from multiallelic affinity inference. `--allow-untyped-ms` permits panel-affinity
inference for untyped samples. An observed longer peptide never supports a
different unobserved nested peptide. Human-host DLA transfection is preserved as
a separate context from endogenous canine MS.

Assembly reuses Tsarina's native interval/segment objects, padding, linker and
bounded junction search. The entire product is checked for background 8-mers,
including added methionine and linker-only windows; such overlap rejects the
construct. `ranked` screens top-k candidates; `supported` backfills to k
contributors subject to MS and whole-segment length constraints. A shortfall is
reported and the CLI exits nonzero. Canine `budget` allocation is not yet
implemented: its species-scoped marginal objective belongs with
[#201](https://github.com/pirl-unc/tsarina/issues/201) and
[#202](https://github.com/pirl-unc/tsarina/issues/202).

RNA/DNA, AA/total-nt caps, stop codon and polyA constraints apply. Canine coding
defaults to deterministic standard-code reverse translation, labelled `generic`
and unoptimized. The installed codon-table inventory has no canine table.
Other explicit codon species remain recorded choices; `h_sapiens` is rejected as
a canine default. UTRs require explicit sequences or `none`; human HBB presets
are not implicit canine elements. Human-only Pepsickle is not called; absent
canine cleavage assessment contributes no search preference.

## Data integration and prediction limits

Frozen affinity archives and explicit Python affinity callbacks are supported.
The canine CLI does not implicitly run a human presentation-percentile model.
Every selected allele must declare capability for every requested peptide length.
Sequence-extrapolated selection requires `--allow-exploratory-dla`; final canine
presentation/processing remains unassessed. `--require-clean-junctions` cannot
certify a canine construct. Existing fixed-membership/search limits remain in
[#191](https://github.com/pirl-unc/tsarina/issues/191) and
[#192](https://github.com/pirl-unc/tsarina/issues/192).

Canvax's current real-data inventory includes NCBI canine coding sequences,
Bgee/BarkBase expression screens, canine tumor RNA inventories and curated DLA
ligands. It has not produced an admitted real sequence-level CTA/prevalence
bundle: osteosarcoma units and annotation/donor crosswalks require reconciliation,
and melanoma cell lines cannot provide independent primary-dog prevalence.
This importer supplies the downstream contract rather than inventing that table.

Upstream dependencies and limits:

| Owner | Tracked work |
| --- | --- |
| [OncoRef #571](https://github.com/pirl-unc/oncoref/issues/571) | Species/reference-scoped exact groups preserving every source occurrence |
| [Hitlist #660](https://github.com/pirl-unc/hitlist/issues/660) | Reviewed offline canine/DLA supplements with host/source/MHC attribution |
| [Hitlist #661](https://github.com/pirl-unc/hitlist/issues/661) | Canine reference and normal-policy evidence bundles |
| [mhctools #544](https://github.com/openvax/mhctools/issues/544) | Actual DLA model-key mismatch; no human-allele substitution |
| [MHCflurry #490](https://github.com/openvax/mhcflurry/issues/490) | Release-bound capability and training-disjoint canine evaluation |
| [Vaxrank #596](https://github.com/openvax/vaxrank/issues/596) | Canine evidence admission through native antigen/construct objects |
| [Tsarina #209](https://github.com/pirl-unc/tsarina/issues/209) | Remaining human vaccine gate admits unknown MS outcomes; canine policy requires positive results |
| [Tsarina #210](https://github.com/pirl-unc/tsarina/issues/210) | Human assembled-product background check; canine constructor checks initiation, joins and linkers |

The [2026 canine/human tumor-antigen ligand study](https://pubmed.ncbi.nlm.nih.gov/42199926/)
includes different host and restriction contexts; those facets must survive
import. [Canine class-I genotyping](https://pmc.ncbi.nlm.nih.gov/articles/PMC9926114/)
and [breed-cohort DLA haplotypes](https://pmc.ncbi.nlm.nih.gov/articles/PMC10001263/)
support explicit genotype/cohort treatment rather than renaming HLA loci.
[NetMHCpan 4.2](https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/)
includes dog EL data; canine support is not absent from all existing models.
Model execution, training inclusion and independent validation are separate
claims, assessed per allele, length and release.

The portable Python route uses the same checks:

```python
from tsarina.canine_inputs import load_canine_vaccine_inputs
from tsarina.vaccine import design_vaccine
from tsarina.vaccine_construct import VaccineConfig

inputs = load_canine_vaccine_inputs("reviewed-canine-bundle.json")
config = VaccineConfig.for_canine(
    "osteosarcoma", top_k=10, selection_mode="supported",
    allow_exploratory_dla=True, max_length_nt=2500,
)
result = design_vaccine(config, inputs, output_dir="canine-out")
# An explicit affinity_fn(peptides, alleles) can supply a runtime backend.
# Its returned DataFrame needs peptide, canonical allele and affinity_nm.
```

Human command defaults, selection formulas and saved scientific examples are
preserved. The canine extension does not change the existing human coverage
proxy or fix its separate allocation limitations.
