# CTA panel design

The panel workflow designs an off-the-shelf set of CTA peptide-MHC targets
across a population HLA panel. It selects biologically suitable CTAs, finds
evidence-supported peptides for each allele, and reports both the resulting
matrix and estimated HLA coverage.

Use this workflow for cohort-level vaccine or TCR-program design. For one
patient with known tumor findings, use
[Personalized target selection](personalized-targets.md) instead.

## Quick start

Build the default panel as a readable report:

```bash
tsarina panel
```

Write one row per selected pMHC for analysis:

```bash
tsarina panel --format long --output panel-long.csv
```

Write a CTA × HLA matrix whose cells contain selected peptides:

```bash
tsarina panel --format wide --output panel-wide.csv
```

The default run targets up to 25 non-empty CTAs across `global54_abc`, uses
8–11mer CTA-exclusive peptides, requires public-MS evidence, scores with
MHCflurry, and retains up to three peptides per CTA × HLA cell.

## Selection pipeline

Panel construction is intentionally staged. Early stages define and screen CTA
sources; later stages select pMHCs without changing CTA definitions.

### 1. Rank candidate CTAs

Automatic selection begins with oncoref's canonical CTA universe. The default
rank, `tumor_prevalence_panel_score`, combines:

- HPA cancer RNA prevalence at pTPM ≥ 2.0;
- cancer-type breadth using a 5% within-type sample-prevalence floor; and
- HPA cancer IHC as a weak tie-breaker.

These bundled prevalence features establish the scan order. Live public-MS
support and safety are recalculated from the current observations index in a
later stage.

Change the initial ranking with:

```bash
tsarina panel \
  --cta-rank-by tumor_prevalence_panel_score \
  --cancer-rna-threshold 2.0 \
  --cancer-type-prevalence-floor 0.05
```

### 2. Apply source-level safety gates

Automatic selection keeps `HIGH` and `MODERATE` oncoref
`restriction_confidence` by default. It also excludes:

- CTAs with RNA above 2.0 nTPM in brain/CNS/cerebellum, heart, lung, liver, or
  pancreas;
- CTAs with peptide evidence that maps uniquely to that CTA in those healthy
  vital tissues; and
- MAGE-family CTAs other than MAGEA4.

`PRAME`, `CTAG1A/CTAG1B`, and `MAGEA4` are allowlisted by default. The
allowlist is applied explicitly and its members are placed ahead of
lower-ranked automatic candidates.

Explicit `--ctas` requests bypass the automatic family and source-selection
gates. Aliases such as `NY-ESO-1` and `MAGE-A4` are normalized.

### 3. Enumerate and group peptides

The workflow enumerates 8–11mer CTA peptides and removes peptides that also
occur in non-CTA human proteins. A single displayed target can expand to more
than one Ensembl gene: CTAG1A and CTAG1B both resolve to one target, reported
as `NY-ESO-1`.

Identical-protein groups are reported under oncoref's preferred symbol — a
curated alias where one exists (`NY-ESO-1`), else the prefix-contracted
members (`XAGE1A/B`, `SSX2/B`, `GAGE12C/D/E`) — which is the same naming
[`tsarina personalize`](personalized-targets.md) uses. Long format keeps the
full membership in `cta_members`, and `panel_summary()["cta_groups"]` lists
it for either format. Every one of these symbols is also accepted as `--ctas`
input, as are the member symbols and the full members label.

Pass `--proteoform-labels members` to report the full members label
(`CTAG1A/CTAG1B`) instead. Selection, filtering, ranking and grouping are
identical either way; only the reported label changes.

Two forms of redundancy are grouped by default:

- CTAs with identical enumerated peptide sets; and
- CTAs with identical final selected pMHC panels.

This prevents paralogs from consuming multiple automatic panel slots.

### 4. Add MS evidence and score HLA presentation

For each candidate batch, Tsarina queries the current hitlist observations
index. It separates monoallelic, sample-allele/deconvolved, unrestricted, and
prediction-only evidence, then applies an evidence-specific presentation
percentile cutoff.

Within each CTA × HLA cell, candidates are ordered by MS source count, MS hit
count, prediction percentile, and deterministic tie-breakers.

### 5. Backfill non-empty targets

The default `--cta-count 25` means up to 25 downstream CTAs with at least one
selected pMHC. When a highly ranked CTA becomes empty after peptide,
exclusivity, MS, or prediction gates, Tsarina scans lower-ranked candidates to
backfill the slot.

Use `--show-empty-ctas` to inspect the original top-ranked candidates,
including downstream failures. Explicit `--ctas` requests are always retained
in the audit view, even when they produce no selected pMHCs.

## Evidence tiers

Lower presentation percentiles are better. The default pMHC evidence tiers
are:

| Evidence tier | Default cutoff | Interpretation |
|---|---:|---|
| `monoallelic_ms` | ≤ 2.0 | Peptide observed in monoallelic MS for the selected HLA |
| `sample_allele_ms` | ≤ 1.0 | Multi-allelic MS with an exact restriction or donor allele set; the selected HLA is best predicted within that set |
| `unrestricted_ms` | ≤ 0.5 | Class-I MS evidence without a usable allele assignment; prediction assigns the panel allele |
| `predicted_only` | ≤ 0.1 | No MS support; disabled unless explicitly requested |

Rows reported only as `HLA class I` can enter `sample_allele_ms` when a usable
donor allele set is available. Otherwise they use the stricter unrestricted
tier.

Tune the cutoffs with:

```bash
tsarina panel \
  --monoallelic-ms-max-percentile 2.0 \
  --sample-allele-ms-max-percentile 1.0 \
  --unrestricted-ms-max-percentile 0.5
```

Prediction-only candidates are opt-in:

```bash
tsarina panel \
  --include-predicted-only \
  --predicted-only-max-percentile 0.1
```

## Output and progress

| Format | Intended use |
|---|---|
| `table` | Human-readable pMHC report plus coverage summary |
| `long` | One row per selected peptide-HLA pair with evidence and score provenance |
| `wide` | CTA rows × HLA columns, with selected peptides in each cell |

The table summary orders CTAs by selected peptide count, HLA-hit count, and
estimated population coverage. It reports monoallelic MS support separately
from sample-genotype/deconvolved support.

Progress messages cover peptide enumeration, MS loading, scoring, evidence-tier
construction, and final selection. Interactive terminals also show a scoring
bar. Control this behavior with:

```bash
tsarina panel --progress off
tsarina panel --progress on
```

MHCflurry scores all alleles in one batch by default because chunking repeats
allele-independent processing-model work. Use `--score-chunk-size` only when
that tradeoff is desirable.

## HLA panels

| Panel | Alleles | Purpose |
|---|---:|---|
| `iedb27_ab` | 27 | IEDB global baseline for HLA-A/B |
| `iedb36_abc` | 36 | IEDB baseline extended with HLA-C |
| `global44_abc` | 44 | Adds representation for East Asia, South Asia, and Sub-Saharan Africa |
| `global48_abc` | 48 | Adds representation for Latin America and MENA |
| `global51_abc_ssa` | 51 | Legacy Global-48 extension for Sub-Saharan Africa |
| `global51_abc` | 51 | Reference A/B backbone plus frequent HLA-C and common-A/B complements |
| `global53_abc` | 53 | Legacy CTA-MS extension; fixed membership omits C*14:03 |
| `global54_abc` | 54 | Default Global-51 extension with CTA-MS-supported alleles, including C*14:03 |

Use a named panel or provide an explicit list:

```bash
tsarina panel --panel iedb27_ab
tsarina panel --alleles 'HLA-A*02:01,HLA-A*24:02,HLA-B*07:02'
```

### Why Global-54 is the default

`global51_abc` contains:

- the 27 IEDB/TepiTool class-I HLA-A/B reference alleles;
- 21 frequent HLA-C allotypes from the Sarkizova HLA-C peptidome coverage set;
  and
- `B*18:01`, `B*40:02`, and `B*46:01`, the highest-frequency calibrated
  alleles needed to complement the IEDB/Paul common HLA-A/B set.

`global54_abc` adds `A*29:02`, `B*15:02`, and `B*27:05`, which were the top
missing alleles in a public CTA-MS audit. It retains all 21 reference HLA-C
allotypes, including both `C*14:02` and `C*14:03`.

The previous default, `global53_abc`, omitted `C*14:03` because the older
MHCflurry model bundles encoded both C*14 alleles identically and local CTA-MS
support favored `C*14:02`. Its membership remains unchanged; use
`--panel global53_abc` to reproduce that panel.

The MHCflurry 2.3.0 presentation bundle changes the trained representation from
37 to 39 residues and distinguishes these alleles at pseudosequence position
5 (R for `C*14:02`, H for `C*14:03`). The 34-residue representation is the
NetMHCpan-derived reference, not the old trained MHCflurry encoding. The
standalone sequence download already contained 39 residues in download release
2.2.0, although that release's trained models still used 37. Inspect the
sequences bundled with the actual model, not just its software version.

The 2.3.0 presentation archive's `affinity_predictor_train_data.csv.bz2`
contains 171 `C*14:03` rows and unique peptides, all qualitative MS observations
(151 Keskin; 20 PMID 31844290). Thus this allele has direct training evidence;
training inclusion and differing percentile ranks alone do not establish
independently validated biological specificity. The panel restores an existing
reference allotype whose encoding-based exclusion no longer applies.

All 54 default alleles have numeric affinity calibration in the audited 2.3.0
bundle. Tsarina selection uses **presentation** percentiles; it does not require
affinity percentiles for scoring. `HLA-C*15:05`, historically omitted because it
lacked affinity calibration in older bundles, can still be supplied explicitly.
The integration tests verify finite scores without freezing an upstream
allele's calibration availability or assuming permanent sequence equivalence.

Model references:

- [MHCflurry pseudosequence definitions](https://github.com/openvax/mhcflurry/blob/2.3.0/mhcflurry/pseudosequences.py)
- [MHCflurry 2.3.0 release](https://github.com/openvax/mhcflurry/releases/tag/2.3.0)
- [Audited presentation model archive](https://github.com/openvax/mhcflurry/releases/download/2.3.0/models_class1_presentation.20260928.tar.bz2)

Panel references:

- [IEDB allele frequencies and reference sets](https://help.iedb.org/hc/en-us/articles/114094151851-HLA-allele-frequencies-and-reference-sets-with-maximal-population-coverage)
- [TepiTool allele-selection description](https://pmc.ncbi.nlm.nih.gov/articles/PMC4981331/)
- [IEDB/Paul common HLA-A/B thresholds](https://help.iedb.org/hc/en-us/articles/114094151811-Selecting-thresholds-cut-offs-for-MHC-class-I-and-II-binding-predictions)
- [Sarkizova et al. HLA-C peptidome study](https://doi.org/10.1038/s41587-019-0322-9)

## Population coverage

Regional allele frequencies from seven geographic regions support the coverage
estimate. Subpopulation proxy rows remain separate from published global
averages and use the same 0–1 allele-frequency scale.

For each allele, Tsarina uses a regional weighted value when a numeric proxy is
available and falls back to the published global average otherwise. For each
CTA:

1. sum covered allele frequencies within each HLA locus;
2. convert each locus frequency \(f\) to carrier probability
   \(1 - (1 - f)^2\); and
3. combine carrier probabilities across loci.

All default `global54_abc` alleles have source, proxy, resolution, and nonzero
frequency provenance. The result is a panel-design estimate, not a
clinical-grade population-genetics analysis.

## Advanced selection controls

### Source selection

- `--ctas` supplies an explicit target list.
- `--cta-count` changes the maximum automatic non-empty target count.
- `--min-restriction-confidence` and `--restriction-levels` refine oncoref
  restriction axes.
- `--selection-allowlist` replaces the automatic safety allowlist.
- `--no-vital-tissue-filter` or `--vital-tissue-max-ntpm` changes the vital
  tissue gate.
- `--allow-non-magea4-mage-family` permits other MAGE-family candidates during
  automatic selection.

### Peptides and grouping

- `--lengths` changes enumerated peptide lengths.
- `--no-require-cta-exclusive` permits peptide matches in non-CTA proteins.
- `--peptides-per-cell` changes the maximum retained per CTA × HLA cell.
- `--no-group-identical-cta-peptide-sets` keeps peptide-identical sources
  separate.
- `--no-group-identical-cta-pmhcs` keeps sources with duplicate final panels
  separate.
- `--proteoform-labels {symbol,members}` chooses how identical-protein groups
  are labeled in the output (default `symbol`; see
  [Enumerate and group peptides](#3-enumerate-and-group-peptides)).

### Prediction

- `--predictor` selects a supported presentation backend.
- `netmhcpan` and `netmhcpan_el` both use mhctools' version-detecting
  NetMHCpan adapter. Presentation columns contain EL scores/ranks; affinity
  columns contain BA nM values/ranks. The EL selector remains an alias for
  compatibility and does not require a separate `NetMHCpanEL` Python class.
- `--netmhcpan-affinity` adds NetMHCpan BA affinity nM and percentile
  annotations in a second scoring pass. It requires the external NetMHCpan
  backend.

Run `tsarina panel --help` for the complete current option list.

The separate Python API `build_panel_matrix(metric="peptide_count")` counts
unique peptides with a presentation percentile of at most 1.0 per source
and allele. Set `max_presentation_percentile` to change this inclusive cutoff
(0–100). Missing or invalid ranks never count. This parameter does not change
the definitions of `best_percentile`, `ms_peptide_count`, or `has_peptide`,
or the evidence-tier cutoffs used by `tsarina panel`.
