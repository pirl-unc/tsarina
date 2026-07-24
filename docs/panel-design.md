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

The default run targets up to 25 non-empty CTAs across `global53_abc`, uses
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
than one Ensembl gene; for example, `CTAG1A/CTAG1B` and its `NY-ESO-1` alias
expand to both genes while retaining one output label.

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
| `global53_abc` | 53 | Default Global-51 extension with CTA-MS-supported alleles |

Use a named panel or provide an explicit list:

```bash
tsarina panel --panel iedb27_ab
tsarina panel --alleles 'HLA-A*02:01,HLA-A*24:02,HLA-B*07:02'
```

### Why Global-53 is the default

`global51_abc` contains:

- the 27 IEDB/TepiTool class-I HLA-A/B reference alleles;
- 21 frequent HLA-C allotypes from the Sarkizova HLA-C peptidome coverage set;
  and
- `B*18:01`, `B*40:02`, and `B*46:01`, the highest-frequency calibrated
  alleles needed to complement the IEDB/Paul common HLA-A/B set.

`global53_abc` adds `A*29:02`, `B*15:02`, and `B*27:05`, which were the top
missing alleles in a public CTA-MS audit. It keeps `C*14:02` but not
`C*14:03`: MHCflurry uses the same pseudosequence and percentile calibration
for both, while local CTA-MS support favored `C*14:02`.

All 53 default alleles resolve through MHCflurry's calibrated
percentile-rank lookup. `HLA-C*15:05` is not in the default panel because
MHCflurry can score its raw affinity and presentation but lacks an affinity
percentile-rank calibration.

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

All default `global53_abc` alleles have source, proxy, resolution, and nonzero
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

### Prediction

- `--predictor` selects a supported presentation backend.
- `--netmhcpan-affinity` adds NetMHCpan BA affinity nM and percentile
  annotations in a second scoring pass. It requires the external NetMHCpan
  backend.

Run `tsarina panel --help` for the complete current option list.
