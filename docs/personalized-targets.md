# Personalized target selection

The personalized workflow ranks shared peptide-MHC targets for one patient. It
uses the patient's HLA type and tumor findings to choose candidates, removes
known healthy-tissue risks by default, and returns explicit evidence tiers.

Use this workflow when HLA typing is available and at least one of the
following is known:

- tumor CTA expression in TPM;
- a supported recurrent hotspot mutation; or
- the presence of a supported oncogenic virus.

This is target prioritization from curated shared antigens. It is not
whole-exome discovery of private neoantigens.

## Inputs

| Input | Required | Format | Effect |
|---|---|---|---|
| HLA class I alleles | Yes | HLA names such as `HLA-A*02:01` | Defines the presentation-scoring space |
| CTA expression | No | Gene symbol → tumor RNA TPM | Retains expressed, restriction-qualified CTAs |
| Mutations | No | Supported labels such as `KRAS G12D` | Adds mutation-spanning shared neoantigens |
| Viruses | No | Keys such as `hpv16` or `ebv` | Adds peptides from the corresponding viral proteome |
| IEDB/CEDAR data | Recommended | Registered dataset or explicit path | Adds public cancer and healthy-tissue MS evidence |

On the command line, `--hla` and `--cta` accept comma- and/or space-separated
entries, quoted or not (`--hla HLA-A*02:01 HLA-B*07:02` works as well as
`--hla 'HLA-A*02:01,HLA-B*07:02'`). An allele's `*` is optional —
`HLA-A0201` and `A0201` both resolve to `HLA-A*02:01` — which means `--hla`
never strictly needs quoting (`*` is a shell glob character).

A `--cta` entry's `=TPM` half is optional too: `--cta MAGEA4` (no TPM) means
"include this gene regardless of `--min-cta-tpm`" rather than "assume zero
expression." Mix bare and `GENE=TPM` entries freely.

At least one target source—CTA expression, mutations, or viruses—must produce
candidates for the result to be non-empty.

## Quick start

### Command line

```bash
tsarina personalize \
  --hla 'HLA-A*02:01,HLA-A*24:02,HLA-B*07:02,HLA-B*44:02' \
  --cta 'MAGEA4=142.5,PRAME=87.3,CTAG1B=215.0' \
  --mutations 'KRAS G12D,TP53 R175H' \
  --viruses hpv16 \
  --output patient-targets.csv
```

IEDB and CEDAR paths resolve from the data registry by default. For a dry run
without public MS evidence, pass `--skip-ms-evidence`.

### Python

```python
from tsarina import personalized_targets

targets = personalized_targets(
    hla_alleles=[
        "HLA-A*02:01",
        "HLA-A*24:02",
        "HLA-B*07:02",
        "HLA-B*44:02",
    ],
    cta_expression={
        "MAGEA4": 142.5,
        "PRAME": 87.3,
        "CTAG1B": 215.0,
    },
    mutations=["KRAS G12D", "TP53 R175H"],
    viruses=["hpv16"],
)
```

## What the workflow does

The default workflow applies the following stages in order:

1. **Choose tumor-relevant sources.** CTAs must have at least 2.0 TPM and
   `HIGH` or `MODERATE` oncoref restriction confidence. Mutations and viruses
   must match supported identifiers.
2. **Enumerate source peptides.** The default lengths are 8, 9, 10, and 11
   amino acids.
3. **Enforce source specificity.** CTA peptides must not occur in a non-CTA
   human protein. Viral peptides must not occur anywhere in the human proteome.
   Mutation peptides must span the altered residue and differ from wild type.
4. **Attach public MS evidence.** Registered IEDB/CEDAR observations are
   aggregated by peptide.
5. **Apply the healthy-tissue gate.** Peptides observed on healthy,
   non-reproductive tissue are removed by default.
6. **Score HLA presentation.** Each peptide is scored against all supplied
   alleles; the best allele and its score are retained.
7. **Assign tiers and sort.** Results are ordered deterministically by tier, MS
   support, presentation percentile, peptide, and best allele.

See [Data and evidence](data-and-evidence.md) for the exact observation
classification used by the MS and healthy-tissue stages.

## Evidence tiers

Lower presentation percentiles are better. The default MHCflurry-calibrated
tiers are:

| Tier | Label | Requirements |
|---:|---|---|
| 1 | `STRONG` | Percentile ≤ 0.25 and either cancer MS evidence or a viral/mutant source |
| 2 | `MODERATE` | Percentile ≤ 0.50 and either any MS evidence or a viral/mutant source |
| 3 | `CANDIDATE` | Percentile ≤ 1.0 |
| 4 | `WEAK` | Does not meet the above, including unscored rows |

Tier 4 is excluded by default. Use `--keep-weak-tier` or
`drop_weak_tier=False` for diagnostics. The `tier_reason` column records which
condition established each tier.

## Output schema

The result is one row per source peptide, annotated with its best patient HLA
allele.

| Column | Meaning |
|---|---|
| `peptide`, `length` | Peptide sequence and length |
| `category` | `cta`, `cta_flagged`, `viral`, or `mutant` |
| `source` | CTA symbol, virus name, or mutation label |
| `source_detail` | Ensembl gene ID, viral protein accession, or mutation string |
| `source_tpm` | Patient tumor RNA expression for CTA rows |
| `flag_reason` | Why oncoref excludes this gene from the strict CTA set, for `cta_flagged` rows only |
| `ms_hit_count` | Number of aggregated public MS observations |
| `ms_alleles`, `ms_allele_count` | Observed HLA restrictions and their count |
| `ms_in_cancer` | Whether public cancer MS evidence exists |
| `ms_in_healthy_tissue` | Whether disallowed healthy-tissue MS evidence exists |
| `best_allele` | Supplied patient allele with the best predicted presentation |
| `presentation_percentile` | Best presentation percentile; lower is better |
| `presentation_score` | Presentation score from the selected predictor |
| `affinity_nm` | Predicted binding affinity in nM |
| `tier`, `tier_label`, `tier_reason` | Prioritization result and provenance |

## Supported shared targets

### Cancer-testis antigens

CTA definitions come directly from
[oncoref](https://github.com/pirl-unc/oncoref). Tsarina does not maintain a
second CTA list. The patient workflow starts from the canonical oncoref set,
then applies expression, restriction-confidence, peptide-exclusivity, and
downstream evidence gates.

A `--cta` gene oncoref excludes from that strict set but still tracks as a
known clinical target (CTAG2/LAGE-1 is the motivating example: excluded for a
low-level HPA heart RNA signal, kept because it shares the NY-ESO-1 157-165
epitope targeted by the approved TCR-T afami-cel) is not silently dropped. It
appears with `category="cta_flagged"` and the exclusion reason in
`flag_reason`, so a caller who named the gene explicitly sees it and its
caveat. Its peptides are generated directly from the protein sequence and are
not screened for overlap with other proteins the way a strict CTA's are — see
`flag_reason`. A `--cta` gene that is neither a recognized CTA nor a known
clinical target (a typo, or a gene with no CTA evidence at all) is dropped
with a warning naming it.

#### Identical-protein groups

CTAs that translate to a byte-identical protein are reported as one group
by default, under oncoref's preferred symbol for it: CTAG1A + CTAG1B is
reported as `NY-ESO-1`, XAGE1A + XAGE1B as `XAGE1A/B`, SSX2 + SSX2B as
`SSX2/B`, and the SSX4, MAGEA2, MAGEA9, CT45A, CT47A and GAGE12 families
the same way. Both halves come from oncoref: `proteoform_symbol_map` for
membership (the registry the panel workflow uses) and `proteoform_symbol`
for the name, which is a curated alias where one exists and the
prefix-contracted members otherwise. Neither is restated in tsarina.

This matters in both directions. Naming one member (`--cta CTAG1B=215`)
still reports `CTAG1A/CTAG1B`, because the peptides are not unique to the
member you named. Naming both collapses them to a single row instead of
double-counting one finding: the group keeps the highest TPM any member
reported (identical proteins, so an RNA quantifier splits reads between
the loci more or less arbitrarily, and summing would inflate one real
signal) and the union of their Ensembl gene IDs in `source_detail`.

Pass `--no-proteoform-rollup` for one row per gene symbol.

See [CTA ownership and downstream evidence](curation.md) for the definition
boundary.

### Oncogenic viruses

| Virus | Associated cancers | Representative oncoproteins |
|---|---|---|
| HPV-16, HPV-18 | Cervical, oropharyngeal, anal | E6, E7 |
| EBV/HHV-4 | Burkitt lymphoma, nasopharyngeal carcinoma, Hodgkin lymphoma | LMP1, EBNA1, LMP2A |
| HTLV-1 | Adult T-cell leukemia/lymphoma | Tax, HBZ |
| HBV | Hepatocellular carcinoma | HBx |
| HCV | Hepatocellular carcinoma, B-cell lymphoma | Core, NS3, NS5A |
| KSHV/HHV-8 | Kaposi sarcoma, primary effusion lymphoma | vFLIP, vCyclin, LANA |
| MCPyV | Merkel cell carcinoma | Large T, small T |
| HIV-1 | Kaposi sarcoma, non-Hodgkin lymphoma | Tat, Nef |

The clinical helper uses human-exclusive viral peptides. For exploratory use,
`cancer_specific_viral_peptides()` permits CTA overlaps while excluding
non-CTA human overlaps.

### Recurrent mutations

| Gene | Supported hotspots | Common contexts |
|---|---|---|
| KRAS | G12C, G12D, G12V, G12R, G13D | Pancreatic, colorectal, NSCLC |
| BRAF | V600E, V600K | Melanoma, colorectal, thyroid |
| TP53 | R175H, R248W, R273H, G245S, R249S | Pan-cancer |
| PIK3CA | H1047R, E545K | Breast, endometrial |
| IDH1 | R132H | Glioma, AML |
| NRAS | Q61R, Q61K | Melanoma |
| EGFR | L858R, T790M | NSCLC |

Use `HOTSPOT_MUTATIONS` for the executable list rather than parsing this table.

## Advanced controls

- Change CTA inclusion with `--min-cta-tpm` and
  `--min-restriction-confidence`.
- Add an mTEC expression gate with `--mtec-matrix-path` and
  `--mtec-max-tpm`.
- Change peptide lengths with `--lengths`.
- Select another supported predictor with `--predictor`. Tier thresholds are
  calibrated to MHCflurry, so interpret other predictors cautiously.
- Retain healthy-tissue observations only for investigation with
  `--no-enforce-tumor-specificity`.
- Relax viral human exclusivity only for investigation with
  `--no-require-human-exclusive-viral`.
- Disable presentation scoring with `--no-score`.
- Choose output shape with `--format {table,csv,tsv}`. When omitted it is
  inferred: from `--output`'s extension (`.csv`, `.tsv`/`.tab`, `.txt`),
  else CSV for any other `--output` path, else a compact fixed-width
  table when printing to the terminal.
- Report one row per gene symbol instead of per identical-protein group
  with `--no-proteoform-rollup` (see below).
- Suppress the stage-progress messages personalize prints to stderr by
  default (peptide generation, MS evidence lookup, presentation scoring)
  with `--quiet`.

Run `tsarina personalize --help` for the complete current option list.
