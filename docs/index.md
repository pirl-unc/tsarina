# Tsarina documentation

Tsarina ranks peptide-MHC candidates from cancer-testis antigens (CTAs),
oncogenic viruses and recurrent mutations. It supports patient target
selection, population CTA panels and shared CTA vaccine constructs.

## Workflows

### [Personalized target selection](personalized-targets.md)

Patient HLA type, CTA expression, hotspot mutations and viral status determine
the candidates. Results include ranked pMHCs with source abundance, public MS
support, healthy-tissue safety flags and predicted presentation.

### [CTA panel design](panel-design.md)

Automatic selection applies CTA safety filters and cancer-prevalence ranking;
an explicit CTA list is also supported. Results include a CTA × HLA matrix,
evidence tiers and population-coverage estimates.

### [CTA vaccine design](vaccine-design.md)

Rank proteoforms by mortality-weighted p95 prevalence, retain CTA-specific
MS-supported sequence, and audit a single DNA/RNA construct.

See: [Vaccine Atlas](vaccine-results/index.html) — strict/loose antigen
comparisons, source-resolved tissue maps, HLA coverage and cumulative MS evidence.

### [Data and evidence](data-and-evidence.md)

Register IEDB and CEDAR exports, query peptide observations, and interpret
their tissue context, HLA restrictions and prediction scores.

## Target selection

Target selection follows five stages:

1. **Define candidates.** CTA definitions come from oncoref; viral proteins and
   recurrent mutation hotspots come from Tsarina's target modules.
2. **Apply biological context.** Patient workflows retain expressed CTAs and
   detected mutations or viruses. Panel workflows rank CTAs by population-level
   cancer prevalence.
3. **Enforce specificity and safety.** Candidate peptides must be exclusive to
   their intended source rules. Healthy-tissue MS observations are kept
   separate from expected reproductive-tissue and thymus observations.
4. **Add presentation evidence.** Public immunopeptidomics observations and HLA
   presentation predictions establish evidence tiers.
5. **Rank with provenance.** Results retain target identity, observation
   evidence, prediction scores, and the reason for their rank.

[OncoRef](https://github.com/pirl-unc/oncoref) owns CTA membership, aliases,
HPA restriction calls and proteoform groups. Tsarina adds downstream evidence
for selection and ranking.

See: [CTA ownership and downstream evidence](curation.md).

## Install and prepare data

Install the package:

```bash
pip install "tsarina[all]"
```

See which external datasets Tsarina can use:

```bash
tsarina data available
```

Viral proteomes can be fetched automatically. IEDB and CEDAR ligand exports
must be downloaded under their source terms and registered locally:

```bash
tsarina data fetch hpv16
tsarina data register iedb /data/mhc_ligand_full.csv
tsarina data register cedar /data/cedar-mhc-ligand-full.csv
tsarina data list
```

See: [Data and evidence](data-and-evidence.md) — storage, observation
classification and the evidence model.

## Target types

### Cancer-testis antigens

CTAs are normally restricted to reproductive tissues but can be reactivated in
tumors. Oncoref supplies the candidate universe and HPA-derived restriction
axes. Tsarina adds tumor prevalence, public-MS safety evidence, peptide
enumeration, and HLA scoring after membership is known.

### Oncogenic viruses

Viral proteins are foreign to the human proteome. Clinical helpers default to
human-exclusive viral peptides, removing viral k-mers that also occur in human
proteins.

### Recurrent mutations

Hotspot driver mutations produce shared mutant peptides in patients carrying
the same alteration. Tsarina enumerates mutation-spanning peptides rather than
discovering private passenger mutations from whole-exome sequencing.

## Guides

| Guide | Use it for |
|---|---|
| [Personalized target selection](personalized-targets.md) | Patient inputs, CLI and Python use, output columns, and ranking |
| [CTA panel design](panel-design.md) | Automatic CTA selection, safety gates, evidence tiers, HLA panels, and coverage |
| [Data and evidence](data-and-evidence.md) | Dataset setup, source classification, positive/negative evidence, scoring, and naming |
| [CTA ownership and downstream evidence](curation.md) | The oncoref/Tsarina boundary, API behavior, and definition updates |

For the current command-line surface, use:

```bash
tsarina --help
tsarina personalize --help
tsarina panel --help
tsarina hits --help
tsarina data --help
```
