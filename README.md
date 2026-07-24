# tsarina

[![Tests](https://github.com/pirl-unc/tsarina/actions/workflows/tests.yml/badge.svg)](https://github.com/pirl-unc/tsarina/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/tsarina.svg)](https://pypi.org/project/tsarina/)

Tsarina selects shared cancer immunotherapy targets for an individual patient
or a population HLA panel.

It starts with reusable cancer-testis antigen (CTA), oncogenic-virus, and
recurrent-mutation targets. It then combines tumor context, public
immunopeptidomics evidence, healthy-tissue safety evidence, and predicted HLA
presentation to produce ranked peptide-MHC (pMHC) candidates.

## Choose a workflow

| Goal | Entry point | Result |
|---|---|---|
| Prioritize targets for one patient | `tsarina personalize` or `personalized_targets()` | Ranked CTA, viral, and mutant pMHCs for that patient's HLA type and tumor |
| Design an off-the-shelf CTA panel | `tsarina panel` | CTA × HLA matrix with evidence tiers and population-coverage estimates |
| Inspect public peptide observations | `tsarina hits` | Cancer, healthy-tissue, and restriction evidence for specified peptides |

Start with the [documentation guide](docs/index.md) for inputs, data setup, and
the workflow-specific guides.

## Install

```bash
pip install tsarina
```

Install the optional peptide-generation and partitioning dependencies for full
functionality:

```bash
pip install "tsarina[all]"
```

## Quick start

Prioritize targets for a patient:

```bash
tsarina personalize \
  --hla 'HLA-A*02:01,HLA-A*24:02,HLA-B*07:02' \
  --cta 'MAGEA4=142.5,PRAME=87.3' \
  --mutations "KRAS G12D" \
  --viruses hpv16 \
  --output patient-targets.csv
```

Build the default global CTA × HLA panel:

```bash
tsarina panel --format long --output panel.csv
```

Both workflows can use registered IEDB or CEDAR ligand exports as public
immunopeptidomics evidence. See [Data and evidence](docs/data-and-evidence.md)
for setup and interpretation.

## Target categories

| Category | Candidate source | Tumor-specific context |
|---|---|---|
| CTA | The canonical [oncoref](https://github.com/pirl-unc/oncoref) CTA set | Tumor RNA expression and reproductive-tissue restriction |
| Viral | Proteomes from nine oncogenic viruses | Virus detected in the tumor |
| Mutant | Nineteen recurrent hotspots across seven driver genes | Matching mutation detected in the tumor |

Oncoref is the single authority for CTA membership, aliases, HPA restriction
calls, and proteoform groups. Tsarina adds downstream evidence and scoring
without maintaining a second CTA definition library. The exact boundary is
documented in [CTA ownership and downstream evidence](docs/curation.md).

## How ranking works

Tsarina applies the same high-level sequence across workflows:

1. choose candidates from the relevant shared-target sets;
2. enforce tumor context and exclude non-target human peptide matches;
3. annotate public cancer and healthy-tissue MS observations;
4. predict presentation by the requested HLA alleles; and
5. rank pMHCs with explicit evidence and safety provenance.

Public MS evidence strengthens a candidate but does not redefine CTA
membership. Healthy, direct-ex-vivo observations outside reproductive tissues
and thymus are treated as safety evidence.

## Documentation

- [Documentation guide](docs/index.md) — workflow selection and shared concepts
- [Personalized target selection](docs/personalized-targets.md) — patient
  inputs, prioritization, and output schema
- [CTA panel design](docs/panel-design.md) — automatic selection, evidence
  tiers, HLA panels, and coverage
- [Data and evidence](docs/data-and-evidence.md) — data registry, observation
  classification, scoring, and naming
- [CTA ownership and downstream evidence](docs/curation.md) — oncoref/Tsarina
  responsibilities and maintenance

## Development

```bash
./develop.sh
./format.sh
./lint.sh
./test.sh
```
