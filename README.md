# tsarina

[![Tests](https://github.com/pirl-unc/tsarina/actions/workflows/tests.yml/badge.svg)](https://github.com/pirl-unc/tsarina/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/tsarina.svg)](https://pypi.org/project/tsarina/)

Tsarina selects shared cancer immunotherapy targets for an individual patient
or a population HLA panel.

It starts with reusable cancer-testis antigen (CTA), curated cancer fusion,
oncogenic-virus, and recurrent-mutation targets. It then combines tumor
context, biological sex eligibility, public immunopeptidomics evidence,
healthy-tissue safety evidence, and predicted MHC presentation to produce
ranked peptide-MHC (pMHC) candidates.

## Workflows

| Goal | Entry point | Result |
|---|---|---|
| Check setup & browse built-in catalogs | `tsarina status` / `tsarina catalog` | Read-only environment check and CTA, HLA panel, mutation, fusion, and virus catalogs |
| Prioritize targets for one patient | `tsarina personalize` or `personalized_targets()` | Ranked CTA, fusion, viral, and mutant pMHCs for a human, mouse, or dog host (optional `--sex`) |
| Design an off-the-shelf CTA panel | `tsarina panel` or `spanning_pmhc_set()` | CTA × HLA matrix with evidence tiers, MS safety flags, and population coverage |
| Assemble a shared CTA vaccine | `tsarina vaccine` or `design_vaccine()` | Mortality-ranked proteoforms, sequence funnel, junction audit, and DNA/RNA construct |
| Inspect public peptide observations | `tsarina hits` or `cta_healthy_tissue_ms_hits()` | Cancer, healthy-somatic, sex-reproductive, and restriction evidence for a gene |

See: [documentation](docs/index.md) and [Ecosystem & Python API](docs/ecosystem-and-api.md).

## How Tsarina connects to the stack

| Package | Owns |
|---|---|
| [oncoref](https://github.com/pirl-unc/oncoref) | Gene-level CTA & 4-way sex-restricted selection (evaluated human, mouse, and dog), clinical clearance, 1-mm near-self metrics, proteoform groups, 188 curated cancer fusions, TCGA/HPA/GLOBOCAN burden |
| [hitlist](https://github.com/pirl-unc/hitlist) | Public MS observations (IEDB/CEDAR/studies), cancer vs. healthy vs. sex-reproductive tissue classification, `EvidenceStore`, dataset registry |
| [pyensembl](https://github.com/openvax/pyensembl) / [orfdata](https://github.com/openvax/orfdata) | Reference coding proteomes (Ensembl 112 human, 115 mouse, 98 dog) and optional IARC GDC / RefSeq viral ORFs |
| [mhcgnomes](https://github.com/openvax/mhcgnomes) | HLA, H-2, and DLA allele, serotype, and species parsing |
| [mhcflurry](https://github.com/openvax/mhcflurry) / [mhctools](https://github.com/openvax/mhctools) / [topiary](https://github.com/openvax/topiary) | MHC presentation and affinity prediction (`pepsickle` + `dnachisel` for vaccine cleavage & mRNA design) |

## Install

Requires Python 3.11 or newer, `oncoref>=1.55.0`, `hitlist>=2.15.0`, and `datacache>=1.38.0`.

```bash
pip install tsarina
```

Install the optional peptide-generation and partitioning dependencies for full
functionality:

```bash
pip install "tsarina[all]"
```

## Quick start

Check local readiness and browse built-in catalogs (read-only, downloads nothing):

```bash
tsarina status
tsarina catalog
tsarina catalog ctas --tier core --top 10
tsarina catalog fusions --cancer SARC_EWS
```

Prioritize targets for a patient:

```bash
tsarina personalize \
  --hla 'HLA-A*02:01,HLA-A*24:02,HLA-B*07:02' \
  --cta 'MAGEA4=142.5,PRAME=87.3' \
  --mutations "KRAS G12D" \
  --fusion 'EWSR1-FLI1:SYGQQSSY|NPSYDSVR' \
  --viruses hpv16 \
  --sex female \
  --output patient-targets.csv
```

Build the default global CTA × HLA panel:

```bash
tsarina panel --format long --output panel.csv
```

Or use the Python API directly:

```python
from tsarina import CTA_selection, cancer_fusions, personalized_targets, spanning_pmhc_set

ctas = CTA_selection(category="CTA")                 # human (or species="mouse")
ewing = cancer_fusions(cancer_code="SARC_EWS")
targets = personalized_targets(
    hla_alleles=["HLA-A*02:01", "HLA-B*07:02"],
    cta_expression={"PRAME": 87.3, "MAGEA4": 142.5},
    mutations=["KRAS G12D"],
    patient_sex="female",
)
```

Fetch public immunopeptidomics data with `tsarina data fetch iedb` /
`tsarina data fetch cedar` (or inspect grouped catalogs via
`tsarina data available --grouped`). See [Data and evidence](docs/data-and-evidence.md)
for setup and interpretation.

## Target categories

| Category | Candidate source | Tumor-specific context |
|---|---|---|
| CTA | Current-policy [oncoref](https://github.com/pirl-unc/oncoref) CTA & 4-way sex-restricted sets (evaluated human, mouse, and dog) | Tumor RNA expression, reproductive-tissue restriction, and optional biological sex eligibility |
| Fusion | 188 curated [oncoref](https://github.com/pirl-unc/oncoref) cancer fusions + host-exclusive breakpoint k-mer enumeration | Fusion breakpoint sequence detected in the tumor |
| Viral | 9 built-in UniProt oncogenic viruses plus optional `orfdata` IARC GDC subsets (15/20/37 viruses) & RefSeq viral ORFs | Virus detected in the tumor |
| Mutant | 19 recurrent hotspots across 7 driver genes | Matching mutation detected in the tumor |

Oncoref is the single authority for CTA membership, aliases, HPA restriction
calls, and proteoform groups. Tsarina adds downstream evidence and scoring
without maintaining a second CTA definition library. Default human membership
uses `oncoref.target_selection.results()`; historical HPA annotations remain
available through `CTA_hpa_evidence()`. The exact boundary is
documented in [CTA ownership and downstream evidence](docs/curation.md).

## How ranking works

Target selection follows five stages:

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
- [Ecosystem & Python API](docs/ecosystem-and-api.md) — how Tsarina connects to
  `oncoref`, `hitlist`, `pyensembl`, `orfdata`, `mhcgnomes`, and predictors, plus
  the curated Python API map
- [Personalized target selection](docs/personalized-targets.md) — patient
  inputs, sex-restricted eligibility, fusions, mouse CTAs, and output schema
- [CTA panel design](docs/panel-design.md) — automatic selection, evidence
  tiers, HLA panels, and coverage
- [CTA vaccine design](docs/vaccine-design.md) — mortality-weighted p95 selection,
  current-policy selection, historical strict/loose comparisons, MS-supported stretches, construct search and full audit
- [Canine evidence and DLA design](docs/canine-vaccine.md) — frozen dog inputs,
  independent-dog RNA prevalence, exact MS support and exploratory DLA constructs
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

`develop.sh` uses the active virtualenv (or creates `.venv`) and requires
pip >=23.0. It resolves Tsarina and present sibling checkouts together, checks
the proposed versions against installed consumers before changing packages,
and runs `pip check` after installation. If a checkout is incompatible, update
that checkout or use a separate virtualenv; the refused install leaves the
installed packages unchanged. Missing siblings use their released packages.

Real MHCflurry model integration tests run in a separate CI job. To run them
locally (Python 3.11+):

```bash
pip install -e ".[test-mhcflurry]"
mhcflurry downloads fetch models_class1_presentation
pytest tests/test_mhcflurry_integration.py --run-mhcflurry
```

These tests check the downloaded presentation bundle, including its affinity
calibration, through Tsarina's scoring interface. Enabling them requires the
dependency and models; missing prerequisites fail rather than skip. Use
`./test.sh --run-mhcflurry` to include them in the full suite.
`./deploy.sh` requires this extra and model bundle and always runs these tests.
