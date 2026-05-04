# tsarina

[![Tests](https://github.com/pirl-unc/tsarina/actions/workflows/tests.yml/badge.svg)](https://github.com/pirl-unc/tsarina/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/tsarina.svg)](https://pypi.org/project/tsarina/)

Personalized cancer immunotherapy target selection from curated shared antigen data.

Perseus weaves patient-specific tumor characteristics (mutations, CTA expression, viral infections, HLA type) together with curated public mass spectrometry evidence to produce prioritized lists of targetable peptide-MHC complexes. The name reflects the goal: using shared, public knowledge to personalize cancer immunotherapy -- like Perseus using borrowed divine weapons to slay Medusa.

## Concept

The core insight is that many cancer-targetable peptides are **shared across patients**: cancer-testis antigens are recurrently activated in tumors, oncogenic viruses produce the same foreign proteins in every infected cell, and hotspot driver mutations generate identical mutant peptides across thousands of patients. Unlike private passenger-mutation neoantigens that require individual whole-exome sequencing, these shared targets can be curated once and reused.

Perseus combines curated shared targets with per-patient tumor data to produce a prioritized list of peptide-MHC complexes.

**Shared targets** (curated once, reused across patients):
- **CTA genes** — 358 curated from 5 databases, 257 after HPA tissue expression filtering
- **Viral proteomes** — 9 oncogenic viruses (HPV, EBV, HBV, HCV, HTLV-1, HIV, HHV-8, MCPyV, MCV)
- **Hotspot mutations** — 19 recurrent mutations across 8 driver genes

**Public annotation data** (used to score and filter targets):
- **Mass spec evidence** — IEDB/CEDAR immunopeptidomics observations
- **Tissue expression** — HPA RNA (50 tissues) + IHC protein (63 tissues)
- **HLA allele panels** — population-representative panels (27–51 alleles per region)

**Patient data** (per-individual):
- HLA type (Class I alleles)
- Tumor RNA-seq (CTA expression in TPM)
- Detected mutations (cross-referenced against hotspot list)
- Viral status (HPV, EBV, etc.)

Perseus filters shared targets through the patient's HLA type and tumor profile, then ranks the results into a **prioritized target list** annotated with:
- **Public MS evidence** — number of independent IEDB/CEDAR references, source context (cancer vs. healthy tissue)
- **Source protein abundance** — RNA expression in TPM, estimated protein abundance where HPA data permits
- **Predicted presentation** — MHCflurry presentation percentile, NetMHCpan binding affinity
- **Target category** — CTA, viral, or mutant, with full provenance

## Install

```bash
pip install tsarina

# With full functionality (pyensembl for peptide generation + gene partition):
pip install tsarina[all]
```

## Three target categories

### CTA (cancer-testis antigens)

Proteins normally restricted to reproductive tissues (testis, ovary, placenta) that become aberrantly expressed in tumors. Their tissue restriction means immune responses against them should not damage normal somatic tissues. Thymus expression is expected (AIRE-mediated central tolerance) and excluded from restriction checks.

**358 genes** from 5 source databases, systematically filtered using Human Protein Atlas v23 tissue expression data to **257 expressed CTAs** with predominantly reproductive-restricted expression (some pass the adaptive filter with minor somatic RNA signal).

```python
from tsarina import CTA_gene_names, CTA_gene_ids, CTA_evidence

genes = CTA_gene_names()    # recommended default set (257 expressed CTAs)
df = CTA_evidence()          # full evidence table with HPA columns + 3-axis tiers
```

**Per-modality restriction** classifies each CTA independently by protein (IHC), RNA, and MS evidence, then synthesizes a unified restriction with confidence:

| Modality | Column | Values |
|----------|--------|--------|
| Protein IHC | `protein_restriction` | TESTIS / PLACENTAL / OVARIAN |
| RNA | `rna_restriction` | TESTIS / PLACENTAL / OVARIAN / REPRODUCTIVE |
| RNA quality | `rna_restriction_level` | STRICT / MODERATE / PERMISSIVE |
| MS (runtime) | `ms_restriction` | CANCER_ONLY / EXPECTED_TISSUE / SINGLETON_HEALTHY / RECURRENT_HEALTHY |
| **Synthesized** | `restriction` | TESTIS / PLACENTAL / OVARIAN / REPRODUCTIVE |
| **Confidence** | `restriction_confidence` | HIGH / MODERATE / LOW |

```python
from tsarina import CTA_testis_restricted_gene_names, CTA_by_axes

testis = CTA_testis_restricted_gene_names()  # 229 genes (synthesized TESTIS)
strict_testis = CTA_by_axes(restriction="TESTIS", rna_restriction_level="STRICT")
high_conf = CTA_by_axes(restriction="TESTIS", restriction_confidence="HIGH")
```

| Source | Genes | Reference |
|---|---|---|
| [CTpedia](http://www.cta.lncc.br/) | 167 | [Almeida et al. 2009](https://doi.org/10.1093/nar/gkn673), *NAR* |
| [CTexploreR/CTdata](https://www.bioconductor.org/packages/release/bioc/html/CTexploreR.html) | 62 new | [Loriot et al. 2025](https://doi.org/10.1371/journal.pgen.1011734), *PLOS Genetics* |
| Protein-level CT genes | 89 new | [da Silva et al. 2017](https://doi.org/10.18632/oncotarget.21715), *Oncotarget* |
| EWSR1-FLI1 CT gene binding sites | 12 | [Gallegos et al. 2019](https://doi.org/10.1128/MCB.00138-19), *Mol Cell Biol* |
| Meiosis, piRNA, spermatogenesis genes | 28 | Multiple sources (see [docs](docs/curation.md)) |

#### CTA curation pipeline

**Step 1: Collect.** Take the union of protein-coding CT genes across all 5 source databases → **358 genes**. Duplicates are merged; non-protein-coding genes (e.g., lncRNAs) are excluded.

**Step 2: Annotate for tissue restriction.** The goal is to answer: "is this gene's expression restricted to reproductive tissues?" We use two independent data modalities from Human Protein Atlas v23:

- *RNA expression* (50 tissues): What fraction of total expression comes from reproductive tissues (testis, ovary, placenta)? Raw fractions are misleading because many genes have low-level basal transcription (< 1 nTPM) across dozens of tissues, which inflates the denominator. The **deflated reproductive fraction** fixes this by zeroing out sub-1 nTPM values before computing the ratio, so only tissues with meaningful expression count. Example: CTCFL has testis nTPM = 10.8 but ~40 other tissues at 0.1–0.9 nTPM each. Raw reproductive fraction: 54%. Deflated fraction: 100%, because only testis exceeds 1 nTPM.

- *Protein expression* (63 tissues): Does IHC staining detect protein outside reproductive tissues? Each antibody carries a reliability tier — Enhanced (orthogonal validation), Supported, Approved, or Uncertain — which indicates how much to trust the staining result.

Thymus is excluded from all restriction checks because AIRE drives ectopic expression of tissue-restricted antigens in medullary thymic epithelial cells (mTECs) as part of central tolerance. CTA expression in thymus is expected and does not indicate somatic tissue leakage.

**Step 3: Filter.** Two rules determine whether a gene passes:

1. **Protein exclusion (hard):** If protein is detected in any non-reproductive somatic tissue (excluding thymus), the gene fails — regardless of RNA data.
2. **RNA threshold (tiered):** The required deflated reproductive fraction scales with protein data confidence. When high-quality protein data confirms reproductive restriction, we can tolerate more RNA noise in other tissues. When protein data is absent or unreliable, we demand near-perfect RNA restriction:

   | Protein evidence | Min. deflated RNA reproductive fraction |
   |---|---|
   | Enhanced + reproductive only | >= 80% |
   | Supported + reproductive only | >= 90% |
   | Approved + reproductive only | >= 95% |
   | Uncertain or no protein data | >= 99% |

Result: **257 of 358 genes** pass the filter.

See [full curation documentation](docs/curation.md) for the deflated fraction formula, never-expressed flag, and figures.

### Viral (oncogenic virus proteins)

Foreign proteins from oncogenic viruses -- entirely absent from normal human tissue, making them ideal immunotherapy targets when the virus is present in the tumor.

```python
from tsarina.viral import (
    human_exclusive_viral_peptides,
    viral_peptides,
)

peps = viral_peptides("hpv16")                    # all viral peptides
human_exclusive = human_exclusive_viral_peptides("hpv16")  # default clinical helper
```

`personalize()` and `target_peptides()` use the human-exclusive viral helper by
default, dropping viral k-mers that also occur anywhere in the human proteome.
`cancer_specific_viral_peptides()` is available for exploratory workflows that
allow overlaps with CTA proteins while excluding non-CTA overlaps.

| Virus | Cancers | Key oncoproteins |
|---|---|---|
| HPV-16, HPV-18 | Cervical, oropharyngeal, anal | E6, E7 |
| EBV/HHV-4 | Burkitt lymphoma, NPC, Hodgkin lymphoma | LMP1, EBNA1, LMP2A |
| HTLV-1 | Adult T-cell leukemia/lymphoma | Tax, HBZ |
| HBV | Hepatocellular carcinoma | HBx |
| HCV | HCC, B-cell lymphoma | Core, NS3, NS5A |
| KSHV/HHV-8 | Kaposi sarcoma, primary effusion lymphoma | vFLIP, vCyclin, LANA |
| MCPyV | Merkel cell carcinoma | Large T, small T |
| HIV-1 | Kaposi sarcoma, non-Hodgkin lymphoma | Tat, Nef |

### Mutant (recurrent somatic hotspot mutations)

Shared neoantigens from driver mutations that recur across thousands of patients. Unlike private passenger mutations, these produce the same mutant peptide in every patient carrying the same hotspot mutation.

```python
from tsarina.mutations import HOTSPOT_MUTATIONS, mutant_peptides

df = mutant_peptides()  # all mutation-spanning 8-11mer peptides
```

| Gene | Mutations | Cancer types |
|---|---|---|
| KRAS | G12C, G12D, G12V, G12R, G13D | Pancreatic, colorectal, NSCLC |
| BRAF | V600E, V600K | Melanoma, colorectal, thyroid |
| TP53 | R175H, R248W, R273H, G245S, R249S | Pan-cancer |
| PIK3CA | H1047R, E545K | Breast, endometrial |
| IDH1 | R132H | Glioma, AML (peptidomics-validated vaccine target) |
| NRAS | Q61R, Q61K | Melanoma |
| EGFR | L858R, T790M | NSCLC |

## Positive and negative peptide sets

Perseus constructs both **positive sets** (cancer-specific peptides from the three target categories) and **negative sets** (peptides observed on normal non-reproductive, non-thymic tissues) using the same IEDB/CEDAR scanning infrastructure with consistent tissue classification.

### Tissue source classification

Every IEDB/CEDAR mass spec observation is classified by biological context:

| Category | IEDB criteria | Meaning |
|---|---|---|
| `src_cancer` | Process Type = "Occurrence of cancer" | Peptide detected on tumor MHC |
| `src_healthy` | Process Type = "No immunization", Disease = "healthy" or empty | Peptide detected on normal tissue |
| `src_reproductive` | Source Tissue in {testis, ovary, placenta, ...} | Expected for CTAs |
| `src_thymus` | Source Tissue = thymus | Expected for CTAs (AIRE-mediated) |
| `src_cell_line` | Culture Condition = "Cell Line / Clone" | In vitro, not direct tissue |
| `src_ebv_lcl` | Culture Condition contains "EBV transformed, B-LCL" | EBV-immortalized B cells (special case) |
| `src_ex_vivo` | Culture Condition = "Direct Ex Vivo" | Highest confidence tissue evidence |

**Positive set criteria**: peptide has `src_cancer` evidence AND is exclusive to CTA/viral/mutant source proteins (not found in non-target human proteins).

**Negative set criteria**: peptide has `src_healthy` + `src_ex_vivo` evidence from non-reproductive, non-thymic tissues. These are peptides confirmed to be presented on normal somatic tissue -- targeting them would cause on-target, off-tumor toxicity.

## Patient personalization

The main entry point for clinical use:

```python
from tsarina.personalize import personalize

targets = personalize(
    # Patient HLA type
    hla_alleles=["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02", "HLA-B*44:02",
                 "HLA-C*07:02", "HLA-C*05:01"],

    # CTA expression (gene symbol -> TPM from RNA-seq)
    cta_expression={"MAGEA4": 142.5, "PRAME": 87.3, "CTAG1B": 215.0},

    # Detected mutations (match against hotspot list)
    mutations=["KRAS G12D", "TP53 R175H"],

    # Viral status
    viruses=["hpv16"],

    # Data sources
    iedb_path="mhc_ligand_full.csv",
)
```

Returns a DataFrame with columns:

| Column | Description |
|---|---|
| `peptide` | Peptide sequence |
| `category` | `cta`, `viral`, or `mutant` |
| `source` | Gene name, virus, or mutation label |
| `source_abundance_tpm` | RNA expression in tumor (CTAs only) |
| `ms_hit_count` | Number of IEDB/CEDAR MS observations |
| `ms_alleles` | MHC restrictions observed in public data |
| `ms_in_cancer` | Detected in cancer samples |
| `ms_in_healthy_somatic` | Detected in normal non-reproductive, non-thymic tissue (safety flag) |
| `presentation_percentile` | MHCflurry presentation percentile for best patient allele |
| `best_allele` | Patient HLA allele with best predicted presentation |
| `binding_affinity_nm` | Predicted binding affinity (nM) |

Prioritization is by: (1) public MS evidence strength, (2) source protein abundance, (3) predicted presentation quality, (4) absence of healthy-tissue MS evidence.

## CTA x HLA panel matrices

Build a CTA x HLA pMHC matrix for off-the-shelf panel design:

```bash
tsarina panel
```

Defaults:

- top 25 CTAs ranked by cancer MS peptide count
- automatic safety gates remove CTAs with vital-tissue RNA / unique healthy-MS
  evidence, while allowlisting `PRAME`, `NY-ESO-1`, and `MAGEA4`
- `NY-ESO-1` is treated as one grouped CTA target backed by `CTAG1A` and
  `CTAG1B`
- `global51_abc_calibrated` HLA-A/B/C panel
- 8-11mer CTA-exclusive peptides
- MHCflurry presentation scoring
- MS-evidence-first cell selection
- up to 3 peptides per CTA x HLA cell, ranked by MS source count, then prediction
- readable terminal table plus coverage summary

Use CSV formats for scripts:

```bash
tsarina panel --format long -o panel-long.csv
tsarina panel --format wide -o panel-wide.csv
```

The CLI prints progress for peptide enumeration, public-MS evidence loading,
scoring, evidence-tier construction, and final selection. Interactive terminals
also get a `tqdm` scoring progress bar; use `--no-progress` or
`--no-progress-bars` to suppress it.

Use `--selection-allowlist`, `--no-vital-tissue-filter`, and
`--vital-tissue-max-ntpm` to tune automatic CTA safety filtering. The default
vital RNA cutoff is 2.0 nTPM; public healthy-MS observations in vital tissues
remain exclusionary only when the peptide evidence maps uniquely to that CTA,
unless allowlisted. Explicit `--ctas` accepts aliases such as `NY-ESO-1` and
`MAGE-A4`.

Evidence tiers use configurable presentation percentile cutoffs:

| Evidence tier | Default cutoff | Meaning |
|---|---:|---|
| `monoallelic_ms` | < 2.0 | Peptide observed in mono-allelic MS for that HLA |
| `sample_allele_ms` | < 1.0 | Peptide observed in a multi-allelic sample where this HLA is best among sample alleles |
| `unrestricted_ms` | < 0.5 | Peptide observed by MS with no usable allele assignment |
| `predicted_only` | < 0.1 | No MS support; excluded unless `--include-predicted-only` is passed |

Available HLA panels:

| Panel | Alleles | Coverage |
|---|---|---|
| `iedb27_ab` | 27 | Global baseline (HLA-A/B) |
| `iedb36_abc` | 36 | + HLA-C |
| `global44_abc` | 44 | + East Asia, South Asia, Sub-Saharan Africa |
| `global48_abc` | 48 | + Latin America, MENA |
| `global51_abc_ssa` | 51 | Legacy Global-48 + additional Sub-Saharan Africa |
| `global51_abc_calibrated` | 51 | Default calibrated panel: IEDB A/B backbone, frequent HLA-C allotypes, and IEDB/Paul common-A/B complements |

Regional allele frequency data from 7 geographic regions supports population-weighted coverage calculations.
The default calibrated panel keeps all 27 IEDB/TepiTool class-I A/B reference alleles,
adds all 21 frequent HLA-C allotypes from the Sarkizova HLA-C peptidome coverage set,
and fills the remaining 51-panel slots with the highest-frequency calibrated alleles
missing from the IEDB/Paul 38 common HLA-A/B threshold set
(`B*18:01`, `B*40:02`, `B*46:01`). References: IEDB reference set
<https://help.iedb.org/hc/en-us/articles/114094151851-HLA-allele-frequencies-and-reference-sets-with-maximal-population-coverage>,
TepiTool allele-selection description <https://pmc.ncbi.nlm.nih.gov/articles/PMC4981331/>,
IEDB/Paul 38 common A/B thresholds
<https://help.iedb.org/hc/en-us/articles/114094151811-Selecting-thresholds-cut-offs-for-MHC-class-I-and-II-binding-predictions>,
and Sarkizova et al. <https://doi.org/10.1038/s41587-019-0322-9>.
Audit notes: all 51 default alleles resolve through MHCflurry's
`percent_rank_calibrated_allele` lookup and produce numeric affinity percentile
ranks. `HLA-C*15:05` remains excluded because MHCflurry supports raw affinity and
presentation predictions for it but does not have an affinity percentile-rank
calibration.

## Data management

Tsarina uses the shared hitlist data registry for external datasets:

```bash
# See what data is available
tsarina data available

# Auto-download viral proteomes from UniProt
tsarina data fetch hpv16
tsarina data fetch ebv

# Register manually downloaded IEDB/CEDAR exports
tsarina data register iedb /data/mhc_ligand_full.csv
tsarina data register cedar /data/cedar-mhc-ligand-full.csv

# Inspect what's installed
tsarina data list

# Resolve paths for use in scripts
tsarina data path iedb
```

### Data sources

| Dataset | Source | Size | How to get |
|---|---|---|---|
| IEDB MHC ligand | [iedb.org](https://www.iedb.org/) | ~2 GB | Manual download (terms of use) |
| CEDAR MHC ligand | [cedar.iedb.org](https://cedar.iedb.org/) | ~1 GB | Manual download |
| HPV-16 proteome | [UniProt UP000006729](https://www.uniprot.org/proteomes/UP000006729) | ~3 KB | `tsarina data fetch hpv16` |
| EBV proteome | [UniProt UP000153037](https://www.uniprot.org/proteomes/UP000153037) | ~50 KB | `tsarina data fetch ebv` |
| *(9 viral proteomes total)* | UniProt | varies | `tsarina data fetch <name>` |

Storage location: `~/.hitlist/` (override with `HITLIST_DATA_DIR` env var).
`tsarina data` delegates registry and cache management to hitlist.

IEDB column indices are resolved dynamically from CSV headers, with fallback to known defaults -- robust to IEDB schema changes.

## Tissue definitions

Three tiers of reproductive tissue sets for CTA restriction analysis:

```python
from tsarina.tissues import (
    CORE_REPRODUCTIVE_TISSUES,       # {testis, ovary, placenta}
    EXTENDED_REPRODUCTIVE_TISSUES,   # + cervix, endometrium, prostate, ...
    PERMISSIVE_REPRODUCTIVE_TISSUES, # + breast
    is_tissue_restricted,
    adaptive_rna_threshold,
)
```

## MHCflurry scoring

```python
from tsarina.scoring import score_presentation
from tsarina.alleles import get_panel

scores = score_presentation(
    peptides=["SLYNTVATL", "GILGFVFTL"],
    alleles=get_panel("iedb27_ab"),
)
```

## Target naming convention

Perseus uses a unified naming scheme across all target categories:

| Category | `source` column | `source_detail` column | Example |
|---|---|---|---|
| CTA | Gene symbol | Ensembl gene ID | `MAGEA4` / `ENSG00000147381` |
| Viral | Virus short name | UniProt protein accession | `HPV-16` / `P03126` |
| Mutant | Mutation label | Mutation string | `KRAS G12D` / `G12D` |

## Development

```bash
./develop.sh    # install in dev mode
./format.sh     # ruff format
./lint.sh       # ruff check + format check
./test.sh       # pytest with coverage
```
