# hitlist

Shared cancer immunotherapy targets: CTAs, viral oncoproteins, and recurrent mutant peptides.

Three categories of recurrently targetable MHC ligands that are shared across patients -- the off-the-shelf targets for cancer immunotherapy.

| Category | What makes it cancer-specific | Modules |
|---|---|---|
| **CTA** | Restricted to reproductive tissue, aberrantly expressed in tumors | `gene_sets`, `evidence`, `partition`, `peptides` |
| **Viral** | Foreign proteins from oncogenic viruses, not in normal human tissue | `viral` |
| **Mutant** | Recurrent somatic hotspot mutations shared across many patients | `mutations` |

Plus shared infrastructure: HLA allele panels, IEDB/CEDAR mass spec scanning, MHCflurry scoring, tissue definitions, data management.

## Install

```bash
pip install hitlist

# With partition, peptide, and mutation support (requires pyensembl):
pip install hitlist[all]
```

## Quick start

```python
from hitlist import CTA_gene_names, CTA_evidence

# Default: expressed, reproductive-restricted CTAs
genes = CTA_gene_names()    # set of ~257 gene symbols

# Full evidence table with HPA tissue restriction columns
df = CTA_evidence()
```

## CTA targets

Cancer-testis antigens from five source databases, filtered using Human Protein Atlas tissue expression data.

```python
from hitlist import (
    CTA_gene_names,                # expressed + filtered CTAs (recommended default)
    CTA_gene_ids,                  # same, as Ensembl gene IDs
    CTA_never_expressed_gene_names,# filter-passing but no HPA expression
    CTA_filtered_gene_names,       # all filter-passing (= expressed + never_expressed)
    CTA_excluded_gene_names,       # CTAs that FAIL filter (somatic expression)
    CTA_unfiltered_gene_names,     # full CTA universe (all source databases)
    CTA_evidence,                  # full DataFrame with all evidence columns
)
```

### Set relationships

```
CTA_unfiltered (358 genes)
+-- CTA_filtered (278 genes)
|   +-- CTA_gene_names (257 genes) <-- recommended default
|   +-- CTA_never_expressed (21 genes)
+-- CTA_excluded (80 genes)
```

### Gene partitioning

For pMHC analysis, partition all protein-coding genes into CTA / non-CTA. Requires `pyensembl`:

```python
from hitlist.partition import CTA_partition_gene_names

p = CTA_partition_gene_names()
assert "MAGEA4" in p.cta
assert "TP53" in p.non_cta
```

### CTA peptides

```python
from hitlist.peptides import cta_peptides, cta_exclusive_peptides, build_pmhc_table

exclusive = cta_exclusive_peptides()  # not found in any non-CTA protein
pmhc = build_pmhc_table(exclusive, alleles=get_panel("iedb27_ab"))
```

### Source databases

| Source | Genes | Reference |
|---|---|---|
| [CTpedia](http://www.cta.lncc.br/) | 167 | [Almeida et al. 2009](https://doi.org/10.1093/nar/gkn673), *NAR* |
| [CTexploreR/CTdata](https://www.bioconductor.org/packages/release/bioc/html/CTexploreR.html) | 62 new | [Loriot et al. 2025](https://doi.org/10.1371/journal.pgen.1011734), *PLOS Genetics* |
| Protein-level CT genes (136 total, 46 overlap) | 89 new | [da Silva et al. 2017](https://doi.org/10.18632/oncotarget.21715), *Oncotarget* |
| EWSR1-FLI1 CT gene binding sites | 12 | [Gallegos et al. 2019](https://doi.org/10.1128/MCB.00138-19), *Mol Cell Biol* |
| Meiosis, piRNA, spermatogenesis genes | 28 | Multiple sources (see [docs](docs/curation.md)) |

### Curation figures

![CTA Source Venn Diagram](docs/cta-source-venn.png)
![CTA Filter Funnel](docs/cta-filter-funnel.png)

See [full curation documentation](docs/curation.md) for filter logic, deflated RNA fractions, and protein reliability thresholds.

## Viral targets

Peptides from 9 oncogenic virus proteomes, with human-exclusivity filtering and IEDB overlap:

```python
from hitlist.viral import (
    ONCOGENIC_VIRUSES,
    viral_peptides,
    human_exclusive_viral_peptides,
    cancer_specific_viral_peptides,
    viral_iedb_overlap,
)

# All 8-11mer peptides from HPV-16 (auto-fetches from UniProt)
peps = viral_peptides("hpv16")

# Exclude peptides found in ANY human protein
exclusive = human_exclusive_viral_peptides("hpv16")

# Exclude only non-CTA human proteins (keeps CTA overlaps, adds in_cta_protein column)
cancer_specific = cancer_specific_viral_peptides("hpv16")

# Cross-reference with IEDB mass spec data
hits = viral_iedb_overlap("hpv16", iedb_path="mhc_ligand_full.csv")
```

| Virus | Cancers | Key oncoproteins |
|---|---|---|
| HPV-16, HPV-18 | Cervical, oropharyngeal | E6, E7 |
| EBV/HHV-4 | Burkitt lymphoma, NPC, Hodgkin | LMP1, EBNA1 |
| HTLV-1 | Adult T-cell leukemia/lymphoma | Tax, HBZ |
| HBV | Hepatocellular carcinoma | HBx |
| HCV | HCC, B-cell lymphoma | Core, NS3, NS5A |
| KSHV/HHV-8 | Kaposi sarcoma | vFLIP, LANA |
| MCPyV | Merkel cell carcinoma | LT, sT |
| HIV-1 | Kaposi sarcoma, lymphoma | Tat, Nef |

## Mutant targets

Peptides spanning recurrent somatic hotspot mutations -- shared neoantigens present across many patients:

```python
from hitlist.mutations import HOTSPOT_MUTATIONS, mutant_peptides, mutant_iedb_overlap

# Generate all mutation-spanning 8-11mer peptides
df = mutant_peptides()

# Cross-reference with IEDB
hits = mutant_iedb_overlap(iedb_path="mhc_ligand_full.csv")
```

19 hotspot mutations across 8 genes:

| Gene | Mutations | Cancer types |
|---|---|---|
| KRAS | G12C, G12D, G12V, G12R, G13D | Pancreatic, colorectal, NSCLC |
| BRAF | V600E, V600K | Melanoma, colorectal, thyroid |
| TP53 | R175H, R248W, R273H, G245S, R249S | Pan-cancer |
| PIK3CA | H1047R, E545K | Breast, endometrial |
| IDH1 | R132H | Glioma, AML |
| NRAS | Q61R, Q61K | Melanoma |
| EGFR | L858R, T790M | NSCLC |

## HLA allele panels

Population-spanning MHC class I allele panels, organized as nested tiers:

```python
from hitlist.alleles import get_panel, panel_names, PANEL_SOURCE_CATEGORIES

alleles = get_panel("iedb36_abc")  # 36 alleles (A/B/C)
```

| Panel | Alleles | Coverage |
|---|---|---|
| `iedb27_ab` | 27 | Global baseline (IEDB/TepiTool) |
| `iedb36_abc` | 36 | + HLA-C |
| `global44_abc` | 44 | + East Asia, South Asia, Sub-Saharan Africa |
| `global48_abc` | 48 | + Latin America, MENA |
| `global51_abc_ssa` | 51 | + additional Sub-Saharan Africa |

Regional allele frequencies from 7 geographic regions:

```python
from hitlist.regions import region_allele_frequencies, REGION_POPULATIONS
freq_df = region_allele_frequencies()
```

## IEDB/CEDAR mass spec scanning

Scan IEDB and CEDAR MHC ligand exports with MHC class filtering and disease/tissue context classification:

```python
from hitlist.iedb import scan_public_ms, profile_dataset

# Scan for specific peptides (MHC class I, with source context)
hits = scan_public_ms(
    peptides=my_peptides,
    iedb_path="mhc_ligand_full.csv",
    mhc_class="I",
    classify_source=True,  # adds src_cancer, src_healthy, src_cell_line, ...
)

# Profile entire dataset for summary statistics
df = profile_dataset(iedb_path="mhc_ligand_full.csv")
```

## Data management

```bash
hitlist data available                          # show all known datasets
hitlist data fetch hpv16                        # auto-download viral proteome
hitlist data register iedb /path/to/file.csv    # register manual download
hitlist data list                               # show installed datasets
hitlist data path iedb                          # resolve to file path
```

## Development

```bash
./develop.sh    # install in dev mode
./format.sh     # ruff format
./lint.sh       # ruff check + format check
./test.sh       # pytest with coverage
```
