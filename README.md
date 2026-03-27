# ctabase

Curated cancer-testis antigen (CTA) gene sets, HLA allele panels, and pMHC utilities for cancer immunotherapy research.

## Install

```bash
pip install ctabase

# With partition and peptide support (requires pyensembl):
pip install ctabase[all]
```

## Quick start

```python
from ctabase import CTA_gene_names, CTA_evidence

# Default: expressed, reproductive-restricted CTAs
genes = CTA_gene_names()    # set of ~257 gene symbols
print(f"{len(genes)} CTAs: {sorted(genes)[:5]}...")

# Full evidence table with HPA tissue restriction columns
df = CTA_evidence()
print(df[["Symbol", "source_databases", "rna_deflated_reproductive_frac", "filtered"]].head())
```

## Gene sets

ctabase provides the CTA gene universe from five source databases, systematically filtered using Human Protein Atlas tissue expression data to retain only genes with reproductive-restricted expression.

```python
from ctabase import (
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

### Evidence table

```python
df = CTA_evidence()

# Strict CTAs from CTpedia with Enhanced protein evidence
strict = df[
    (df["filtered"] == True) &
    (df["source_databases"].str.contains("CTpedia")) &
    (df["protein_reliability"] == "Enhanced") &
    (~df["never_expressed"])
]

# Genes with tumor mass spec evidence
tumor_protein = df[df["source_databases"].str.contains("daSilva2017_protein", na=False)]
```

## Source databases

| Source | Genes | Reference |
|---|---|---|
| [CTpedia](http://www.cta.lncc.br/) | 167 | [Almeida et al. 2009](https://doi.org/10.1093/nar/gkn673), *NAR* |
| [CTexploreR/CTdata](https://www.bioconductor.org/packages/release/bioc/html/CTexploreR.html) | 62 new | [Loriot et al. 2025](https://doi.org/10.1371/journal.pgen.1011734), *PLOS Genetics* |
| Protein-level CT genes (136 total, 46 overlap) | 89 new | [da Silva et al. 2017](https://doi.org/10.18632/oncotarget.21715), *Oncotarget* |
| EWSR1-FLI1 CT gene binding sites | 12 | [Gallegos et al. 2019](https://doi.org/10.1128/MCB.00138-19), *Mol Cell Biol* |
| Meiosis, piRNA, spermatogenesis genes | 28 | Multiple sources (see [docs](docs/curation.md)) |

## Curation pipeline

### Step 1: Collect

Union of protein-coding CT genes from five source databases (358 genes). Only protein-coding genes (Ensembl biotype) are included. Gene symbols are updated to current HGNC nomenclature.

### Step 2: Annotate with HPA

Each gene is scored against [Human Protein Atlas](https://www.proteinatlas.org/) v23:

- **RNA**: nTPM across 50 normal tissues. A *deflated reproductive fraction* suppresses sub-1 nTPM noise: `(1 + sum_repro max(0, nTPM-1)) / (1 + sum_all max(0, nTPM-1))`. Thymus excluded from denominator (AIRE-mediated expression is expected for CTAs).
- **Protein**: IHC staining across 63 tissues with antibody reliability scores (Enhanced > Supported > Approved > Uncertain).

### Step 3: Filter

Tiered deflated RNA thresholds that scale with protein data confidence:

| Protein evidence | Required deflated RNA fraction |
|---|---|
| Enhanced + reproductive only | >= 80% |
| Supported + reproductive only | >= 90% |
| Approved + reproductive only | >= 95% |
| Uncertain or no protein data | >= 99% |

Genes with protein detected in non-reproductive tissues always fail. See [curation docs](docs/curation.md) for full details.

### Figures

#### Source overlap
![CTA Source Venn Diagram](docs/cta-source-venn.png)

#### Filter funnel
![CTA Filter Funnel](docs/cta-filter-funnel.png)

#### Filter outcome by source
![CTA Filter Outcome](docs/cta-filter-outcome.png)

#### Deflated fraction distribution
![Deflated Fraction Distribution](docs/cta-deflated-frac-dist.png)

## Tissue definitions

Three tiers of reproductive tissue sets for CTA restriction analysis:

```python
from ctabase.tissues import (
    CORE_REPRODUCTIVE_TISSUES,      # {testis, ovary, placenta}
    EXTENDED_REPRODUCTIVE_TISSUES,  # + cervix, endometrium, prostate, ...
    PERMISSIVE_REPRODUCTIVE_TISSUES,# + breast
    HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS,
    is_tissue_restricted,
    adaptive_rna_threshold,
)

# Check if a gene's detected tissues are restricted to reproductive
is_tissue_restricted({"testis", "thymus"})        # True (thymus excluded)
is_tissue_restricted({"testis", "liver"})          # False

# Adaptive RNA thresholds scale with protein evidence quality
adaptive_rna_threshold("Enhanced")   # 0.80
adaptive_rna_threshold("Uncertain")  # 0.99
```

## Gene partitioning

For pMHC discovery, every protein-coding gene needs to go into exactly one bucket. Requires `pyensembl`:

```python
from ctabase.partition import CTA_partition_gene_names, CTA_partition_dataframes

# As sets of gene symbols
p = CTA_partition_gene_names()
assert "MAGEA4" in p.cta
assert "TP53" in p.non_cta
assert len(p.cta & p.non_cta) == 0  # no overlap

# As DataFrames (CTAs include full evidence columns)
frames = CTA_partition_dataframes()
print(frames.cta[["Symbol", "rna_deflated_reproductive_frac"]].head())
```

## HLA allele panels

Population-spanning MHC class I allele panels for vaccine design, organized as nested tiers:

```python
from ctabase.alleles import get_panel, panel_names, PANEL_SOURCE_CATEGORIES

# Available panels (each is a superset of the previous)
print(panel_names())
# ['iedb27_ab', 'iedb36_abc', 'global44_abc', 'global48_abc', 'global51_abc_ssa']

# Get a specific panel
alleles = get_panel("iedb36_abc")  # 36 alleles (A/B/C)
print(f"{len(alleles)} alleles: {alleles[:3]}...")

# See where each allele comes from
print(PANEL_SOURCE_CATEGORIES["HLA-B*46:01"])
# "East Asia / Southeast Asia add-on"
```

| Panel | Alleles | Coverage |
|---|---|---|
| `iedb27_ab` | 27 | Global baseline (IEDB/TepiTool published panel) |
| `iedb36_abc` | 36 | + HLA-C coverage |
| `global44_abc` | 44 | + East Asia, South Asia, Sub-Saharan Africa |
| `global48_abc` | 48 | + Latin America, MENA |
| `global51_abc_ssa` | 51 | + additional Sub-Saharan Africa |

### Regional allele frequencies

Allele panels can be augmented or replaced with custom allele sets. The `regions` module provides per-region allele frequency data from published population studies:

```python
from ctabase.regions import region_allele_frequencies, REGION_POPULATIONS

# DataFrame with allele frequencies by region, proxy population, and locus
freq_df = region_allele_frequencies()
print(freq_df[freq_df["region"] == "East Asia"][["allele", "frequency", "proxy"]].head())

# Regional population estimates (millions, 2024)
print(REGION_POPULATIONS)
# {'Europe': 743.9, 'MENA / Arab': 599.9, 'East Asia': 1603.2, ...}
```

## IEDB/CEDAR mass spec scanning

Scan downloaded IEDB and CEDAR MHC ligand exports for validated peptide observations:

```python
from ctabase.iedb import scan_public_ms, peptide_ms_support

# Basic scan
hits = scan_public_ms(
    peptides={"SLYNTVATL", "GILGFVFTL"},
    iedb_path="mhc_ligand_full.csv",
    cedar_path="cedar-mhc-ligand-full.csv",
)

# With MHC class I filtering and disease/tissue context classification
hits = scan_public_ms(
    peptides=my_peptides,
    iedb_path="mhc_ligand_full.csv",
    mhc_class="I",             # MHC class I only
    classify_source=True,       # adds context columns
)
# Context columns: src_cancer, src_healthy, src_cell_line, src_ex_vivo, src_ebv_lcl
# Additional: disease, source_tissue, cell_name, process_type, culture_condition

# Quick lookup: peptide -> set of MHC restrictions
support = peptide_ms_support(
    peptides=my_peptides,
    iedb_path="mhc_ligand_full.csv",
    mhc_class="I",
)
```

## MHCflurry scoring

Score peptide-MHC binding and presentation (requires `mhcflurry`):

```python
from ctabase.scoring import score_presentation, AFFINITY_THRESHOLDS_NM
from ctabase.alleles import get_panel

scores = score_presentation(
    peptides=["SLYNTVATL", "GILGFVFTL"],
    alleles=get_panel("iedb27_ab"),
)
# Columns: peptide, allele, presentation_score, presentation_percentile, affinity_nm

# Filter to strong binders
strong = scores[scores["presentation_percentile"] <= 1.0]
```

## Peptide generation

Generate peptides from CTA protein sequences and build pMHC tables. Requires `pyensembl`:

```python
from ctabase.peptides import cta_peptides, cta_exclusive_peptides, build_pmhc_table

# All 8-11mer peptides from expressed CTA proteins
peptides = cta_peptides()
print(f"{len(peptides)} total peptide occurrences")

# Only peptides exclusive to CTA proteins (not found in any non-CTA protein)
exclusive = cta_exclusive_peptides()
print(f"{len(exclusive)} CTA-exclusive peptide occurrences")

# Cross with an allele panel to build a pMHC table
pmhc = build_pmhc_table(exclusive, alleles=get_panel("iedb27_ab"))
print(f"{len(pmhc)} pMHCs")
```

## Development

```bash
# Install in development mode
./develop.sh

# Format, lint, test
./format.sh
./lint.sh
./test.sh
```
