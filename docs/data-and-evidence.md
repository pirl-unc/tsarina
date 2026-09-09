# Data and evidence

Tsarina keeps three questions separate:

1. **What is the target?** CTA, viral, and mutation definitions establish the
   candidate source.
2. **Where has its peptide been observed?** IEDB/CEDAR observations provide
   cancer, healthy-tissue, and HLA-restriction evidence.
3. **Could this HLA present it?** Prediction fills gaps and compares candidate
   peptide-allele pairs.

This separation is important: downstream MS or prediction evidence can select,
rank, or exclude a candidate, but it cannot add a gene to the canonical CTA
universe.

## Evidence model

### Target definitions

| Target category | Definition authority | Patient or panel context |
|---|---|---|
| CTA | oncoref membership, aliases, HPA restriction, and proteoform groups | Tumor expression or population cancer prevalence |
| Viral | Tsarina's supported oncogenic-virus proteomes | Tumor viral status |
| Mutant | Tsarina's recurrent hotspot registry | Matching tumor mutation |

### Public observations

IEDB and CEDAR ligand exports contribute observed peptides, sample context,
reported HLA restrictions, donor allele sets, tissues, diseases, cell-line
status, publications, and other provenance. Tsarina uses hitlist to normalize
and index these records.

### Prediction

Presentation predictors score candidate peptides against patient alleles or a
population panel. Prediction is used both to rank candidates and, under
stricter thresholds, to assign an allele to peptide-level MS evidence that
lacks a usable exact restriction.

## Observation classification

Each IEDB/CEDAR MS observation is classified by biological source before
aggregation:

| Flag | Source criteria | Interpretation |
|---|---|---|
| `src_cancer` | Process type is occurrence of cancer | Peptide observed in a tumor context |
| `src_healthy` | No immunization and disease is healthy or empty | Peptide observed in a normal context |
| `src_reproductive` | Tissue belongs to the reproductive definition | Expected source context for many CTAs |
| `src_thymus` | Source tissue is thymus | Tracked separately because of AIRE-mediated expression |
| `src_cell_line` | Culture condition is cell line or clone | In-vitro evidence rather than direct tissue |
| `src_ebv_lcl` | EBV-transformed B lymphoblastoid cell line | Special-case transformed-cell evidence |
| `src_ex_vivo` | Culture condition is direct ex vivo | Highest-confidence direct-tissue context |

The flags are not mutually exclusive. For example, an observation can be both
healthy and reproductive.

## Positive and negative evidence

Tsarina applies the same source classification to evidence supporting a target
and evidence warning against it:

- A **positive cancer peptide** has `src_cancer` evidence and is exclusive to
  the applicable CTA, viral, or mutation source rules.
- A **negative safety peptide** has healthy, direct-ex-vivo evidence from
  non-reproductive, non-thymic tissue.

The second category is an on-target, off-tumor warning. Personalized selection
removes such peptides by default. Panel selection applies source- and
peptide-level safety gates while retaining provenance for audit.

Reproductive-tissue and thymus observations are not silently relabeled as
healthy somatic risk. Cell-line and EBV-LCL observations remain distinguishable
from direct-ex-vivo evidence.

## Install and register datasets

List all datasets Tsarina knows about, then inspect what is installed:

```bash
tsarina data available
tsarina data list
```

Fetch supported remotely available assets:

```bash
tsarina data fetch hpv16
tsarina data fetch ebv
```

IEDB and CEDAR exports must be downloaded under their source terms, then
registered:

```bash
tsarina data register iedb /data/mhc_ligand_full.csv
tsarina data register cedar /data/cedar-mhc-ligand-full.csv
```

Inspect metadata or resolve a registered path:

```bash
tsarina data info iedb
tsarina data path iedb
```

Build the normalized observations index once:

```bash
tsarina build observations
```

The index is built automatically on first use if needed. Rebuild it after
replacing an input export:

```bash
tsarina build observations --force
```

Freshness is hitlist's call, not Tsarina's. Every query routes through hitlist,
which compares the stored artifact version and the fingerprints of its curation
files against `observations.parquet`, and the builder contract against
`peptide_mappings.parquet`. Either artifact that no longer matches is rebuilt
automatically, so a curation fix upstream reaches your results on the next
command instead of waiting for a manual `--force`. Validation costs
milliseconds when both artifacts are current, and Tsarina reports on stderr
when a rebuild actually happened.

An index copied in without its IEDB/CEDAR exports cannot be validated — hitlist
needs the sources to fingerprint — so it is used as found.

### Data sources

| Dataset | Upstream source | Acquisition |
|---|---|---|
| IEDB MHC ligand | [IEDB](https://www.iedb.org/) | Manual export and `tsarina data register iedb …` |
| CEDAR MHC ligand | [CEDAR](https://cedar.iedb.org/) | Manual export and `tsarina data register cedar …` |
| Supported viral proteomes | UniProt | `tsarina data fetch <virus>` |
| Mirrored Tsarina assets | Project data releases | `tsarina data fetch-all` |

The data registry is shared with hitlist. Its default storage location is
`~/.hitlist/`; set `HITLIST_DATA_DIR` to use another location.

IEDB columns are resolved from CSV headers, with known-index fallbacks for
compatible historical exports.

## Inspect peptide observations

Query all indexed MS hits attributed to a gene:

```bash
tsarina hits --gene PRAME
```

Filter the observations and include reference-level aggregation:

```bash
tsarina hits \
  --gene PRAME \
  --allele 'HLA-A*24:02' \
  --mhc-class I \
  --format refs
```

Inspect healthy somatic tissue observations separately:

```bash
tsarina hits --gene PRAME --healthy-tissue
```

The cached observations path preserves multi-gene attribution for shared
peptides. Use `--mono-allelic-only` for direct single-allele evidence or
`--include-binding-assays` to union binding-assay rows with MS observations.

Run `tsarina hits --help` for output formats, resolution filters, raw-export
overrides, and the theoretical peptide-enumeration mode.

## Tissue definitions

Tsarina re-exports oncoref's canonical reproductive-tissue sets for downstream
analysis:

```python
from tsarina.tissues import (
    CORE_REPRODUCTIVE_TISSUES,
    EXTENDED_REPRODUCTIVE_TISSUES,
    PERMISSIVE_REPRODUCTIVE_TISSUES,
    adaptive_rna_threshold,
    is_tissue_restricted,
)
```

The sets progress from core reproductive tissues to broader definitions. The
helper excludes thymus by default when deciding whether detected tissues are
confined to the allowed set.

These re-exports do not create an independent tissue definition. CTA
restriction decisions and their underlying HPA evidence remain owned by
oncoref; see [CTA ownership and downstream evidence](curation.md).

## Presentation scoring

Score an explicit peptide set against a named HLA panel:

```python
from tsarina.alleles import get_panel
from tsarina.scoring import score_presentation

scores = score_presentation(
    peptides=["SLYNTVATL", "GILGFVFTL"],
    alleles=get_panel("iedb27_ab"),
)
```

The main selection workflows interpret lower presentation percentile as
better. Personalized tier thresholds are documented in
[Personalized target selection](personalized-targets.md); panel evidence
thresholds are documented in [CTA panel design](panel-design.md).

## Output identity and naming

All target categories use the same `source`/`source_detail` relationship:

| Category | `source` | `source_detail` | Example |
|---|---|---|---|
| CTA | Gene symbol | Ensembl gene ID | `MAGEA4` / `ENSG00000147381` |
| Viral | Virus short name | UniProt protein accession | `HPV-16` / `P03126` |
| Mutant | Mutation label | Mutation string | `KRAS G12D` / `G12D` |

Keep both fields when joining or exporting results. `source` is the
human-readable target label; `source_detail` preserves the underlying
identifier needed for provenance.
