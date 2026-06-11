# Cancer-Testis Antigen (CTA) Curation

This document describes how the tsarina CTA gene set is built, filtered, and maintained.

## Overview

Cancer-testis antigens (CTAs) are proteins normally restricted to reproductive tissues (testis, ovary, placenta) that become aberrantly expressed in tumors. Their tissue restriction makes them attractive targets for cancer immunotherapy because immune responses against them should not damage normal somatic tissues.

The tsarina CTA set is built as an **unbiased union** from multiple published CT antigen databases, then **systematically filtered** using Human Protein Atlas tissue expression data to retain only genes with reproductive-restricted expression.

## Figures

### Source overlap
![CTA Source Venn Diagram](cta-source-venn.png)

### Filter funnel by source
![CTA Filter Funnel](cta-filter-funnel.png)

### Filter outcome by source
![CTA Filter Outcome](cta-filter-outcome.png)

### Deflated reproductive fraction distribution
![Deflated Fraction Distribution](cta-deflated-frac-dist.png)

### Protein reliability vs RNA fraction
![Protein vs RNA](cta-protein-vs-rna.png)

## Source databases

### CTpedia (167 genes)

The [CTdatabase/CTpedia](http://www.cta.lncc.br/) is the foundational CT antigen reference, maintained by the Ludwig Institute for Cancer Research.

- **Publication**: [Almeida et al. 2009, *Nucleic Acids Research*](https://doi.org/10.1093/nar/gkn673)
- **Content**: ~279 CT gene families curated from literature, with expression data, immune response data, and PubMed links
- **Our subset**: 167 genes from CTpedia that had HPA tissue antibody staining restricted to testis and placenta (original pirlygenes CTA list)

### CTexploreR / CTdata (62 new genes)

[CTexploreR](https://www.bioconductor.org/packages/release/bioc/html/CTexploreR.html) is a modern Bioconductor R package providing a curated, actively maintained CT gene list based on GTEx, CCLE, TCGA, and ENCODE data.

- **Publication**: [Loriot et al. 2025, *PLOS Genetics*](https://doi.org/10.1371/journal.pgen.1011734)
- **Content**: 280 genes classified as CT_gene (146, testis-specific) or CTP_gene (134, testis-preferential)
- **Classification criteria**: Testis-specific genes must have low expression in all somatic tissues and >=10x higher expression in testis. Testis-preferential genes must have low expression in >=75% of somatic tissues with >=10x testis enrichment. Both must show activation in cancer cell lines (CCLE) and tumors (TCGA).
- **Our subset**: 62 protein-coding CTexploreR genes not already in CTpedia that pass our HPA filter. Source tagged as `CTexploreR_CT` or `CTexploreR_CTP`.
- **GitHub**: [UCLouvain-CBIO/CTdata](https://github.com/UCLouvain-CBIO/CTdata)

### da Silva et al. 2017 protein-level CT genes (89 new genes)

[da Silva et al. 2017](https://doi.org/10.18632/oncotarget.21715) performed genome-wide identification of CT genes using RNA-seq across multiple normal tissue datasets and validated a subset at the protein level using tumor mass spectrometry proteomics.

- **Publication**: [da Silva et al. 2017, *Oncotarget*](https://doi.org/10.18632/oncotarget.21715)
- **Full predicted set**: 1,103 testis-biased genes (proportional score >= 0.9)
- **Protein-level subset**: 136 genes detected by mass spectrometry in tumor samples (melanoma, colorectal, breast, prostate, ovary — 209 samples total). These are genes with direct evidence of protein expression in tumors.
- **Our subset**: 135 of the 136 protein-level genes (1 excluded: SPATA8 reclassified as lncRNA in Ensembl 112). 46 overlap existing genes, 89 new. Source tagged as `daSilva2017_protein`.
- **Note**: The full 1,103-gene set uses a broad 0.9 proportional score threshold. The paper also describes a stricter 0.99 threshold (478 genes) and S-score tumor expression filter (418 genes), but these subsets are not provided as downloadable supplementary data. We use only the 136 protein-level genes, which represent the highest-confidence subset with direct proteomic evidence.

### EWSR1-FLI1 CT gene binding sites (12 genes)

CT genes identified as having EWSR1-FLI1 transcription factor binding sites in Ewing sarcoma.

- **Publication**: [Gallegos et al. 2019, *Molecular and Cellular Biology*](https://doi.org/10.1128/MCB.00138-19)
- **Content**: 40 CT genes from Table 2. 15 already in CTpedia, 12 new additions, 13 broadly expressed (excluded).

### Literature scan (28 genes)

Testis-specific genes from meiosis, piRNA pathway, and spermatogenesis literature that pass the HPA reproductive-tissue filter. These include well-characterized testis-specific gene families not captured by the databases above:

- Synaptonemal complex: SYCP3, SMC1B, RAD21L1, SYCE2, MEIOB, FKBP6
- piRNA pathway: PIWIL1, MAEL, DDX4
- Spermatogenesis: BRDT, LDHC, BOLL, NANOS2, ZPBP, ZPBP2, CALR3, ACTL7A, ACTL7B, DMRTB1
- Pluripotency: DPPA3, DPPA5, UTF1
- Known CT antigens: MAGEA8, MAGEA12, MAGEB10, GAGE1, PASD1, TEX14

## Other CT databases considered

Several additional databases were evaluated but not used as primary sources:

| Database | Genes | Why not included directly |
|---|---|---|
| [da Silva et al. 2017](https://doi.org/10.18632/oncotarget.21715) full set | 1,103 | Broad 0.9 threshold; only the 136 protein-level subset used |
| [Wang et al. 2016](https://doi.org/10.1038/ncomms10499) | 876 | Genome-wide screen; supplementary gene lists overlap with other sources |
| [Bruggeman et al. 2018](https://doi.org/10.1038/s41388-018-0357-2) | 756 | Germ cell-specific genes; broader than classic CTA definition |
| [MSigDB YOKOE set](https://www.gsea-msigdb.org/gsea/msigdb/cards/YOKOE_CANCER_TESTIS_ANTIGENS) | 35 | Small, fully covered by CTpedia |
| [HPA testis-elevated](https://www.proteinatlas.org/humanproteome/tissue/testis) | 1,994 | Testis-elevated, not CTA-specific |

Genes from these sources are cross-referenced in the `source_databases` column (e.g., `daSilva2017` tag indicates the gene appears in the da Silva full 1,103-gene set).

## Scope and provenance (tsarina#79)

**Protein-coding only.** Every CTA in the table is `biotype == protein_coding`;
the filter rejects non-coding genes outright (see [Filter logic](#filter-logic)).
So the lncRNA / antisense CT genes in CTexploreR's `CT_genes` list
(LINC*, *-AS1, *-DT, …) are **deliberately out of scope** — tsarina curates
protein-coding pMHC targets, not the full CT transcriptome. This is the
resolution of the lncRNA/antisense scope question in tsarina#79.

**Every gene carries source provenance.** The `source_databases` column is
populated for all genes — there are **no sourceless or "literature-only"
entries**. Each token records membership in a published source:

| `source_databases` token | meaning |
|---|---|
| `CTpedia` | CTdatabase / CTpedia (curated CT antigen reference) |
| `CTexploreR_CT` / `CTexploreR_CTP` | CTexploreR testis-specific / testis-preferential |
| `daSilva2017` | da Silva 2017 full 1,103-gene RNA-seq screen |
| `daSilva2017_protein` | da Silva 2017 tumor-MS protein-level subset |
| `paralog:<sibling>` | near-identical paralog copy of the named curated CTA, added by tsarina's paralog reasoning (not in a source DB) — see below |

Provenance is therefore **queryable directly**: a gene in a hand-curated CT
database vs. one surfaced only by the high-throughput da Silva screen are
distinguishable by token. Current split (382-gene universe): **283** genes are
in a curated CT database (CTpedia and/or CTexploreR); **71** are screen-only
(`daSilva2017`); **10** are `paralog:` copies; and **18** are `placental_antigen`
(see below). To select by provenance:

```python
from tsarina import CTA_evidence
df = CTA_evidence()
src = df["source_databases"]
curated = df[src.str.contains("CTpedia|CTexploreR")]               # 283
paralog = df[src.str.contains("paralog")]                         # 10
placental = df[src.str.contains("placental_antigen")]            # 18
screen_only = df[src.str.contains("daSilva") & ~src.str.contains("CTpedia|CTexploreR")]  # 71
```

(The literature/EWSR1-FLI1 sources listed above confirmed genes already present
in the curated databases — none entered the panel by literature scan alone.)

### Non-CTA exclusions (conserved/multicopy families)

A few genes enter the CTA universe via a source database but are **not** cancer-
testis antigens — conserved/multicopy housekeeping families. They can even pass
the reproductive-restriction filter on a testis-enriched copy, yet they make
poor targets (ubiquitous/essential proteins) and pollute any CTA-keyed
sequence/expression analysis. They are dropped from the universe at load time
via `tsarina.tissues.NON_CTA_EXCLUDED_GENE_IDS` (tsarina#92):

| excluded | family | why not a CTA |
|---|---|---|
| `H4C6`, `H2BC1`, `H2BC3` | core histones (H4, H2B) | conserved, multicopy nucleosome core |
| `H1-1` | somatic linker histone H1.1 (HIST1H1A) | conserved/multicopy; somatic, not testis-specific |
| `TUBA3C`, `TUBA3E` | alpha-tubulins | essential cytoskeleton |

Testis-specific histone *variants* (e.g. `H1-6` / HIST1H1T) are genuine CTAs and
are **kept** — the exclusion is a curated per-gene deny-list, not a blanket
family filter. (`CGB8` / hCG-beta was formerly on this deny-list as a placental
gene, but is now curated as a **placental onco-germline antigen** — see below.)

### Paralog copies (tsarina#93)

Some CTA families have near-identical paralog copies that aren't in any source
database. Left out of the universe, their (sequence-identical) protein lands in
the *non-CTA* set of downstream sequence analyses and cancels the family's
specificity. The testis-restricted copies are added to the universe tagged
`paralog:<curated sibling>` (e.g. `paralog:CT47A1`): `CT47A8/A9/A10`,
`CT45A8/A9` (`paralog:CT45A1`), `GAGE12D` (`paralog:GAGE12C`), `GAGE10`
(`paralog:GAGE2A`). Most are `never_expressed` (very low HPA RNA) so they sit in
the universe/filtered sets but not the expressed default.

A paralog copy is added only if its annotated protein is real: the
`add_cta_gene.py` guard refuses **fragment gene models** (tsarina#108). For
example `GAGE12B` (`ENSG00000236737`) is held out — in Ensembl it is a degenerate
117-bp, single-exon locus encoding a 7-residue ORF (the shared GAGE C-terminus),
not a GAGE protein, so its high HPA testis nTPM is quantification noise on a
malformed model rather than evidence of a real antigen.

Three more — `MAGEA2B`, `CT45A5`, `SSX4B` — were already in the universe **by
Ensembl ID** under a different symbol (`MAGEA2`, `CT45A10`, `SSX4`); they are
added as aliases so symbol lookups resolve. Three somatic-on-bulk-HPA copies
(`CT45A6`, `DAZ2`, `DAZ4`) are held out — their bulk signal (likely DAZ
cross-mapping at a hard-to-map multicopy Y locus) fails the restriction filter.
`DAZ1` (already curated) is testis-restricted on bulk; its reported gastric
single-cell signal is below the bulk resolution tsarina curates on, and it is
already `never_expressed`.

### Placental onco-germline antigens (tsarina#111)

Cancer-testis antigens are the *testicular* members of the cancer-germline
family; the **placenta** is the other immune-privileged reproductive tissue with
the same logic. 18 placenta-restricted antigen genes — absent from the
testis-oriented source databases but documented onco-placental tumor antigens —
are added to the universe tagged `placental_antigen`:

| family | genes | placental product |
|---|---|---|
| hCG-beta | `CGB1`, `CGB2`, `CGB3`, `CGB5`, `CGB7` | chorionic gonadotropin beta |
| pregnancy-specific glycoproteins | `PSG2`, `PSG4`, `PSG6`, `PSG7` | PSG family |
| endogenous retrovirus envelopes | `ERVW-1` (syncytin-1), `ERVFRD-1` (syncytin-2), `ERVH48-1` (suppressyn), `ERVV-1`, `ERVV-2` | placental fusogens |
| placental galectins | `LGALS13` (PP13), `LGALS14` (PP14) | placental protein 13/14 |
| placental lactogen / GH | `CSH1`, `GH2` | placental hormones |

They enter the universe and are then **judged by the same reproductive-restriction
filter** — the placenta is a core reproductive tissue, so genuinely placenta-
restricted members pass, while the somatically-broad members (e.g. `CSH1`, `GH2`)
self-exclude on the data rather than by hand. All are full-length protein-coding
(139–538 aa; the fragment-model guard above applies here too). `CGB8` was moved
out of the non-CTA deny-list into this set.

### Identical-protein groups (proteoforms)

Several CTA families have multiple gene loci that encode a **byte-identical
protein** — e.g. `CTAG1A`/`CTAG1B` (NY-ESO-1), `SSX4`/`SSX4B`, the 12-member
`CT47A` family. At the pMHC level these are one antigen: their peptides are
indistinguishable, so scoring them separately would double-count the same target
and waste panel slots.

tsarina ships a registry of these groups in `tsarina/data/proteoform-groups.csv`,
**keyed by Ensembl gene ID** (`proteoform_id`, `member_symbol`, `member_gene_id`,
`protein_length`, `n_members`). It is **computed from Ensembl protein sequences**
via `cancerdata.proteoform_groups()` and synced with
`scripts/sync_proteoform_groups.py` — currently **15 families**. The four anchor
groups (`CTAG1A/CTAG1B`, `XAGE1A/XAGE1B`, `SSX4/SSX4B`, `MAGEA2/MAGEA2B`) are
hard-coded as an offline fallback so grouping never silently disappears. Panel
construction collapses each group to a single target
(`--group-identical-cta-peptide-sets` / `--group-identical-cta-pmhcs`, on by
default), and the group label resolves common aliases (e.g. `NY-ESO-1`).

This identical-protein grouping is **CTA-internal**. Cross-mapping against
*normal non-CTA* proteins is handled separately and at the peptide level: a CTA
peptide that also occurs in **any** non-CTA protein is dropped by the
CTA-exclusivity filter (`cta_exclusive_peptides()`, which streams the candidate
peptides against the non-CTA proteome k-mer set), rather than by grouping
non-CTA genes.

## HPA tissue expression annotation

Every gene is scored against [Human Protein Atlas](https://www.proteinatlas.org/) v23 using two data modalities.

> **Reproducibility / provenance.** Both HPA modalities are pinned to release
> **v23** — the most recent HPA version whose download mirror serves *both*
> `rna_tissue_consensus` (RNA) and `normal_tissue` (IHC) as a matched pair, so
> the whole table is derived from one release. Fetch and inspect the sources
> with `tsarina reference {list,fetch}`; regenerate every RNA/protein
> column for the bundled table with `python scripts/regenerate_table.py`
> (`--apply` to write). A small curated `_CROSS_REACTIVE_IHC` override in that
> script forces the IHC to "no data" for sequence-near-identical paralog
> antigens (MAGEB6, XAGE2, CT45A8, CT45A9) whose shared antibody cross-reacts
> (Low-level protein in tissues with ~0 nTPM RNA).

### RNA expression

**Data source**: [HPA RNA tissue consensus](https://www.proteinatlas.org/about/download) (`rna_tissue_consensus.tsv`)

This dataset provides normalized transcripts per million (nTPM) values across **50 normal human tissues**, representing a consensus of RNA-seq data from HPA, GTEx, and FANTOM5.

**Core reproductive tissues**: testis, ovary, placenta

**Thymus exclusion**: Thymus is excluded from all restriction calculations because AIRE (autoimmune regulator) drives ectopic expression of tissue-restricted antigens in medullary thymic epithelial cells (mTECs) as part of central immune tolerance. CTA expression in thymus is expected and does not indicate somatic tissue leakage.

**Deflated reproductive fraction**: To suppress low-level basal transcription noise, we compute a deflated metric:

```
deflated_fraction = (1 + sum_reproductive(max(0, nTPM - 1))) / (1 + sum_all(max(0, nTPM - 1)))
```

- `max(0, nTPM - 1)` zeros out sub-1 nTPM values (below HPA's own detection threshold)
- The `+1` pseudocount on numerator and denominator prevents 0/0 for very-low-expression genes where all tissues have nTPM < 1
- Thymus is excluded from the denominator (sum_all)

**Example**: CTCFL/BORIS has testis nTPM = 10.8 but ~40 other tissues at 0.1-0.9 nTPM each. Raw reproductive fraction: 54%. Deflated fraction: 100% (only testis exceeds 1 nTPM, so all other tissues contribute 0 after deflation).

### Protein expression

**Data source**: [HPA normal tissue IHC](https://www.proteinatlas.org/about/download) (`normal_tissue.tsv`)

This dataset provides immunohistochemistry (IHC) staining levels across **63 normal human tissues** with antibody reliability scores.

**Detection levels**: Not detected, Low, Medium, High

**Antibody reliability** (highest to lowest confidence):
- **Enhanced**: Orthogonal validation (mass spectrometry, Western blot, or similar)
- **Supported**: Staining consistent with gene/protein characterization
- **Approved**: Basic validation passed (at least one cell type detected as expected)
- **Uncertain**: Contradictory or unreliable staining pattern

A gene's `protein_reproductive` flag is True when all tissues with detected protein (excluding thymus) are in {testis, ovary, placenta}. The `protein_reliability` column reports the best (highest confidence) reliability score across all antibodies for that gene.

## Filter logic

The `passes_filters` column uses tiered deflated RNA reproductive fraction
thresholds that scale with protein data confidence. Higher-confidence protein
data in reproductive tissues provides corroborating evidence, allowing a more
permissive RNA threshold:

| Protein evidence | Required deflated RNA fraction |
|---|---|
| Enhanced (orthogonal validation) + reproductive only | >= 80% |
| Supported (consistent characterization) + reproductive only | >= 90% |
| Approved (basic validation) + reproductive only | >= 95% |
| Uncertain or no protein data available | >= 97% |

The thresholds live in `tsarina.tissues.HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS`.
The no-protein/`Uncertain` gate was relaxed 0.99 → 0.98 (tsarina#51) and then
0.98 → 0.97 to admit borderline reproductive-restricted genes whose dominant
tissue is strongly expressed but which carry a low-level somatic signal — e.g.
**XAGE2** (placenta 187 nTPM, lung 5.3; tsarina#79), **CT83**/KK-LC-1 (62 nTPM,
salivary 2.7), and **PRM3** (86 nTPM, bone marrow 2.6). Such genes keep a
visible `safety_flags` entry for the somatic tissue.

**Additional filter criteria**:
- Gene must be protein-coding (Ensembl biotype = `protein_coding`)
- Genes with protein detected in **non-reproductive tissues** (excluding thymus) always fail, regardless of RNA fraction
- Thymus is excluded from both RNA and protein restriction checks

## Three-axis CTA tier system

`tiers.assign_all_axes` classifies every gene along three independent axes:
**`restriction`** (what tissues), **`restriction_confidence`** (how confident),
and **`ms_restriction`** (MS evidence in healthy tissue). Counts below are over
the full 382-gene table.

### Axis 1: `restriction` — tissue restriction category

`tiers.synthesize_restriction` takes the IHC protein restriction
(`protein_restriction`) when protein data exists, otherwise falls back to the
RNA restriction (`rna_restriction`).

| Value | Count | Meaning |
|-------|-------|---------|
| `TESTIS` | 265 | Only testis detected (the most clinically actionable — blood-testis barrier sequesters these from circulation) |
| `REPRODUCTIVE` | 25 | Multiple reproductive tissues (testis + ovary/placenta) |
| `PLACENTAL` | 10 | Only placenta (± testis) |
| `SOMATIC` | 73 | Detected in a non-reproductive somatic tissue — a leakage flag; these typically **fail** the restriction filter |
| `NO_DATA` | 9 | No protein and no usable RNA restriction call |

### Axis 2: `restriction_confidence` — confidence in the restriction call

`synthesize_restriction` averages a per-source evidence score (protein IHC, RNA,
MS) and bins the mean: **≥ 1.2 → HIGH, ≥ 0.8 → MODERATE, else LOW**. So
`MODERATE` ≈ one solid uncorroborated source.

| Value | Count | Meaning |
|-------|-------|---------|
| `HIGH` | 146 | Strong / corroborated evidence for the restriction call |
| `MODERATE` | 163 | One solid source, uncorroborated |
| `LOW` | 64 | Weak evidence (e.g. healthy-tissue MS) |
| `NO_DATA` | 9 | No usable evidence |

The confidence is **capped at MODERATE** when the only evidence is RNA below the
expression floor (`never_expressed`, no protein, no MS) — the scorer otherwise
grants the same STRICT credit to near-noise RNA as to robust expression, so such
a call would over-state HIGH (tsarina#114). Genes with protein or MS evidence are
untouched. Note the axis measures restriction *cleanliness*, not target validity
or abundance — use it as a ranking signal, not a hard filter (see below).

### Axis 3: `ms_restriction` — MS evidence classification

Panel construction recomputes public MHC-ligand MS support from the current
hitlist observations index. In the default cached path, hitlist builds that
index from registered IEDB and CEDAR exports plus hitlist's manually curated
supplementary MS rows. These panel-time counts do not use hitlist's bulk
proteomics or line-expression indexes. Passing explicit IEDB/CEDAR paths uses
the raw export scanner for those files instead of the cached observations index.

| Value | Criteria |
|-------|----------|
| `CANCER_ONLY` | Gene's peptides found only in cancer MS; zero healthy somatic tissue hits |
| `EXPECTED_TISSUE` | MS hits only in reproductive/thymic tissue (expected for CTAs) |
| `SINGLETON_HEALTHY` | 1 peptide in ≤ 1 healthy somatic tissue (possible noise) |
| `RECURRENT_HEALTHY` | Multiple peptides or tissues in healthy somatic MS (genuine off-target) |
| `UNCLASSIFIED_MS` | MS evidence exists but every observation is from an unclassifiable source (e.g. cell-line only) — distinct from no data |
| `NO_MS_DATA` | No MS evidence available for this gene's peptides |

### Published CTAs across the axes

The tiers describe evidence **quality**, not target **validity**.  A naive hard
gate of `restriction ∈ {TESTIS, PLACENTAL}` × `confidence ≥ MODERATE` would
reject clinically-validated targets that are in the expressed (POSITIVE) set:

| gene | in `CTA_gene_names()` | restriction | confidence | ms_restriction | hard gate |
|---|:---:|---|---|---|:---:|
| NY-ESO-1 (CTAG1B) | ✅ | TESTIS | HIGH | CANCER_ONLY | ✅ admit |
| NUTM1 | ✅ | TESTIS | HIGH | CANCER_ONLY | ✅ admit |
| HORMAD1 | ✅ | TESTIS | MODERATE | CANCER_ONLY | ✅ admit |
| PAGE2 | ✅ | TESTIS | MODERATE | CANCER_ONLY | ✅ admit |
| PAGE5 | ✅ | TESTIS | MODERATE | CANCER_ONLY | ✅ admit |
| **MAGE-A4** | ✅ | **REPRODUCTIVE** | MODERATE | RECURRENT_HEALTHY | ❌ reject |
| **PRAME** | ✅ | TESTIS | **LOW** | SINGLETON_HEALTHY | ❌ reject |
| **XAGE1A** | ✅ | **SOMATIC** | MODERATE | CANCER_ONLY | ❌ reject |

MAGE-A4 (FDA-approved afami-cel/Tecelra), PRAME (e.g. IMA203), and XAGE1A are all
real targets the gate would drop.  (CTAG2, the third NY-ESO-1 paralog, correctly
lands in the somatic-leak held-out class — `CTA_excluded_*`.)  **Use the tiers as
ranking signals, not hard filters.**

### Peptide-level healthy-tissue screening

The gene-level `ms_restriction` collapses detail that matters for safety.
`cta_healthy_tissue_ms_hits(gene)` surfaces the underlying healthy-*somatic* MS
hits at peptide × tissue × allele granularity (reproductive/thymic hits, which
are expected for CTAs, are excluded):

```python
from tsarina import cta_healthy_tissue_ms_hits

cta_healthy_tissue_ms_hits("PRAME")
#   peptide tissue       allele provenance      pmid  vital_organ
# SQLTTLSFY  Blood  HLA class I  unmatched  29557506        False
# -> cautious gene tier (SINGLETON_HEALTHY), but clean: 1 blood hit, 0 vital organs.

cta_healthy_tissue_ms_hits("MAGEA4")
# -> 3 heart-eluted peptides (AETSYVKV family, HLA-B*49:01; PMID 33858848,
#    vital_organ=True) + 1 blood hit. The clinical epitope GVYDGREHTV (A*02:01)
#    is seen only in cancer, never on heart.
```

`vital_organ` flags tissues matching a `SAFETY_TISSUE_GROUPS` vital organ
(brain / heart / lung / liver / pancreas) — the "which organ / which allele /
which peptide" detail a per-patient screen needs.

The vital-organ vocabulary is exported for downstream consumers (e.g. vaxrank
safety scoring): `tsarina.SAFETY_TISSUE_GROUPS` (HPA tissue names, grouped) and
`tsarina.VITAL_TISSUE_MS_NAMES` (the MS source-tissue spellings).

The same screen is available from the CLI:

```bash
tsarina hits --gene MAGEA4 --mhc-class I --healthy-tissue
# peptide,tissue,allele,allele_set,provenance,pmid,vital_organ
# AETSYVKV,Heart,HLA-B*49:01,...,33858848,True
```

The packaged CTA evidence table intentionally does not include MS count columns.
Use `CTA_detailed_evidence()` or `tsarina panel` to recompute current
hitlist-derived counts. The runtime `ms_cta_exclusive_*` counts use the stricter
CTA-exclusive peptide set: peptides found in any non-clean-CTA protein are not
counted there.

Automatic panel ranking now uses bundled HPA cancer RNA prevalence by default
(`tumor_prevalence_panel_score`): cancer-type breadth at pTPM >= 2.0 and a 5%
sample-prevalence floor, then sample-level prevalence, with HPA cancer IHC as a
weak tie-breaker. Public-MS support/safety is then recomputed from hitlist for
candidate batches before pMHC scoring. This keeps MS evidence live without
bundling derived count snapshots.

The bundled HPA cancer prevalence tables can be regenerated for an updated CTA
candidate CSV with:

```bash
python scripts/regenerate_hpa_cancer_prevalence.py \
  --cta-csv tsarina/data/cancer-testis-antigens.csv \
  --output-dir tsarina/data
```

The script streams HPA's large `rna_cancer_sample.tsv.gz` table and writes only
per-gene/cancer prevalence summaries for genes in the CTA CSV, plus a compact
subset of HPA's cancer IHC count table.

### Axis-aware API

```python
from tsarina import (
    CTA_testis_restricted_gene_names,     # 247 genes
    CTA_placental_restricted_gene_names,  # 10 genes
    CTA_by_axes,                          # flexible multi-axis filter
    RESTRICTION_VALUES, CONFIDENCE_VALUES, MS_RESTRICTION_VALUES,
)

# Testis-restricted with high-confidence restriction call
strict_testis = CTA_by_axes(restriction="TESTIS", restriction_confidence="HIGH")

# All testis + placental genes (serum biomarkers for non-pregnant)
biomarkers = CTA_by_axes(restriction={"TESTIS", "PLACENTAL"})

# Evidence table — filter by any combination
from tsarina import CTA_evidence
df = CTA_evidence()
testis_high = df[(df["restriction"] == "TESTIS") & (df["restriction_confidence"] == "HIGH")]

# Full per-tissue detail (requires HPA data files)
from tsarina import CTA_detailed_evidence
detailed = CTA_detailed_evidence(hpa_bulk_path="proteinatlas.tsv")
# Adds: rna_testis_ntpm, rna_ovary_ntpm, rna_placenta_ntpm,
#        rna_max_somatic_tissue, rna_max_somatic_ntpm, rna_somatic_detected_count
```

## Never-expressed flag

The `never_expressed` column flags genes where:
- No HPA protein (IHC) data is available, AND
- Maximum RNA nTPM across all tissues is below the expression floor

The floor is the explicit, parameterized constant
`tsarina.tissues.HPA_EXPRESSION_FLOOR_NTPM` (default **2 nTPM**), rather than a
magic number baked into the flag (tsarina#78).

These genes pass the filter (because the `+1` pseudocount gives a 1.0 deflated fraction when all nTPMs are below 1), but the evidence for their tissue restriction is weak — HPA simply doesn't have enough signal to confirm or deny reproductive specificity. They are typically very low-abundance transcripts below HPA's detection sensitivity. Many are still legitimate CTAs supported by other evidence (e.g., CTpedia listing, tumor mass spectrometry detection), but users should be aware of the limited HPA evidence.

Currently **29 genes** are flagged as low evidence.

**Manual rescue**: rather than lowering the global floor (which would admit
paralog cross-mapping noise; see tsarina#78), individual borderline-but-real
CTAs are rescued into the expressed set via
`tsarina.tissues.MANUALLY_EXPRESSED_CTA`. These genes keep `never_expressed = True`
in the table (HPA truth) but are still returned by `CTA_gene_names()`. Currently
rescued: **XAGE5** (`ENSG00000171405`; testis 1.1 nTPM; CTpedia/CTexploreR/daSilva2017).

## Aliases and name resolution

The `Aliases` column is populated for every gene from NCBI `gene_info` synonyms
(matched by Ensembl gene ID), merged with curated colloquial names NCBI omits
(e.g. `NY-ESO-1`).  Regenerate with `python scripts/backfill_aliases.py` (it
downloads + caches the NCBI source; see tsarina#77).

`cta_symbol_for_alias(name)` resolves a CTA name or synonym to its official
symbol, case- and punctuation-insensitive:

```python
from tsarina import cta_symbol_for_alias
cta_symbol_for_alias("NY-ESO-1")   # -> "CTAG1B"
cta_symbol_for_alias("ESO1")       # -> "CTAG1B"
cta_symbol_for_alias("CT12.2")     # -> "XAGE2"
cta_symbol_for_alias("not-a-gene") # -> None
```

So callers querying CTAs by common antigen names rather than HGNC symbols no
longer silently miss hits.

## Gene symbol maintenance

Gene symbols are updated to current HGNC nomenclature, with old symbols preserved in the `Aliases` column. Known renames:

| Old symbol | Current symbol | Reason |
|---|---|---|
| TSPY9P | TSPY9 | Reclassified from pseudogene to protein-coding |
| ODF3 | ODF3 (alias: CIMAP1A) | Renamed in Ensembl 112 |
| TEX33 | TEX33 (alias: CIMIP4) | Renamed in Ensembl 112 |
| TEX37 | TEX37 (alias: SPMIP9) | Renamed in Ensembl 112 |
| THEG | THEG (alias: SPMAP2) | Renamed in Ensembl 112 |
| C17orf104 | TLCD3A | HGNC rename |
| CCDC155 | KASH5 | HGNC rename |
| FAM71E2 | GARIN4 | HGNC rename |
| HIST1H1A | H1-1 | Histone nomenclature update |
| HIST1H1T | H1-6 | Histone nomenclature update |
| HIST1H2BA | H2BC1 | Histone nomenclature update |
| HIST1H2BB | H2BC3 | Histone nomenclature update |
| HIST1H4F | H4C6 | Histone nomenclature update |

All Ensembl Gene IDs are validated against Ensembl release 112. Canonical transcript IDs (longest protein-coding transcript) are provided in the `Canonical_Transcript_ID` column.

## Column reference

| Column | Description |
|---|---|
| `Symbol` | Current HGNC gene symbol |
| `Aliases` | Previous/alternative gene symbols (semicolon-separated) |
| `Full_Name` | Gene full name |
| `Function` | Functional annotation |
| `Ensembl_Gene_ID` | Ensembl gene ID (validated against release 112) |
| `source_databases` | Source databases (CTpedia, CTexploreR_CT, CTexploreR_CTP, daSilva2017, daSilva2017_protein) |
| `protein_reproductive` | IHC detected only in {testis, ovary, placenta} (excl. thymus), or `"no data"` |
| `protein_thymus` | IHC detected in thymus |
| `protein_reliability` | Best HPA antibody reliability (Enhanced / Supported / Approved / Uncertain / `"no data"`) |
| `rna_reproductive` | All tissues with >=1 nTPM (excl. thymus) are in {testis, ovary, placenta} |
| `rna_thymus` | Thymus nTPM >= 1 |
| `protein_strict_expression` | Semicolon-separated tissues with IHC detection (excl. thymus) |
| `rna_reproductive_frac` | Fraction of total nTPM (excl. thymus) in core reproductive tissues |
| `rna_reproductive_and_thymus_frac` | Same, with thymus added to numerator and denominator |
| `rna_deflated_reproductive_frac` | `(1 + sum_repro(max(0, nTPM-1))) / (1 + sum_all(max(0, nTPM-1)))` |
| `rna_deflated_reproductive_and_thymus_frac` | Same, with thymus added to reproductive numerator |
| `Canonical_Transcript_ID` | Longest protein-coding transcript (Ensembl 112) |
| `biotype` | Ensembl gene biotype (must be `protein_coding` to pass filter) |
| `rna_max_ntpm` | Maximum nTPM across all tissues |
| `rna_80_pct_filter` | Deflated reproductive fraction >= 80% |
| `rna_90_pct_filter` | Deflated reproductive fraction >= 90% |
| `rna_95_pct_filter` | Deflated reproductive fraction >= 95% |
| `rna_97_pct_filter` | Deflated reproductive fraction >= 97% (the active no-protein/`Uncertain` gate) |
| `rna_98_pct_filter` | Deflated reproductive fraction >= 98% |
| `rna_99_pct_filter` | Deflated reproductive fraction >= 99% |
| `passes_filters` | Final inclusion flag (see filter logic above) |
| `never_expressed` | No HPA protein data AND max RNA nTPM < 2 |
| `restriction` | Tissue restriction axis: `TESTIS`, `PLACENTAL`, `REPRODUCTIVE`, `SOMATIC`, or `NO_DATA` |
| `rna_restriction` / `rna_restriction_level` | RNA-only restriction call + quality (`STRICT`/`MODERATE`/`PERMISSIVE`/`LEAKY`) |
| `protein_restriction` | IHC-only restriction call (`TESTIS`/`PLACENTAL`/`REPRODUCTIVE`/`SOMATIC`/`NO_DATA`) |
| `restriction_confidence` | Confidence axis: `HIGH`, `MODERATE`, `LOW`, or `NO_DATA` |
| `protein_testis` | IHC protein detected in testis (`True`/`False`/empty if no data) |
| `protein_ovary` | IHC protein detected in ovary (`True`/`False`/empty if no data) |
| `protein_placenta` | IHC protein detected in placenta (`True`/`False`/empty if no data) |
| `ms_restriction` | MS safety axis: `CANCER_ONLY`, `EXPECTED_TISSUE`, `SINGLETON_HEALTHY`, `RECURRENT_HEALTHY`, `UNCLASSIFIED_MS`, or `NO_MS_DATA` |
| `safety_flags` | Semicolon-separated somatic tissue groups with RNA above the safety threshold (visible leakage markers) |

## Python API

```python
from tsarina import (
    CTA_gene_names,                # expressed + filter-passing CTAs (recommended default)
    CTA_gene_ids,                  # same, as Ensembl gene IDs
    CTA_never_expressed_gene_names,# filter-passing but no HPA expression
    CTA_filtered_gene_names,       # all filter-passing (= expressed + never_expressed)
    CTA_excluded_gene_names,       # CTAs that FAIL filter (somatic expression)
    CTA_unfiltered_gene_names,     # full CTA universe (all source databases)
    CTA_evidence,                  # full DataFrame with all evidence columns
)

# Default: expressed, reproductive-restricted CTAs
cta_genes = CTA_gene_names()

# Full CTA universe (for excluding from non-CTA comparison sets)
all_ctas = CTA_unfiltered_gene_names()

# Evidence table — filter however you like
df = CTA_evidence()

# Example: strict CTAs from CTpedia with Enhanced protein evidence
strict = df[
    (df['passes_filters'] == True) &
    (df['source_databases'].str.contains('CTpedia')) &
    (df['protein_reliability'] == 'Enhanced') &
    (~df['never_expressed'])
]

# Example: genes with tumor mass spec evidence
tumor_protein = df[df['source_databases'].str.contains('daSilva2017_protein', na=False)]
```

## Gene partitioning for pMHC analysis

When comparing CTA pMHCs against non-CTA pMHCs, every protein-coding gene needs to go into exactly one bucket. The `CTA_partition_*` functions handle this:

```python
from tsarina.partition import (
    CTA_partition_gene_ids,       # sets of Ensembl gene IDs
    CTA_partition_gene_names,     # sets of gene symbols
    CTA_partition_dataframes,     # DataFrames with evidence columns
)

# Each returns a dataclass with .cta, .cta_never_expressed, .non_cta
p = CTA_partition_gene_ids()
p.cta                   # set of Ensembl IDs for expressed CTAs
p.cta_never_expressed   # set of Ensembl IDs for never-expressed CTAs
p.non_cta               # set of Ensembl IDs for everything else

p = CTA_partition_gene_names()
"MAGEA4" in p.cta       # True
"TP53" in p.non_cta     # True

p = CTA_partition_dataframes()
p.cta.columns           # full evidence columns for CTAs
p.non_cta.columns       # Symbol, Ensembl_Gene_ID
```

| Partition | Description | Typical count |
|---|---|---|
| `p.cta` | Expressed, reproductive-restricted CTAs. Source of CTA pMHCs. | ~272 |
| `p.cta_never_expressed` | CTAs from databases but no meaningful HPA expression (max nTPM < 2, no protein data). Pass filter on a technicality (pseudocount). Separate from analysis. | ~29 |
| `p.non_cta` | All other protein-coding genes, **including** CTAs that fail the reproductive-tissue filter (somatic expression). Clean non-CTA comparison set. | ~19,800 |

These three sets are **non-overlapping** and their union covers all protein-coding genes from Ensembl.

**Why three partitions instead of two?**
- **Never-expressed CTAs** pass our filter because the +1 pseudocount gives them a 1.0 deflated fraction when all nTPMs are below 1. They are in CT antigen databases but HPA has no real signal. Including them in pMHC analysis would add noise — you can't target a protein that's never made.
- **Excluded CTAs** (those that fail the filter due to somatic expression) are folded into `non_cta`. They express in healthy tissue, so their peptides would appear in the non-CTA proteome anyway. Keeping them in `non_cta` gives a realistic comparison set.
