# Ecosystem and Python API

Tsarina is the decision and design layer for shared cancer immunotherapy
targets. Rather than duplicating reference genomes, antigen-restriction calls,
mass-spec curation, or MHC predictors, it composes focused packages from the
PIRL / OpenVax stack.

## How Tsarina connects to the stack

```text
  oncoref (gene-level specificity, human/mouse/dog CTAs, 4-way sex restriction, fusions, cancer burden)
       +
  pyensembl / orfdata (reference coding proteomes; IARC GDC & RefSeq viral ORFs)
       +
  hitlist (public MS observations, healthy/cancer/sex-repro tissue flags, EvidenceStore)
       +
  mhcgnomes (HLA / H-2 / DLA allele, serotype, and species normalization)
       +
  mhcflurry / mhctools / topiary / pepsickle / dnachisel (presentation, cleavage, mRNA encoding)
       |
       v
    tsarina
       ├── personalize  --> ranked CTA, fusion, mutant, viral pMHCs for one patient
       ├── panel        --> population CTA x HLA matrix with MS tiers & coverage
       ├── hits         --> gene/protein MS & healthy-tissue safety lookup
       └── vaccine      --> mortality-ranked CTA segment selection & DNA/RNA construct
```

| Package | Role in the stack | What Tsarina reads from it |
|---|---|---|
| [oncoref](https://github.com/pirl-unc/oncoref) | Single authority for gene-level cancer specificity | Evaluated human, mouse & dog `CTA_selection`, 4-way sex-restricted eligibility (`CTA-M`/`F`/`MF`/`P` & `RA-*-only`), clinical clearance, 1-mm near-self coverage, proteoform groups, 188 curated cancer fusions, TCGA/HPA/GLOBOCAN burden |
| [hitlist](https://github.com/pirl-unc/hitlist) | Public immunopeptidomics MS curation & indexing | Normalized IEDB/CEDAR/study observations, cancer vs. somatic vs. sex-specific reproductive tissue flags, verified `EvidenceStore` packages, dataset registry & preflight |
| [pyensembl](https://github.com/openvax/pyensembl) | Ensembl reference proteomes & species registry | Human (Ensembl 112), mouse (Ensembl 115), and dog (Ensembl 98, CanFam3.1) coding sequences for peptide enumeration and self-proteome exclusivity screening |
| [orfdata](https://github.com/openvax/orfdata) *(optional)* | Viral ORF & GDC reference collections | Three offline IARC GDC oncogenic viral subsets (`cancer-associated` 37, `cancer-causing` 20, `cancer-integrated-driver` 15) and spliced RefSeq viral CDSs (e.g. HPV16 `E1^E4`, `E6*`) |
| [mhcgnomes](https://github.com/openvax/mhcgnomes) | MHC nomenclature & ontology | Canonical parsing of HLA, H-2, and DLA alleles, serotypes, donor sets, and species validation |
| [mhcflurry](https://github.com/openvax/mhcflurry) / [mhctools](https://github.com/openvax/mhctools) / [topiary](https://github.com/openvax/topiary) | MHC presentation & affinity prediction | Presentation percentile, presentation score, and affinity calibration |
| [pepsickle](https://github.com/pirl-unc/pepsickle) / [dnachisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel) *(vaccine extra)* | Proteasomal cleavage & sequence optimization | Junction cleavage audit and constrained DNA/mRNA construct assembly |

## CLI command map

| Command | Purpose |
|---|---|
| `tsarina status` | Read-only check of package versions, OncoRef policy, Hitlist index/EvidenceStore, Ensembl caches, and MHCflurry weights |
| `tsarina catalog {ctas,panels,mutations,fusions,viruses}` | Browse built-in CTA (human/mouse/dog, optional `--sex`), HLA panel, hotspot mutation, OncoRef fusion, and viral catalogs |
| `tsarina data {available,list,fetch,preflight,register,info,path,remove}` | Manage IEDB, CEDAR, viral proteomes, and mirrored datasets (`--grouped` / `--category` supported on `available`) |
| `tsarina build observations` | Build or refresh the local Hitlist observations index |
| `tsarina personalize` | Prioritize CTA, fusion, mutant, and viral pMHCs for one patient (human, mouse, or dog, optional `--sex`) |
| `tsarina hits` | Inspect public MS hits or `--healthy-tissue` safety screens for a gene or UniProt ID |
| `tsarina panel` | Build a CTA × HLA pMHC matrix with evidence tiers, MS safety flags, and population coverage |
| `tsarina vaccine` / `tsarina vaccine-report` | Design a shared CTA vaccine construct and render an interactive report website |

## Python API map

`import tsarina` uses lazy PEP 562 attribute resolution so startup stays fast.
Although `tsarina.__all__` preserves full backward compatibility for historical
gene-set helpers, most workflows only need the functions grouped below.

### Headline workflows

| Entry point | Purpose |
|---|---|
| `personalized_targets(...)` | Rank CTA, fusion, mutant, and viral pMHCs for one patient (human, mouse, or dog, optional `patient_sex` RNA + MS safety) |
| `spanning_pmhc_set(...)` | Build a population CTA × HLA panel with evidence tiers, safety flags, Cys MS-detectability caveats, and coverage |
| `design_vaccine(...)` / `VaccineConfig` | Mortality-weighted CTA segment selection, junction optimization, and audited DNA/RNA construct design |
| `cta_healthy_tissue_ms_hits(...)` | Per-peptide healthy somatic (and optional same-sex adult reproductive) MS safety screen for one gene |

### Target and allele catalogs

| Entry point | Purpose |
|---|---|
| `CTA_selection(category=..., species=..., sex_code=..., ectopic_only=..., replicated_only=..., host_context=..., filter_host=...)` | Current OncoRef CTA catalog across evaluated `human`, `mouse`, and `dog`, with 4-way `sex_code`, ectopic/replicated filters, and `HostContext` assessment/filtering |
| `CTA_selection_coverage(species=...)` | Species readiness (`evaluated`), anatomy, category counts, and 4-way `sex_restricted_counts` |
| `CTA_sex_restricted_antigens(sex=..., eligibility=..., sex_code=..., species=...)` | 4-way biological sex eligibility (`all_adult` or `sex_organ_only`) and `CTA-M`/`F`/`MF`/`P` + `RA-*-only` filtering for `human`, `mouse`, and `dog` |
| `CTA_sex_eligible_gene_names(...)` / `CTA_sex_eligible_gene_ids(...)` | Eligible gene symbols or Ensembl IDs for a given `sex`, `eligibility`, `sex_code`, and `species` |
| `HostContext` / `build_host_context(...)` / `host_compartments(species=...)` / `SEX_RESTRICTION_CODES` | Composable OncoRef host-compartment protection (`pregnant`, `lactating`, `protect_tissues`, `accept_loss_tissues`, `absent_tissues`, `allow_incomplete_scope`) and canonical per-species compartment discovery |
| `annotate_cta_safety_caveats(...)` / `screen_near_self_1mm(...)` / `annotate_near_self_1mm(...)` | OncoRef 1.55+ CTA safety caveats and 1-aa Hamming near-self screening against OncoRef `Unsafe` / `Unsafe\|Other` |
| `cancer_fusions(...)` / `fusion_junction_peptides(...)` | Filter OncoRef's 188 curated cancer fusions and enumerate host-exclusive breakpoint-spanning peptides |
| `HOTSPOT_MUTATIONS` / `mutant_peptides(...)` | Curated driver hotspot mutations and host-exclusive mutant-spanning peptides |
| `viral_peptides(...)` / `orfdata_viral_subsets()` / `orfdata_viral_peptides(...)` | Built-in UniProt oncogenic viruses and optional `orfdata` IARC GDC viral subsets / RefSeq viral ORFs |
| `panel_names()` / `get_panel(name)` | Built-in population HLA class I panels (`iedb27_ab` through `global54_abc`) |

### MS evidence and presentation scoring

| Entry point | Purpose |
|---|---|
| `attribute_ms_evidence(...)` / `summarize_ms_attribution(...)` | Separate exact requested-MHC, donor-set ambiguous, predicted-restriction, and coarse MS observations |
| `annotate_ms_detectability(df)` | Flag zero-MS cysteine / C-terminal Cys ligands (e.g. NY-ESO-1 `SLLMWITQC`) without altering tiers |
| `score_presentation(...)` / `score_affinity(...)` | MHCflurry / NetMHCpan presentation and affinity scoring via `mhctools` / `topiary` |
| `resolve_species_context(...)` / `allele_prediction_capability(...)` | Species-aware Ensembl release defaults and predictor capability checks |

## Python quick start

### Patient target prioritization (with sex eligibility and a fusion breakpoint)

```python
from tsarina import personalized_targets

targets = personalized_targets(
    hla_alleles=["HLA-A*02:01", "HLA-B*07:02"],
    cta_expression={"PRAME": 87.3, "MAGEA4": 142.5},
    mutations=["KRAS G12D"],
    fusions=[
        {
            "label": "EWSR1-FLI1",
            "gene_5prime": "EWSR1",
            "gene_3prime": "FLI1",
            "n_terminal_seq": "SYGQQSSY",
            "c_terminal_seq": "NPSYDSVR",
        }
    ],
    patient_sex="female",
    sex_eligibility="all_adult",
)
```

### Explore built-in catalogs without running predictors

```python
from tsarina import (
    CTA_selection,
    CTA_sex_restricted_antigens,
    cancer_fusions,
    get_panel,
)

human_ctas = CTA_selection(category="CTA")
mouse_ctas = CTA_selection(category="CTA", species="mouse")
female_eligible = CTA_sex_restricted_antigens(sex="female", eligibility="all_adult")
ewing_fusions = cancer_fusions(cancer_code="SARC_EWS", defining_only=True)
alleles = get_panel("global54_abc")
```

### Population CTA × HLA panel

```python
from tsarina import spanning_pmhc_set

panel = spanning_pmhc_set(
    ctas=["PRAME", "MAGEA4", "CTAG1B"],
    alleles=["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02"],
    output_format="long",
)
```
