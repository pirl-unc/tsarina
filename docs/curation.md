# CTA ownership and downstream evidence

Tsarina consumes one executable cancer-testis antigen definition:
[oncoref](https://github.com/pirl-unc/oncoref). Tsarina enriches that definition
with target-selection evidence, but never maintains a second CTA universe.

## Current selection and historical evidence

Tsarina 1.59 requires Python >=3.11, OncoRef >=1.55.0 and Hitlist >=2.15.0.
Default human membership is the `CTA` category from OncoRef's rebuilt
set-based target catalog (`oncoref.target_selection.results()`), enriched with
canonical gene identities from `oncoref.gene_identity_catalog()`. This
protein-level policy owns entity identities, member genes, categories,
`cta_tier` classifications, 4-way sex-restriction codes (`CTA-M`, `CTA-F`,
`CTA-MF`, `CTA-P`, `RA-*-only`), nuclear snRNA-seq caution, clinical
organ/epitope clearance, 1-mm near-self coverage, assessment status, evidence
flags and policy/reference/source receipts. Tsarina does not reproduce its
biological thresholds.

```python
from tsarina import (
    CTA_by_tier,
    CTA_cancer_antigen_gene_ids,
    CTA_evidence,
    CTA_gene_ids,
    CTA_hpa_evidence,
    CTA_reproductive_antigen_gene_ids,
    CTA_selection,
    CTA_target_candidate_gene_ids,
)

selected = CTA_selection(category="CTA")  # Normalized OncoRef target_selection decisions.
assert CTA_gene_ids() == {g for members in selected.gene_ids for g in members}
core_ctas = CTA_by_tier("core")           # core / core_with_extended_exposure / extended_only
evidence = CTA_evidence()                 # One row per current member gene.
historical = CTA_hpa_evidence()           # Original oncoref.cta HPA table.
```

With OncoRef 1.55.0 (policy `oncoref.target-selection.2026-11-26.2`), the
rebuilt human catalog classifies 19,893 reference protein groups (20,099 member
genes) into **322 CTA** (354 gene loci: 266 `core` / 296 genes, 33
`core_with_extended_exposure` / 34 genes, 23 `extended_only` / 24 genes;
214 `ectopic`, 97 `ectopic_and_lineage_retained`, 11 `lineage_retained`;
248 `replicated`, 74 `single_sample`),
**66 Cancer Antigen only** (426 total CA genes), **146 Reproductive Antigen
only** (501 total RA genes; 573 `CA|RA` target candidate genes), **1,110 Other**
(`Other:Review`, `Other:Expression`, `Other:Unscreened tissue`,
`Other:Sequence`) and **18,249 Unsafe** (18,406 genes). Key adjudications in
`2026-11-26.2` include `CTAG2` admitted as `CTA` (`core`, `CTA-MF`) via
epitope-scoped clinical organ clearance (`SLLMWITQC` on LAGE-1a
`ENST00000369585`), `PLAC1` and `DUX4` admitted as `CTA` (`core`, with
`rare_cell_status="nuclear_only"` and `nuclear_rna_caution`), `NUTM1` and
`PRM1` admitted as `CTA`, `LALBA` admitted as `CTA` (`extended_only`, `CTA-F`),
`PRM2` and chemosensory ORs (`OR1N1`, `OR4M2`, `OR10J1`) in `Reproductive Antigen
only`, `GAGE10` held in `Other:Review`, and `SP7`, `MAGEA3`, `MAGEA12`, `SUN3`
in `Unsafe`.
Cross-species catalogs are also genome-wide evaluated for **mouse** (8 CTAs,
311 sex-restricted antigens, Ensembl 115) and **dog** (118 CTAs across 19
cohorts / 935 primary tumors including `PRAME` and `SSX6`, 431 sex-restricted
antigens, Ensembl 98 CanFam3.1), both with `representative_transcripts`,
`tumor_activation_mode`, and `tumor_support`.
OncoRef's `target_selection.metadata()["release_ready"]` is currently `False`
while three listed scientific reviews remain open (tumor-source comparability,
read-level rare-cell mapping-risk confirmation, and hair-follicle IHC replay);
unresolved evidence remains excluded from positive sets and Tsarina preserves
`cta_selection_release_ready` and `cta_selection_open_reviews` in selection
provenance. These counts describe this release, not thresholds enforced by
Tsarina.

`CTA_evidence()` starts from **all current decisions**, projects protein
entities onto member genes, and left-joins historical HPA annotations and
Tsarina's MS overlay by Ensembl gene ID. The `selection_*` columns preserve
the upstream decision (including `selection_cta_tier`,
`selection_classification`, `selection_is_cancer_antigen`,
`selection_is_reproductive_antigen`, `selection_is_target_candidate`,
`selection_sex_restriction_class`, `selection_sex_restriction_code`,
`selection_tumor_activation_mode`, `selection_tumor_support`,
`selection_qualifying_cancer_codes`,
`selection_rare_cell_status`, `selection_nuclear_rna_caution`,
`selection_clinical_clearance_scope`, `selection_clinical_cleared_organs`,
`selection_clinical_cleared_epitopes`, `selection_clinical_hla_alleles`,
`selection_clinical_cleared_transcripts`, `selection_clinical_unverified_isoforms`,
`selection_exact_unique_9mers`, `selection_safe_2mm_9mers`, and
`selection_unsafe_near_self_1mm_pct`);
`is_current_cta` is true only for category `CTA`. Current member identities
supply symbols, with stable gene IDs used when a symbol is missing. Missing
HPA measurements stay missing, with empty text or `NO_DATA` for unavailable
restriction axes; they do not become zero or pass an HPA filter. There are
41 selected CTA genes absent from the historical HPA table in OncoRef 1.55.0.

`Unsafe`, `Other`, `Cancer Antigen only`, and `Reproductive Antigen only`
decisions remain distinct. The rebuilt policy gives every protein a category
and reports incomplete coverage in `selection_assessment_status`; a null
category from an older release is still handled as unassigned, never as
`Unsafe` or `CTA`.
`CTA_excluded_gene_ids/names()` returns only `Unsafe` members;
`CTA_unfiltered_gene_ids/names()` returns the complete current decision
universe. `CTA_filtered_gene_ids/names()` is a compatibility name for current
CTA membership. Protein groups come from the current selected entities.
Non-human `CTA_selection(species=...)` queries delegate to
`oncoref.cta_selection.results()`.

Historical testis/placental/relaxed-reproductive, never-expressed and clinical
helpers remain explicit `oncoref.cta` compatibility views. Historical flags
such as `passes_filters`, `never_expressed` and `specificity_action` annotate
evidence; they do not decide current membership. Clinical helpers can overlap
current CTAs. The vaccine `strict`, `loose` and `both` modes remain available
for historical comparisons; the human default is `current`.

## Sequence background and provenance

For gene-based peptide consumers, current selected CTA member genes form the
CTA set. Audited same-chromosome `HSCHR*` alternate-haplotype annotations of
those curated CTAs (resolved via OncoRef's `alternate_ensembl_gene_ids` or
exact non-empty symbol on the matching chromosome) are excluded from the
non-CTA background and recorded separately (`CTA_partition_alternate_haplotype_aliases()`
and `background_cta_aliases.csv`) so an alternate-contig copy of a CTA does not
self-veto its primary gene. Primary-contig non-CTA loci never inherit CTA
membership. Every other coding gene in the requested Ensembl sequence release
remains in the non-CTA background, including `Unsafe`, `Other`,
`Cancer Antigen only`, `Reproductive Antigen only`, unassigned and genes
absent from the selection table. The compatibility partition field
`cta_never_expressed` is empty. Other or missing assessment does not justify
removing a sequence from the background.

This consumer policy is recorded as
`tsarina.current-selection-alternate-haplotypes.v1`. It is separate from OncoRef's
independently versioned exact-protein background policy. Admission of a
representative protein does not certify every isoform of its member genes.

Outputs retain selection policy/reference/source IDs separately from the
consumer's sequence release. OncoRef's selection reference describes
its Ensembl-115 gene space; it is not the sequence asset used by Tsarina's
default Ensembl-112 enumeration. Vaccine provenance also records the actual
protein reference digest. See the
[upstream selection contract](https://pirl-unc.github.io/oncoref/cta-selection/).

## Downstream evidence flow

After CTA membership is known, Tsarina adds evidence used for selection and
ranking:

| Evidence | Storage or source | Can change CTA membership? |
|---|---|---|
| Gene-level MS safety | `data/gene-ms-safety-evidence.csv` | No |
| Live immunopeptidomics | IEDB/CEDAR through hitlist | No |
| Cancer RNA prevalence | `oncoref.hpa_cancer_rna_prevalence()` (HPA v25.1) | No |
| Cancer IHC prevalence | `oncoref.hpa_cancer_ihc_prevalence()` (HPA v25.1) | No |
| pMHC prediction and scoring | Runtime | No |

The static MS overlay has only this schema:

```text
Ensembl_Gene_ID
ms_restriction
ms_healthy_somatic_tissues
ms_pmids
```

It contains no symbol, membership flag, specificity action, HPA restriction
field, or proteoform definition. Tsarina left-joins these columns onto the
oncoref evidence frame. Therefore, an overlay row for a non-CTA gene cannot
enter the CTA universe.

HPA cancer RNA and IHC prevalence tables are delegated directly to
`oncoref.hpa_cancer_rna_prevalence()` and `oncoref.hpa_cancer_ihc_prevalence()`
(HPA v25.1, genome-wide, covering 353/357 current CTAs for RNA and 195/357 for
IHC) rather than bundled as stale CSVs in Tsarina. They supply downstream
ranking features only and never change CTA membership.

## Restriction axes and API behavior

The historical columns below come from `oncoref.cta` where HPA rows exist:

| Axis | Values or role |
|---|---|
| `protein_restriction` | HPA protein-IHC restriction call |
| `rna_restriction` | HPA RNA restriction call |
| `rna_restriction_level` | `STRICT`, `MODERATE`, or `PERMISSIVE` |
| `restriction` | Synthesized HPA restriction |
| `restriction_confidence` | `HIGH`, `MODERATE`, `LOW`, or `NO_DATA` |

Tsarina's `ms_restriction` is an independent downstream safety axis:

- `CANCER_ONLY`
- `EXPECTED_TISSUE`
- `SINGLETON_HEALTHY`
- `RECURRENT_HEALTHY`
- `UNCLASSIFIED_MS`
- `NO_MS_DATA`

`CTA_by_axes()` applies explicit HPA/MS predicates to the current selected
members. With `filtered_only=False`, it queries all current decisions. Default
panel and personalized selection do not reapply historical HPA confidence,
vital-RNA or MAGE-family exclusions; these remain optional user filters:

```python
from tsarina import CTA_by_axes

strict_testis = CTA_by_axes(
    restriction="TESTIS",
    rna_restriction_level="STRICT",
)
high_confidence = CTA_by_axes(
    restriction="TESTIS",
    restriction_confidence="HIGH",
)
```

Live workflows can explicitly call `synthesize_restriction()` or
`assign_all_axes()` to combine current IEDB/CEDAR MS evidence with the oncoref
HPA axes. These helpers do not derive or modify the upstream protein and RNA
classifications. The combined tissue scope comes directly from oncoref's union
of the observed modalities: testis protein with somatic RNA remains `SOMATIC`,
and testis protein with placental RNA becomes `REPRODUCTIVE`. Only matching
protein/RNA labels receive an agreement bonus. Without MS evidence, both the
scope and confidence are the upstream result; MS evidence contributes only to
the explicitly requested confidence calculation.

## Coverage limits: endogenous retrovirus antigens

**The gene-based CTA panel does not comprehensively cover HERV-K (HML-2)
Env, Rec, or Np9 antigens.** Their expression can arise from multiple proviral
loci and transcript isoforms. Searching HGNC symbols for `ERVK` or evaluating
one gene such as `ERVK3-1` cannot establish family-wide expression or safety.
Gene-level ERV-derived entries, when present in oncoref, describe those genes
only. An absent gene-panel hit is not evidence that the tumor lacks HERV-K
antigens. This is the scope limitation tracked in
[#120](https://github.com/pirl-unc/tsarina/issues/120).

Locus/family expression belongs in a separate antigen-evidence source.
[Telescope](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006453)
addresses ambiguous read assignments to estimate expression at individual
transposable-element loci.
[ERVmap](https://pmc.ncbi.nlm.nih.gov/articles/PMC6294949/)
uses a curated proviral-locus reference and stringent mapping filters. Their
annotations and read-assignment policies differ; resulting counts are not
interchangeable with gene TPM or with each other. Neither method by itself
establishes a translated antigen or HLA presentation.

Tsarina currently has no HERV locus-expression input or adapter for these
outputs. Its `--cta` input selects gene symbols, and its supported `--virus`
proteomes do not provide a HERV-K route. A future integration needs:

- genome assembly, annotation/quantifier versions, locus coordinates or
  explicitly marked family-level identifiers, and read-assignment uncertainty;
- matched tumor and normal-tissue expression with sample provenance;
- locus/transcript/ORF-to-peptide mapping, retaining ambiguous peptide sources;
- separate protein, immunopeptidomic, and HLA-assignment evidence before
  interpretation as a presented target.

HERV-K is a meaningful experimental target class: for example,
[HERV-K Env-directed CAR T-cell work in melanoma](https://pmc.ncbi.nlm.nih.gov/articles/PMC4506228/)
studied recognition of the envelope protein. Such evidence does not validate
every HML-2 locus, demonstrate presentation of a particular peptide-HLA
complex, or establish normal-tissue safety. Endogenous human sequences also
need an explicit specificity assessment; the existing exogenous-virus
human-exclusivity filter is not a substitute for ERV-specific curation.

## Reviewed SUN-domain candidates

SUN-domain membership alone does not establish a testis-restricted CTA. The
[#131](https://github.com/pirl-unc/tsarina/issues/131) audit confirms that
oncoref already supplies both SUN5 and SUN3 to Tsarina's evidence table:

| Gene | HPA v23 RNA/protein evidence | Current handling |
|---|---|---|
| SUN5 / SPAG4L | Testis-restricted RNA and enhanced testis IHC | `Reproductive Antigen only` under `target_selection` (reproductive-restricted without qualifying cancer enrichment; not `Unsafe`) |
| SUN3 / SUNC1 | Testis-restricted RNA, but IHC in bronchus, kidney, skin, and testis | Retained in evidence; excluded from default selection |
| SPAG4 / SUN4 / CT127 | Pancreas 87.3 nTPM versus testis 31.1 nTPM; additional somatic RNA | Outside default selection; candidate-reference omission tracked in [oncoref #548](https://github.com/pirl-unc/oncoref/issues/548) |
| SUN1 / SUN2 | Broad somatic expression | Outside default selection |

The original human [SPAG4 study](https://pubmed.ncbi.nlm.nih.gov/14614621/)
reports pancreas as well as testis expression. Its
[renal-cell carcinoma evidence](https://pubmed.ncbi.nlm.nih.gov/23602831/)
does not remove that normal-tissue caveat. Mouse spermatid-specific SUN4
findings should not override the human tissue data.

[SUN5 colorectal-cancer experiments](https://pmc.ncbi.nlm.nih.gov/articles/PMC9654567/)
support tumor-associated expression and a role in proliferation/migration;
they do not establish peptide presentation or clinical safety. SUN3's
discordant IHC remains visible rather than being presumed cross-reactivity.
Candidate additions or evidence reinterpretation belong in oncoref.

## Updating definitions

### Migrating older installations

The reference package was renamed from `cancerdata` to `oncodata` and then
to `oncoref`. Use `oncoref` for new code and environments; installing the
intermediate name does not update Tsarina's current reference data.

| Older workflow | Current workflow |
|---|---|
| Import `cancerdata` or `oncodata` for CTA definitions | Import `oncoref.target_selection` (or `oncoref.antigen_sets`) for current membership |
| Run `scripts/sync_proteoform_groups.py` to refresh a local mirror | Read current selected entity/member identities directly |
| Skip registry integration checks when the source package is absent | Install the required `oncoref` dependency and run the integration tests |

The sync script and mirrored registry have been removed. No manual sync is
needed after upgrading oncoref. The required integration check in
`tests/test_spanning.py` compares Tsarina's groups with the live oncoref API;
it does not use `pytest.importorskip`.

For development, run `./develop.sh` with the intended virtualenv active. It
installs an adjacent oncoref checkout in editable mode and reports its resolved
path. Verify that both imports point to the intended checkouts before comparing
curation results:

```bash
python -c 'import oncoref, tsarina; print(oncoref.__file__); print(tsarina.__file__)'
pytest tests/test_oncoref_authority.py tests/test_gene_sets.py tests/test_spanning.py
```

### Changing the canonical definitions

Protein evidence must be refreshed along with RNA when adding or reviewing a
candidate. The historical RNA-only additions in
[#130](https://github.com/pirl-unc/tsarina/issues/130) are now corrected in the
oncoref data consumed by Tsarina: CGB2/3/5/7, PSG4/6/7, CT45A5, CSH1, and
GAGE10 all retain their HPA v23 IHC evidence. Recomputing their protein columns
from the pinned `normal_tissue` table reproduces the current values; regression
tests cover all ten genes, including CGB2's expressed status. The corrected
data is also present in Tsarina's minimum supported oncoref 1.43.0.

A protein restriction call is evidence, not an admission decision. For
example, CGB3 has reproductive-tissue IHC, while oncoref still excludes it
from the default CTA set after considering the full evidence. Do not restore
old membership predictions from an issue's snapshot or regenerate a separate
Tsarina table to change that decision.

For changes to CTA membership, aliases, HPA classifications, evidence
provenance, or proteoform groups:

1. file the omission or error in
   [oncoref issues](https://github.com/pirl-unc/oncoref/issues);
2. correct and release oncoref; and
3. raise Tsarina's minimum oncoref dependency when the correction requires a
   new API or data release.

Do not add a local CTA row, exclusion list, alias override, or fallback table
to Tsarina. The ownership audit followed this rule by filing
[oncoref #435](https://github.com/pirl-unc/oncoref/issues/435) for an
RNA/protein restriction-confidence synthesis error instead of compensating
locally.

Tsarina's obsolete CTA regeneration, HPA prevalence regeneration, and full
reference-data commands have been removed; HPA cancer RNA/IHC prevalence is
loaded directly from `oncoref.hpa_cancer_rna_prevalence()` and
`oncoref.hpa_cancer_ihc_prevalence()`.
