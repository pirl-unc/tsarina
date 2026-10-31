# CTA ownership and downstream evidence

Tsarina consumes one executable cancer-testis antigen definition:
[oncoref](https://github.com/pirl-unc/oncoref). Tsarina enriches that definition
with target-selection evidence, but never maintains a second CTA universe.

## Current selection and historical evidence

Tsarina 1.36 requires OncoRef >=1.23.1. Default human membership is the
`CTA` category from `oncoref.cta_selection.results(species="human")`.
This protein-level policy owns entity identities, member genes, categories,
assessment status, evidence flags and policy/reference/source receipts.
Tsarina does not reproduce its biological thresholds.

```python
from tsarina import CTA_selection, CTA_gene_ids, CTA_evidence, CTA_hpa_evidence

selected = CTA_selection(category="CTA")  # Unmodified upstream protein decisions.
assert CTA_gene_ids() == {g for members in selected.gene_ids for g in members}
evidence = CTA_evidence()                # One row per current member gene.
historical = CTA_hpa_evidence()          # Original oncoref.cta HPA table.
```

With OncoRef 1.23.1, the current selection contains 522 protein groups and 558
gene loci, compared with 624 historical HPA default genes: 225 genes enter and
291 leave. CTAG2 is admitted by the current policy. GAGE10 remains unassigned.
These counts describe this release, not thresholds enforced by Tsarina.

`CTA_evidence()` starts from **all current decisions**, projects protein
entities onto member genes, and left-joins historical HPA annotations and
Tsarina's MS overlay by Ensembl gene ID. The `selection_*` columns preserve
the upstream decision; `is_current_cta` is true only for category `CTA`.
Current member identities supply symbols, with stable gene IDs used when a
symbol is missing. Missing HPA measurements stay missing, with empty text or
`NO_DATA` for unavailable restriction axes; they do not become zero or pass
an HPA filter. There are 92 selected genes absent from the historical table
in OncoRef 1.23.1.

`Unsafe`, `Other`, and null/unassigned decisions remain distinct.
`CTA_excluded_gene_ids/names()` returns only `Unsafe` members;
`CTA_unfiltered_gene_ids/names()` returns the complete current decision
universe. `CTA_filtered_gene_ids/names()` is a compatibility name for current
CTA membership. Protein groups come from the current selected entities.

Historical testis/placental/relaxed-reproductive, never-expressed and clinical
helpers remain explicit `oncoref.cta` compatibility views. Historical flags
such as `passes_filters`, `never_expressed` and `specificity_action` annotate
evidence; they do not decide current membership. Clinical helpers can overlap
current CTAs. The vaccine `strict`, `loose` and `both` modes remain available
for historical comparisons; the human default is `current`.

## Sequence background and provenance

For gene-based peptide consumers, current selected member genes form the CTA
set. Every other coding gene in the requested Ensembl sequence release remains
in the non-CTA background, including `Unsafe`, `Other`, unassigned and genes
absent from the selection table. The compatibility partition field
`cta_never_expressed` is empty. Other or missing assessment does not justify
removing a sequence from the background.

This consumer policy is recorded as
`tsarina.current-selection-other-genes.v1`. It is separate from OncoRef's
independently versioned exact-protein background policy. Admission of a
representative protein does not certify every isoform of its member genes.

Outputs retain selection policy/reference/source IDs separately from the
consumer's sequence release. OncoRef 1.23.1's selection reference describes
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
| Cancer RNA prevalence | Derived HPA cancer-prevalence tables | No |
| Cancer IHC prevalence | Derived HPA cancer-prevalence tables | No |
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

The HPA cancer-prevalence tables are feature caches for downstream ranking.
They are restricted to oncoref gene IDs and contain no CTA membership or
specificity fields.

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
| SUN5 / SPAG4L | Testis-restricted RNA and enhanced testis IHC | In the default CTA set |
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
| Import `cancerdata` or `oncodata` for CTA definitions | Import `oncoref.cta_selection` for current membership |
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
data is also present in Tsarina's minimum supported oncoref 1.23.1.

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

Tsarina's obsolete CTA regeneration and full reference-data commands have
been removed. The remaining HPA cancer-prevalence features can be regenerated
with `scripts/regenerate_hpa_cancer_prevalence.py`; that script takes its historical HPA gene
universe from `oncoref.cta.cta_evidence()`; new current members can therefore
lack these cached ranking features.
