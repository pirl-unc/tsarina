# CTA ownership and downstream evidence

Tsarina consumes one executable cancer-testis antigen definition:
[oncoref](https://github.com/pirl-unc/oncoref). Tsarina enriches that definition
with target-selection evidence, but never maintains a second CTA universe.

## At a glance

| Question | Owner |
|---|---|
| Is this gene a CTA candidate, default, filtered, excluded, or low-expression member? | oncoref |
| What is its canonical symbol or alias? | oncoref |
| What do HPA RNA and protein data imply about its tissue restriction? | oncoref |
| Which CTA proteoform group contains it? | oncoref |
| How prevalent is it across cancer samples? | Tsarina downstream features |
| Has a peptide or gene been observed by public MS in cancer or healthy tissue? | Tsarina/hitlist evidence |
| Which pMHCs score best for a patient or HLA panel? | Tsarina runtime selection |

An oncoref release defines the CTA rows. A Tsarina evidence join may add
columns to those rows, but it cannot create a CTA.

## Definition boundary

Oncoref owns:

- the CTA candidate universe and gene identities;
- default, filtered, unfiltered, low-expression, and excluded membership;
- symbols, aliases, source provenance, and specificity decisions;
- HPA RNA and protein evidence, restriction axes, and confidence; and
- CTA proteoform groups.

Tsarina does not ship a CTA definition table or reproduce those decisions.
Foundational helpers such as `CTA_gene_names`, `CTA_filtered_gene_ids`,
`CTA_excluded_gene_names`, `CTA_testis_restricted_gene_names`, and
`cta_symbol_for_alias` directly alias their `oncoref.cta` implementations.
Proteoform grouping comes from
`oncoref.proteoforms.proteoform_symbol_map(scope="cta")`.

The boundary is executable:

```python
import oncoref.cta
import tsarina.gene_sets

assert tsarina.gene_sets.CTA_gene_names is oncoref.cta.cta_gene_names
assert (
    tsarina.gene_sets.CTA_filtered_gene_ids
    is oncoref.cta.cta_filtered_gene_ids
)
```

`CTA_evidence()` starts from `oncoref.cta.cta_evidence()` and preserves its
exact row universe and upstream columns.

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

The static columns below come unchanged from oncoref:

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

`CTA_by_axes()` queries oncoref's row universe and membership tiers. Its only
Tsarina extension is the optional `ms_restriction` predicate:

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
classifications.

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
| Import `cancerdata` or `oncodata` for CTA definitions | Import `oncoref.cta` |
| Run `scripts/sync_proteoform_groups.py` to refresh a local mirror | Read `oncoref.proteoforms.proteoform_symbol_map(scope="cta")` directly |
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
data is also present in Tsarina's minimum supported oncoref 1.8.150.

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
with `scripts/regenerate_hpa_cancer_prevalence.py`; that script takes its gene
universe from `oncoref.cta.cta_evidence()`.
