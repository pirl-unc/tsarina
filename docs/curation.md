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

## Updating definitions

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
