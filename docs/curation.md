# CTA ownership and downstream evidence

## One definition library

[oncoref](https://github.com/pirl-unc/oncoref) is the sole authority for
cancer-testis antigen definitions. It owns:

- the CTA candidate universe and gene identities;
- default, filtered, unfiltered, low-expression, and excluded membership;
- symbols, aliases, source provenance, and specificity decisions;
- HPA RNA and protein evidence, restriction axes, and confidence; and
- CTA proteoform groups.

Tsarina does not ship a CTA definition table or recreate those decisions.
Foundational helpers such as `CTA_gene_names`, `CTA_filtered_gene_ids`,
`CTA_excluded_gene_names`, `CTA_testis_restricted_gene_names`, and
`cta_symbol_for_alias` are direct aliases of their `oncoref.cta`
implementations. Proteoform grouping is loaded from
`oncoref.proteoforms.proteoform_symbol_map(scope="cta")`.

This makes an installed oncoref release the single executable definition:

```python
import oncoref.cta
import tsarina.gene_sets

assert tsarina.gene_sets.CTA_gene_names is oncoref.cta.cta_gene_names
assert tsarina.gene_sets.CTA_filtered_gene_ids is oncoref.cta.cta_filtered_gene_ids
```

`CTA_evidence()` starts from `oncoref.cta.cta_evidence()`. Tsarina preserves
the same CTA row universe and every upstream column.

## Tsarina's evidence layer

Tsarina owns evidence used to select and score targets after CTA membership is
known:

| Evidence | Storage or source | Effect on CTA membership |
|---|---|---|
| Gene-level MS safety | `data/gene-ms-safety-evidence.csv` | None |
| Live immunopeptidomics | IEDB/CEDAR through hitlist | None |
| Cancer RNA prevalence | Derived HPA cancer-prevalence tables | None |
| Cancer IHC prevalence | Derived HPA cancer-prevalence tables | None |
| pMHC predictions and scoring | Runtime | None |

The static MS overlay is deliberately generic. Its schema is only:

```text
Ensembl_Gene_ID
ms_restriction
ms_healthy_somatic_tissues
ms_pmids
```

It contains no symbol, CTA membership flag, specificity action, HPA
restriction field, or proteoform definition. A left join adds these columns
only to genes already present in oncoref's evidence frame. Consequently, an MS
row for a non-CTA gene cannot introduce that gene into Tsarina's CTA universe.

The HPA cancer-prevalence tables are downstream feature caches. They are
restricted to oncoref gene IDs and do not contain CTA membership or
specificity columns.

## Restriction axes

The static `protein_restriction`, `rna_restriction`,
`rna_restriction_level`, `restriction`, and `restriction_confidence` columns
come unchanged from oncoref.

Tsarina's `ms_restriction` is an independent downstream safety axis:

- `CANCER_ONLY`
- `EXPECTED_TISSUE`
- `SINGLETON_HEALTHY`
- `RECURRENT_HEALTHY`
- `UNCLASSIFIED_MS`
- `NO_MS_DATA`

`CTA_by_axes()` queries oncoref's row universe and membership tiers. Its only
Tsarina extension is an optional `ms_restriction` predicate.

Live target-selection workflows can explicitly call Tsarina's
`synthesize_restriction()` or `assign_all_axes()` to combine current
IEDB/CEDAR MS evidence with the oncoref HPA axes. Those helpers do not derive
or modify the upstream protein and RNA classifications.

## Updating definitions

To update CTA membership, evidence provenance, aliases, HPA classifications,
or proteoform groups:

1. make or request the correction in
   [oncoref](https://github.com/pirl-unc/oncoref/issues);
2. release oncoref; and
3. raise Tsarina's minimum oncoref dependency if the change requires a new
   API or data release.

Do not add a local CTA row, exclusion list, alias override, or fallback table
in Tsarina. Upstream omissions and errors must be filed as oncoref issues
instead. For example, the audit that established this ownership boundary
reported [oncoref #435](https://github.com/pirl-unc/oncoref/issues/435) for an
RNA/protein restriction-confidence synthesis error.

Tsarina's old local CTA regeneration and reference-data commands have been
removed. HPA cancer-prevalence features can still be regenerated with
`scripts/regenerate_hpa_cancer_prevalence.py`; the script takes its gene
universe directly from `oncoref.cta.cta_evidence()`.
