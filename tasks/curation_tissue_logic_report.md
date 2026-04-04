# Curation And Tissue Logic Audit Report

Date: 2026-04-04

Scope:
- `tsarina` CTA curation, restriction synthesis, and gene-set helpers
- `hitlist` source classification, tissue categories, and PMID override handling

Method:
- Read the relevant code and docs.
- Queried the shipped `CTA_evidence()` table for internal inconsistencies.
- Ran small synthetic calls against `hitlist.curation.classify_ms_row()` to verify override behavior.

## Findings

### 1. High: `tsarina`'s local IEDB wrapper drops `hitlist` curation inputs, so curated adjacent/healthy calls are wrong

Code path:
- `tsarina/iedb.py:301-304`
- `tsarina/iedb.py:381-382`
- `hitlist/hitlist/curation.py:218-245`
- `hitlist/hitlist/data/pmid_overrides.yaml:37-70`

What is happening:
- `tsarina.iedb.scan_public_ms()` and `tsarina.iedb.profile_dataset()` call `hitlist.curation.classify_ms_row(...)`, but they do not pass `pmid` or `mhc_restriction`.
- `hitlist` needs `pmid` to apply study-specific overrides and `mhc_restriction` for mono-allelic detection.
- As a result, `tsarina` does not actually use the full `hitlist` curation logic even though the package says MS curation is delegated to `hitlist`.

Concrete example:
- `classify_ms_row("No immunization", "healthy", "Direct Ex Vivo", "Lung", pmid=35051231)` returns `src_adjacent_to_tumor=True`.
- The same call without the PMID returns `src_healthy_tissue=True`.
- `tsarina` currently behaves like the second case, because it omits the PMID.

Impact:
- `target_peptides(..., cancer_specific=True)` can treat curated tumor-adjacent samples as healthy tissue or vice versa.
- `healthy_tissue_peptides()` and `healthy_tissue_peptide_details()` can build the wrong negative set.
- The behavior of `tsarina` and `hitlist` diverges for the same underlying IEDB rows.

Recommendation:
- Stop reimplementing scanning in `tsarina/iedb.py`; call `hitlist.scanner.scan()` directly.
- If the wrapper stays, pass through `pmid`, `submission_id`, `mhc_restriction`, and host fields exactly as `hitlist.scanner` does.

### 2. High: `hitlist` defines a `healthy` override type, but the classifier never applies it

Code path:
- `hitlist/hitlist/curation.py:247-276`
- `hitlist/hitlist/data/pmid_overrides.yaml:10-16`
- `hitlist/hitlist/data/pmid_overrides.yaml:26-40`
- `hitlist/hitlist/data/pmid_overrides.yaml:83-100`
- `hitlist/docs/pmid-curation.md:7-13`

What is happening:
- YAML and docs define `healthy` as a valid override.
- `classify_ms_row()` has branches for `cancer_patient`, `adjacent`, `activated_apc`, and `cell_line`, but not for `healthy`.
- Any `healthy` override therefore falls through to the default structured-field logic.

Concrete example:
- `classify_ms_row("Occurrence of cancer", "melanoma", "Direct Ex Vivo", "Liver", pmid=33858848)` still returns `src_cancer=True` even though PMID `33858848` is declared `override: healthy`.
- The same failure happens for rule-level `healthy` overrides such as PMID `36589698`.

Impact:
- Any study that actually needs a `healthy` override to correct bad IEDB fields will still be misclassified.
- The YAML gives a false sense that these studies are protected by curation when they are not.

Recommendation:
- Add an explicit `healthy` branch that forces the row into the healthy-donor path.
- Add a regression test using a deliberately misannotated row plus a `healthy` PMID override.

### 3. Medium: `hitlist` claims submission-ID override support, but the code does not implement it

Code path:
- `hitlist/hitlist/data/pmid_overrides.yaml:21-22`
- `hitlist/hitlist/scanner.py:239`
- `hitlist/hitlist/scanner.py:256-266`
- `hitlist/hitlist/curation.py:150-158`
- `hitlist/hitlist/curation.py:196-244`

What is happening:
- The override YAML says studies without PMIDs can use `submission_id`, and that the curation code checks both.
- `scanner.py` reads `submission_id`, but `classify_ms_row()` has no `submission_id` parameter and only looks up integer PMIDs.

Impact:
- Any dataset rows keyed only by submission ID cannot be curated despite the docs and YAML comment saying they can.
- This is currently a silent capability gap.

Recommendation:
- Either implement `submission_id` lookup end-to-end or remove the claim from the YAML/docs.

### 4. High: `tsarina` axis-based CTA helpers leak filtered-out genes back into the public sets

Code path:
- `tsarina/gene_sets.py:139-147`
- `tsarina/gene_sets.py:150-223`

What is happening:
- `_axis_filter()` and `CTA_by_axes()` operate on the full 358-row CTA universe.
- They do not restrict to `filtered == True`.
- The helper names and docs read like they return CTA subsets, but they currently return filtered plus excluded candidates.

Concrete evidence from the shipped table:
- `CTA_testis_restricted_gene_names()` returns 245 genes, including 18 filtered-out genes:
  `ACRBP, CABYR, CTAG2, EPPIN, GAGE2E, GAPDHS, GOLGA6B, GOLGA6C, GPX5, MAEL, MAGEA11, MAGEB3, PATE1, PATE4, ROPN1B, SKA3, TEX15, XAGE1B`
- `CTA_by_axes(restriction="REPRODUCTIVE")` includes 4 filtered-out genes:
  `BSPH1, CSN3, LALBA, SPAG11A`

Impact:
- Downstream users can accidentally resurrect genes that were explicitly excluded from the curated CTA set.
- README counts are now misleading because the helper outputs no longer match the curated population they describe.

Recommendation:
- Make axis helpers filtered-by-default.
- If access to the unfiltered universe is needed, add an explicit `filtered_only=False` knob instead of the current silent behavior.

### 5. High: the default CTA set is described as “reproductive-restricted”, but 16 shipped default genes are synthesized as `SOMATIC`

Code path:
- `README.md:51-54`
- `README.md:101-113`
- `tsarina/tiers.py:158-214`

What is happening:
- `CTA_gene_names()` returns 257 “expressed CTAs”.
- The README describes these as “257 expressed, reproductive-restricted CTAs”.
- In the shipped `CTA_evidence()` table, 16 genes in that default set have `restriction == "SOMATIC"` and no protein evidence, meaning the synthesized restriction logic itself says they are somatic rather than reproductive-only.

Examples from the shipped table:

| Symbol | Max somatic tissue | Max somatic nTPM | Somatic tissues detected |
|---|---|---:|---:|
| `RNF148` | cerebellum | 9.9 | 2 |
| `TULP2` | cerebral cortex | 8.4 | 15 |
| `SPATA22` | white matter | 4.4 | 9 |
| `XAGE1A` | liver | 1.5 | 1 |
| `MEIOB` | kidney | 1.0 | 1 |

Interpretation:
- This comes from the current “deflated fraction >= 0.99” rule for no-protein genes.
- A gene can therefore pass the CTA filter even when `rna_restriction` is `SOMATIC`, as long as testis-heavy expression keeps the fraction above 0.99.

Impact:
- The public description of the default CTA set overstates how clean the restriction filter currently is.
- Users relying on the helper name alone would assume these genes have no meaningful somatic RNA evidence, which is false.

Recommendation:
- Decide which claim is correct:
  - If the biological intent is truly “reproductive-restricted”, require `rna_reproductive == True` when protein data is missing/unreliable.
  - If the intent is “reproductive-dominant but not necessarily exclusive”, rename the filter/set/docs accordingly.

### 6. Low: public annotations/docs are stale relative to the shipped restriction schema

Code path:
- `README.md:64-71`
- `README.md:167-179`
- `tsarina/evidence.py:74-101`
- `tsarina/__init__.py:53-86`

What is happening:
- Docs still advertise `OVARIAN` restriction values and old source flags like `src_healthy`, `src_reproductive`, and `src_thymus`.
- The shipped CSV actually uses:
  - `TESTIS`, `PLACENTAL`, `REPRODUCTIVE`, `SOMATIC`, `NO_DATA`
  - `src_healthy_tissue`, `src_healthy_reproductive`, `src_healthy_thymus`
- `tsarina.__all__` still exports `CTA_ovarian_restricted_gene_ids` and `CTA_ovarian_restricted_gene_names`, but those symbols do not exist. `from tsarina import *` raises `AttributeError`.

Impact:
- Public-facing documentation and exports no longer describe the actual annotation schema.
- This is easy to trip over during exploratory analysis.

Recommendation:
- Update the docs/docstrings to match the real columns and values.
- Remove or implement the stale ovarian exports.

## Overall Assessment

The biggest practical risk is not the tissue vocabulary itself; it is curation drift between the packages:
- `hitlist` has override machinery that is partially nonfunctional (`healthy`, `submission_id`).
- `tsarina` bypasses part of that machinery anyway by using its own scanner wrapper.

On the CTA side, the most important issue is semantic drift:
- the default/public CTA sets are described as reproductive-restricted,
- but the shipped table and exported helpers do not consistently enforce that meaning.

If I were fixing this next, I would do it in this order:
1. Make `tsarina` use `hitlist.scanner.scan()` directly.
2. Fix `healthy` override handling in `hitlist`.
3. Decide whether CTA filtering is meant to be exclusive or merely dominant, then align the data, helpers, and docs to that decision.
