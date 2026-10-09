# Canine cancer-antigen evidence and exploratory design

Cohort: **osteosarcoma**. Status: **exploratory_unassessed**.

RNA cohort: Synthetic independent-dog fixture; not canine prevalence evidence. Population: three invented dogs; repeated biopsies plus unresolved donor.

Rank exact protein sequences by absolute RNA expression in independent untreated primary-tumor dogs, then retain contiguous regions with exact positive-MS observations and qualifying DLA affinity evidence. The construct is unassessed for canine presentation, processing, safety and efficacy.

Expression threshold: 5 TPM. Repeated biopsies count once: all measured biopsies must pass for the lower bound; any possibly positive or missing biopsy contributes to the upper bound. Missing RNA is unknown. These are evidence bounds, not confidence intervals or clinical coverage.

Allowed normal tissues: testis, placenta. Policy: synthetic reproductive-restriction policy (v1); definition: strict. Somatic RNA counterevidence threshold: 1 TPM. These labels refer to the supplied canine policy; the human HPA definitions are not used.

All rejected and unknown translated source occurrences remain specificity background. Verified healthy-primary MS outside the allowed tissues supplies an 8-mer exclusion. Tumor/adjacent tissue, cell lines and unverified tissue labels do not establish healthy-normal counterevidence; absence of observations does not establish safety.

MS evidence belongs to the exact observed peptide. Monoallelic restriction is separate from affinity-inferred multiallelic/untyped assignment. Human-host DLA transfection establishes a distinct experimental context, not endogenous canine presentation. Raw sample and host/source taxonomy facets remain in the tables.

Genotype frame: synthetic cohort. Three invented matched genotypes; not a breed or global population estimate. Joint coverage uses per-target expression and MS-qualified alleles on full supplied genotypes. Pairing modes: observed. Unsupported genotype mass: 0.533; missing mass: 0.200. Unsupported genotype mass counts dogs with any unassessed allele; it can overlap demonstrated reach through another allele and is not an additional coverage category. These estimates do not imply worldwide or breed coverage. No HWE, human mortality prior or human p95 is used.

Affinity provider: frozen_affinity_archive. Capability and model/sequence hashes: `dla_capabilities.csv`. Cleavage is unassessed unless an explicit callback was supplied; injected predictions retain that provenance and do not validate a model. Coding policy: `generic`; generic coding uses a deterministic standard genetic-code reverse translation without species optimization. UTR and polyA settings and total AA/nt limits are recorded in the manifest.

## Counts

- source occurrences: 4
- exact protein groups: 3
- identified cohort dogs: 3
- cohort samples: 5
- assembled proteins: 2
- assembled native pieces: 2
- construct length aa: 18
- construct length nt: 63
- MS observed peptides: 3
- observed monoallelic peptides: 1
- inferred only peptides: 2
- MS observations: 3
- peptide DLA pairs: 3
- panel alleles: 3
- admitted groups: 2
- rejected groups: 1
- unknown groups: 0
- rejected MS observations: 3
- selection enabled alleles: 2

## Data sources

- synthetic: https://github.com/pirl-unc/tsarina/blob/main/tests/fixtures/canine-vaccine.json — version synthetic-offline-v1; asset SHA256 `9e12300be0dc926273fe72a24b00881e245b837f9a71f466fe09d1bff8f5e186`; license Apache-2.0.

Input file SHA256: `480dd82ebcf52a04ab550f8a42b2aeab2964482afd01e1e43dc53f88a1d79cad`. Reference identity: `3c57c3e384f0568d67ad1ce55826a3a4cbb12e6ef89d9faae1991b8a8d536f23`. The input copy has the same canonical JSON identity; its formatting may differ from the original file.

## Tables and figures

- [ranking.csv](ranking.csv)
- [ligands.csv](ligands.csv)
- [funnel.csv](funnel.csv)
- [occurrences.csv](occurrences.csv)
- [samples.csv](samples.csv)
- [rna_bounds.csv](rna_bounds.csv)
- [donor_prevalence.csv](donor_prevalence.csv)
- [cohort_sample_audit.csv](cohort_sample_audit.csv)
- [ms_raw.csv](ms_raw.csv)
- [ms_modality_rejected.csv](ms_modality_rejected.csv)
- [ms_support_decisions.csv](ms_support_decisions.csv)
- [normal_ms_exclusions.csv](normal_ms_exclusions.csv)
- [dla_capabilities.csv](dla_capabilities.csv)
- [capabilities_raw.csv](capabilities_raw.csv)
- [tumor_genotype_pairs.csv](tumor_genotype_pairs.csv)
- [paired_coverage.csv](paired_coverage.csv)
- [cumulative_coverage.csv](cumulative_coverage.csv)
- [affinity_predictions_used.csv](affinity_predictions_used.csv)
- [native_specific_intervals.csv](native_specific_intervals.csv)
- [protein_ms_map.csv](protein_ms_map.csv)
- [coverage-and-MS.svg](coverage-and-MS.svg)
- [dla-evidence.svg](dla-evidence.svg)
- [sequence-funnel.svg](sequence-funnel.svg)
- [native-MS-maps.svg](native-MS-maps.svg)

See: https://github.com/pirl-unc/tsarina/issues/208; https://github.com/pirl-unc/oncoref/issues/571; https://github.com/pirl-unc/hitlist/issues/660; https://github.com/pirl-unc/hitlist/issues/661; https://github.com/openvax/mhctools/issues/544; https://github.com/openvax/mhcflurry/issues/490.
