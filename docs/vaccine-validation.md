# Real-data CTA vaccine validation (2026-10-06)

This example uses Hitlist **1.64.7**, OncoRef **1.8.207**, Ensembl r112, MHCflurry and human-only in-vivo Pepsickle. It ranks ten individual proteoform scores using 27 expression cohorts mapped to 22 distinct mortality categories, representing **85.95%** of the reference world mortality shares. Identical full protein sequences consume one slot. See the [method guide](vaccine-design.md) for the reference axis, cohort assumptions and options.

Both examples use the `global54_abc` panel, RNA, HBB/HBB_FI UTRs, 120-nt polyA, limits of 1000 aa / 3500 total nt, padding 0–10 aa in steps of 2, beam width 6 and ten search rounds. Direct and AAY joins were considered; these final designs select direct joins. The [54-allele frequency/provenance audit](vaccine-example/strict/hla_panel.csv) records the population proxy inputs.

**Both designs require junction review.** A lower prediction count is a search result, not evidence that the assembled antigen is free of junctional epitopes. `--require-clean-junctions` makes remaining predictions below 1000 nM fail the command after writing its audit.

## Evidence source and validation scope

The global observations index required a freshness rebuild. Its full-corpus lossless contributor export exceeded the available scratch disk ([Hitlist #643](https://github.com/pirl-unc/hitlist/issues/643)). This example instead uses the supported `VaccineInputs` API with **fresh scoped raw scans**, never the stale observations cache. The released Hitlist scanner/collector performed classification, MS/binding partitioning, cross-source and supplementary deduplication, curated exclusions and human MHC-I filtering. The union snapshot contains **588 observations for 217 peptides** and **1,101 retained contributor links**. All raw source snapshots are recorded by SHA-256; [source-snapshots.json](vaccine-example/source-snapshots.json) includes their logical-row locator conventions.

Final checks independently verified translation, native/assembled ligand coordinates, every final junction-window/allele pair, limits including UTRs/stop/polyA and all 60 output-file hashes. Retained native 8-mers were also checked against every independent non-CTA translated source. Formatting/lint, strict documentation build and the full model-enabled test suite passed (654 tests; two optional Topiary backend skips).

The bundled CSVs below contain the **selected ten candidates**; full command output also contains all unselected candidates, source reference tables, raw queried MS observations, contributor audits and complete manifests. Absolute incidence/mortality counts remain missing wherever the upstream reference lacks sourced counts. Shares are percentages, while p95 prevalence is a fraction.

## Strict CTA design

**1000 aa / 3441 total nt; 8 retained proteoforms in 26 native pieces.** Below-1000 nM window/allele predictions: **1481 → 1088**. GAGE1 and GAGE2A have no qualifying panel MS ligand after 8-mer subtraction. Whole-piece length exclusions and all remaining predictions are retained in the tables.

| rank | name | length_aa | mortality_weighted_score | complete_measurement |
| --- | --- | --- | --- | --- |
| 1 | XAGE1A/XAGE1B | 81 | 0.099315 | True |
| 2 | MAGEA4 | 317 | 0.052540 | True |
| 3 | PRAME | 509 | 0.049262 | True |
| 4 | MAGEA3 | 314 | 0.040509 | True |
| 5 | MAGEA6 | 314 | 0.024097 | True |
| 6 | PAGE2 | 111 | 0.020858 | True |
| 7 | CTAG1A/CTAG1B | 180 | 0.019119 | True |
| 8 | GAGE1 | 117 | 0.017259 | True |
| 9 | GAGE2A | 116 | 0.016378 | True |
| 10 | MAGEA9/MAGEA9B | 315 | 0.013851 | True |

[Selection / gene-protein-transcript IDs](vaccine-example/strict/selection.csv) · [Cancer incidence, mortality, p95 and contributions](vaccine-example/strict/cancer_summary.csv) · [Per-cohort denominators](vaccine-example/strict/cancer_cohorts.csv)

![strict p95 prevalence and world mortality](vaccine-example/strict/cancer-priorities.svg)

| name | raw_aa | specific_aa | ms_supported_aa | padded_aa | assembled_aa | assembled_pieces | retained_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| XAGE1A/XAGE1B | 81 | 73 | 69 | 50 | 32 | 1 | 39.5% |
| MAGEA4 | 317 | 252 | 210 | 176 | 129 | 4 | 40.7% |
| PRAME | 509 | 445 | 445 | 430 | 372 | 7 | 73.1% |
| MAGEA3 | 314 | 233 | 213 | 200 | 154 | 5 | 49.0% |
| MAGEA6 | 314 | 232 | 215 | 201 | 150 | 5 | 47.8% |
| PAGE2 | 111 | 97 | 64 | 54 | 36 | 1 | 32.4% |
| CTAG1A/CTAG1B | 180 | 72 | 64 | 63 | 53 | 1 | 29.4% |
| GAGE1 | 117 | 8 | 0 | 0 | 0 | 0 | 0.0% |
| GAGE2A | 116 | 17 | 0 | 0 | 0 | 0 | 0.0% |
| MAGEA9/MAGEA9B | 315 | 288 | 288 | 166 | 73 | 2 | 23.2% |

[Funnel / dropout reasons](vaccine-example/strict/funnel.csv) · [CTA-specific intervals](vaccine-example/strict/specific_intervals.csv) · [Qualified MS ligands / allele support](vaccine-example/strict/ligands.csv)

![strict sequence retention](vaccine-example/strict/sequence-funnel.svg)

![strict native-to-construct coordinates](vaccine-example/strict/construct-map.svg)

[Construct layers](vaccine-example/strict/layers.csv) · [Remaining junction binders](vaccine-example/strict/junction-binders.csv) · [All final junction predictions (gzipped CSV)](vaccine-example/strict/junctions.csv.gz) · [Boundary cleavage](vaccine-example/strict/cleavage.csv) · [Excluded pieces](vaccine-example/strict/excluded_segments.csv) · [Search history](vaccine-example/strict/search_history.csv)

[Protein FASTA](vaccine-example/strict/protein.fasta) · [CDS/stop FASTA](vaccine-example/strict/cds.fasta) · [Full RNA FASTA](vaccine-example/strict/full.fasta) · [Configuration, versions and hashes](vaccine-example/strict/validation.json)

## Loose CTA design

**1000 aa / 3441 total nt; 8 retained proteoforms in 26 native pieces.** Below-1000 nM window/allele predictions: **1444 → 950**. GAGE1 and GAGE2A have no qualifying panel MS ligand after 8-mer subtraction. Whole-piece length exclusions and all remaining predictions are retained in the tables.

| rank | name | length_aa | mortality_weighted_score | complete_measurement |
| --- | --- | --- | --- | --- |
| 1 | XAGE1A/XAGE1B | 81 | 0.099315 | True |
| 2 | MAGEA4 | 317 | 0.052540 | True |
| 3 | PRAME | 509 | 0.049262 | True |
| 4 | MAGEA3 | 314 | 0.040509 | True |
| 5 | MAGEA6 | 314 | 0.024097 | True |
| 6 | PAGE2 | 111 | 0.020858 | True |
| 7 | CTAG1A/CTAG1B | 180 | 0.019119 | True |
| 8 | PAGE4 | 102 | 0.017912 | True |
| 9 | GAGE1 | 117 | 0.017259 | True |
| 10 | GAGE2A | 116 | 0.016378 | True |

[Selection / gene-protein-transcript IDs](vaccine-example/loose/selection.csv) · [Cancer incidence, mortality, p95 and contributions](vaccine-example/loose/cancer_summary.csv) · [Per-cohort denominators](vaccine-example/loose/cancer_cohorts.csv)

![loose p95 prevalence and world mortality](vaccine-example/loose/cancer-priorities.svg)

| name | raw_aa | specific_aa | ms_supported_aa | padded_aa | assembled_aa | assembled_pieces | retained_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| XAGE1A/XAGE1B | 81 | 73 | 69 | 50 | 38 | 1 | 46.9% |
| MAGEA4 | 317 | 264 | 223 | 189 | 143 | 4 | 45.1% |
| PRAME | 509 | 445 | 445 | 430 | 377 | 7 | 74.1% |
| MAGEA3 | 314 | 241 | 223 | 216 | 175 | 5 | 55.7% |
| MAGEA6 | 314 | 240 | 234 | 220 | 161 | 6 | 51.3% |
| PAGE2 | 111 | 103 | 70 | 54 | 36 | 1 | 32.4% |
| CTAG1A/CTAG1B | 180 | 72 | 64 | 63 | 57 | 1 | 31.7% |
| PAGE4 | 102 | 94 | 89 | 22 | 12 | 1 | 11.8% |
| GAGE1 | 117 | 8 | 0 | 0 | 0 | 0 | 0.0% |
| GAGE2A | 116 | 17 | 0 | 0 | 0 | 0 | 0.0% |

[Funnel / dropout reasons](vaccine-example/loose/funnel.csv) · [CTA-specific intervals](vaccine-example/loose/specific_intervals.csv) · [Qualified MS ligands / allele support](vaccine-example/loose/ligands.csv)

![loose sequence retention](vaccine-example/loose/sequence-funnel.svg)

![loose native-to-construct coordinates](vaccine-example/loose/construct-map.svg)

[Construct layers](vaccine-example/loose/layers.csv) · [Remaining junction binders](vaccine-example/loose/junction-binders.csv) · [All final junction predictions (gzipped CSV)](vaccine-example/loose/junctions.csv.gz) · [Boundary cleavage](vaccine-example/loose/cleavage.csv) · [Excluded pieces](vaccine-example/loose/excluded_segments.csv) · [Search history](vaccine-example/loose/search_history.csv)

[Protein FASTA](vaccine-example/loose/protein.fasta) · [CDS/stop FASTA](vaccine-example/loose/cds.fasta) · [Full RNA FASTA](vaccine-example/loose/full.fasta) · [Configuration, versions and hashes](vaccine-example/loose/validation.json)
