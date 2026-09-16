# Sid Sijbrandij personalized neoantigen vaccine overlap table

`vaccine_overlap_summary.tsv` is copied without modification from the
public spreadsheet linked at https://osteosarc.com/vaccines (the
"vaccine overlap summary" Google Sheet), itself derived from public,
openly (CC0) shared personal cancer-genomics research data -- see
https://osteosarc.com and `vaxrank/tests/data/osteosarc/README.md` for the
wider provenance of this dataset and its licensing.

SHA-256 of the committed file: `dbd89d6f3a34a127dc533163cea801e4afec8cb8b773c1f9bd1473eb580cb965`

## What it is

One row per somatic mutation nominated across four independent
personalized-vaccine efforts for one patient: an mRNA vaccine, three
sequential JLF peptide vaccine designs (V1/V2/V3), and a CeGaT vaccine.
Columns record, per construct, whether a target was included and its
peptide sequence (a short "minimal epitope" and/or a longer flanking
"full peptide"), plus ELISPOT immunogenicity testing results where run.

`TECPR1*` is a second row for `TECPR1` because that gene had distinct
long peptides for MHC class I and class II in the JLF vaccines (see the
footnote in the original sheet). `KRT18`'s CeGaT peptide is `?` in the
source and is preserved as-is rather than guessed. `NME1`'s peptide was
"transferred from data shared with PFO" per the original sheet's footnote.

## Why it's here

Real, patient-specific neoantigen sequences with known immunogenicity
outcomes are useful test data for anything that reasons about
mutation-derived (as opposed to shared/CTA) tumor antigens and their
overlap with public mass-spec evidence -- see
`tsarina/neoantigen_evidence_plots.py`. This is deliberately kept as test
fixture data, not packaged production data: unlike `tsarina/data/*.csv`,
nothing in the shipped library depends on this file existing.

This is personal data about one named, real individual. It is used here
only because it was openly and explicitly shared for cancer research
reuse; treat it accordingly, and do not extend this pattern to data that
was not offered on the same terms.
