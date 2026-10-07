"""Exact full-protein MS observations with source context and native coordinates."""

from __future__ import annotations

import pandas as pd

from .vaccine_sequences import peptide_occurrences

TISSUE_COLUMNS = [
    "donor_id",
    "donor_status",
    "source_record_id",
    "sample_alleles",
    "source_dataset",
    "tissue_status",
    "proteoform_key",
    "name",
    "peptide",
    "start",
    "end",
    "source_group",
    "source_tissue",
    "disease",
    "cell_name",
    "src_cell_line",
    "src_ex_vivo",
    "pmid",
    "assay_iri",
    "provenance_id",
    "mhc_restriction",
    "is_monoallelic",
    "restriction_evidence",
    "panel_assignments",
    "specific",
    "retained",
    "cancer_only_in_corpus",
]


def truth(value):
    return str(value).lower() in {"true", "1", "1.0", "yes"}


def sample_group(row):
    """Prefer explicit disease/source curation; unknown is never healthy."""
    if truth(row.get("src_cancer")):
        return "cancer"
    if truth(row.get("src_adjacent_to_tumor")):
        return "tumor_adjacent"
    if truth(row.get("src_healthy_reproductive")):
        return "healthy_reproductive"
    if truth(row.get("src_healthy_tissue")):
        return "healthy_nonreproductive"
    return "unknown_or_other"


def verified_normal_ms(atlas_dir, peptides, min_donors=1):
    """Use primary donor-resolved Atlas HLA-I evidence, not IEDB tissue flags.

    Atlas tissues are nonmalignant primary autopsy tissue, not necessarily
    disease-free. Scope is verified heart/brain/lung. No reproductive tissue
    exceptions apply to these organs under either CTA definition. The explicit
    donor threshold is recorded; one donor is conservative, not replication.
    """
    try:
        from hitlist.tissue_blacklist import build_tissue_blacklist, load_atlas_tissue_evidence
    except ImportError as error:
        raise ImportError("Verified normal-MS evidence requires Hitlist >=1.65.1") from error

    observations, provenance = load_atlas_tissue_evidence(atlas_dir)
    query_kmers = {p[i : i + 8] for p in peptides for i in range(len(p) - 7)}
    # A longer normal ligand need not occur intact in a CTA to veto its shared
    # 8-mer. Search every source ligand length before applying the donor gate.
    matched = {
        p
        for p in observations.loc[observations.mhc_class.eq("I"), "peptide"].unique()
        if any(p[i : i + 8] in query_kmers for i in range(len(p) - 7))
    }
    observations = observations[
        observations.mhc_class.eq("I") & observations.peptide.isin(matched)
    ].copy()
    summary, audit = build_tissue_blacklist(observations, min_donors=min_donors)
    qualified = audit[audit.qualified].copy()
    qualified["src_healthy_tissue"] = True
    qualified["src_cell_line"] = False
    qualified["src_cancer"] = False
    qualified["src_adjacent_to_tumor"] = False
    qualified["src_ex_vivo"] = True
    qualified["species"] = "Homo sapiens"
    qualified["assay_method"] = "cellular MHC/mass spectrometry"
    qualified["source"] = "hla_ligand_atlas"
    qualified["provenance_id"] = qualified.source_record_id
    qualified["is_monoallelic"] = False
    qualified["restriction_evidence"] = "donor_genotype_not_measured_restriction"
    qualified["cell_name"] = "Primary nonmalignant autopsy tissue"
    excluded = qualified[qualified.peptide.isin(summary.loc[summary.blacklisted, "peptide"])].copy()
    provenance.update(
        {
            "min_donors": min_donors,
            "class": "I",
            "n_query_peptides": len(peptides),
            "query_match": "any shared 8-mer, including longer source ligands",
            "n_qualifying_observations": len(qualified),
            "n_blacklisted_peptides": int(summary.blacklisted.sum()),
        }
    )
    return excluded, qualified, summary, audit, provenance


def tissue_map(proteins, observations, intervals, layers, assignments=None):
    """Map each observation to all exact native occurrences, not a gene claim.

    Membership in a source protein does not prove the peptide originated from
    that gene. Shared sequences remain visible. No normal-tissue exclusion is
    performed here; overlaps are an audit result, not a claim of target safety.
    """
    assignments = assignments if assignments is not None else pd.DataFrame()
    assignment_keys = {}
    if not assignments.empty:
        for row in assignments.itertuples(index=False):
            key = (row.peptide, str(row.provenance_id))
            assignment_keys.setdefault(key, set()).add(f"{row.allele} ({row.evidence_tier})")
    rows = []
    for protein in proteins.itertuples(index=False):
        specific = intervals[intervals.proteoform_key.eq(protein.proteoform_key)]
        retained = [
            layer for layer in layers if layer.get("proteoform_key") == protein.proteoform_key
        ]
        for observation in observations.to_dict("records"):
            peptide = str(observation["peptide"])
            for start, end in peptide_occurrences(protein.sequence, peptide):
                row = {key: observation.get(key) for key in TISSUE_COLUMNS}
                row.update(
                    {
                        "proteoform_key": protein.proteoform_key,
                        "name": protein.name,
                        "peptide": peptide,
                        "start": start,
                        "end": end,
                        "source_group": sample_group(observation),
                        "specific": any(
                            r.start <= start and end <= r.end for r in specific.itertuples()
                        ),
                        "retained": any(
                            r["native_start"] <= start and end <= r["native_end"] for r in retained
                        ),
                        "panel_assignments": "; ".join(
                            sorted(
                                assignment_keys.get(
                                    (peptide, str(observation.get("provenance_id"))), set()
                                )
                            )
                        ),
                    }
                )
                rows.append(row)
    frame = pd.DataFrame(rows, columns=TISSUE_COLUMNS)
    if not frame.empty:
        groups = frame.groupby("peptide").source_group.agg(set)
        # Unknown/other or normal observations preclude a cancer-only label.
        frame["cancer_only_in_corpus"] = frame.peptide.map(groups.map(lambda s: s == {"cancer"}))
    return frame
