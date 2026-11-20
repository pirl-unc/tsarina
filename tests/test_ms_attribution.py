"""Tests for reusable restriction-aware MS attribution across target families (#216)."""

from __future__ import annotations

import pandas as pd
import pytest

from tsarina import (
    attribute_ms_evidence,
    personalized_targets,
    summarize_ms_attribution,
    target_peptides,
)
from tsarina.ms_attribution import (
    ANY_MHC_PEPTIDE,
    CLASS_OR_SEROTYPE_ONLY,
    DONOR_SET_AMBIGUOUS,
    EVIDENCE_UNAVAILABLE,
    EXACT_REQUESTED_MHC,
    NO_MS_MATCH,
    PREDICTED_RESTRICTION,
)
from tsarina.mutations import mutant_iedb_overlap
from tsarina.viral import viral_iedb_overlap


def _mixed_family_candidates() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "peptide": "SLLQHLIGL",
                "category": "cta",
                "source": "PRAME",
                "source_detail": "ENSG00000185686",
            },
            {
                "peptide": "YMLDLQPET",
                "category": "viral",
                "source": "hpv16",
                "source_detail": "E7",
            },
            {
                "peptide": "VVVGADGVGK",
                "category": "mutant",
                "source": "KRAS G12D",
                "source_detail": "G12D",
            },
            {
                "peptide": "KLQCVDLHV",
                "category": "cta",
                "source": "MAGEA4",
                "source_detail": "ENSG00000147381",
            },
            {
                "peptide": "ACDEFGHIK",
                "category": "cta",
                "source": "CTAG1B",
                "source_detail": "ENSG00000184033",
            },
            {
                "peptide": "MNPQRSTVW",
                "category": "cta",
                "source": "NOMATCH",
                "source_detail": "ENSG00000000001",
            },
        ]
    )


def _sample_observations() -> pd.DataFrame:
    df = pd.DataFrame(
        [
            # Exact monoallelic HLA-A*02:01 on CTA peptide
            {
                "peptide": "SLLQHLIGL",
                "mhc_restriction": "HLA-A*02:01",
                "is_monoallelic": True,
                "mhc_allele_provenance": "exact",
                "assay_method": "cellular MHC/mass spectrometry",
                "qualitative_measurement": "Positive",
                "source_tissue": "Melanoma",
                "disease": "melanoma",
                "src_cancer": True,
                "src_healthy_tissue": False,
                "src_cell_line": False,
                "src_ex_vivo": True,
                "pmid": "11111111",
                "assay_iri": "http://iedb.org/assay/1",
                "provenance_id": "iedb:1",
                "source": "iedb",
            },
            # Ambiguous multi-allelic donor set on viral peptide (never exact!)
            {
                "peptide": "YMLDLQPET",
                "mhc_restriction": "HLA-A*02:01;HLA-B*07:02",
                "mhc_allele_set": "HLA-A*02:01;HLA-B*07:02",
                "is_monoallelic": False,
                "mhc_allele_provenance": "sample_allele_match",
                "assay_method": "cellular MHC/mass spectrometry",
                "qualitative_measurement": "Positive",
                "source_tissue": "Cervix",
                "disease": "cervical carcinoma",
                "src_cancer": True,
                "src_healthy_tissue": False,
                "src_cell_line": True,
                "src_ex_vivo": False,
                "pmid": "22222222",
                "assay_iri": "http://iedb.org/assay/2",
                "provenance_id": "cedar:2",
                "source": "cedar",
            },
            # Class-only / serotype-only restriction on mutant peptide
            {
                "peptide": "VVVGADGVGK",
                "mhc_restriction": "HLA class I",
                "is_monoallelic": False,
                "mhc_allele_provenance": "class_only",
                "assay_method": "cellular MHC/mass spectrometry",
                "qualitative_measurement": "Positive",
                "source_tissue": "Pancreas",
                "disease": "PDAC",
                "src_cancer": True,
                "src_healthy_tissue": False,
                "pmid": "33333333",
                "assay_iri": "http://iedb.org/assay/3",
                "provenance_id": "cedar:3",
                "source": "cedar",
            },
            # Unrelated allele (HLA-B*44:02) on MAGEA4 -> any_mhc_peptide unless predicted
            {
                "peptide": "KLQCVDLHV",
                "mhc_restriction": "HLA-B*44:02",
                "is_monoallelic": True,
                "mhc_allele_provenance": "exact",
                "assay_method": "cellular MHC/mass spectrometry",
                "qualitative_measurement": "Positive",
                "source_tissue": "Sarcoma",
                "disease": "synovial sarcoma",
                "src_cancer": True,
                "src_healthy_tissue": False,
                "pmid": "44444444",
                "source": "iedb",
            },
            # Longer nested ligand containing ACDEFGHIK: MUST NOT transfer!
            {
                "peptide": "KKACDEFGHIKRR",
                "mhc_restriction": "HLA-A*02:01",
                "is_monoallelic": True,
                "assay_method": "cellular MHC/mass spectrometry",
                "qualitative_measurement": "Positive",
                "pmid": "55555555",
                "source": "iedb",
            },
            # Non-MS / unknown rows on MNPQRSTVW: MUST be rejected by positive-MS gate!
            {
                "peptide": "MNPQRSTVW",
                "mhc_restriction": "HLA-A*02:01",
                "assay_method": "x-ray crystallography",
                "qualitative_measurement": "Positive",
                "source": "iedb",
            },
            {
                "peptide": "MNPQRSTVW",
                "mhc_restriction": "HLA-A*02:01",
                "assay_method": "mass spectrometry",
                "qualitative_measurement": "unknown",
                "source": "iedb",
            },
        ]
    )
    df.attrs["dataset_version"] = "hitlist-test-2026-11"
    return df


def test_attribute_ms_evidence_separates_all_tiers_and_preserves_provenance():
    cands = _mixed_family_candidates()
    obs = _sample_observations()
    predicted = pd.DataFrame([{"peptide": "KLQCVDLHV", "allele": "HLA-A*02:01"}])

    attributed = attribute_ms_evidence(
        cands,
        obs,
        requested_alleles=["HLA-A*02:01"],
        predicted_assignments=predicted,
        evidence_scope="all",
    )

    by_pep = attributed.set_index("peptide")
    assert by_pep.loc["SLLQHLIGL", "ms_attribution_tier"] == EXACT_REQUESTED_MHC
    assert by_pep.loc["SLLQHLIGL", "mhc_resolution"] == "exact_allele"
    assert by_pep.loc["SLLQHLIGL", "matched_requested_alleles"] == "HLA-A*02:01"
    assert by_pep.loc["SLLQHLIGL", "dataset_version"] == "hitlist-test-2026-11"
    assert by_pep.loc["SLLQHLIGL", "category"] == "cta"

    # Multi-allelic donor set must remain donor_set_ambiguous, never exact_requested_mhc.
    assert by_pep.loc["YMLDLQPET", "ms_attribution_tier"] == DONOR_SET_AMBIGUOUS
    assert by_pep.loc["YMLDLQPET", "mhc_resolution"] == "donor_set"
    assert by_pep.loc["YMLDLQPET", "matched_requested_alleles"] == "HLA-A*02:01"
    assert "HLA-B*07:02" in by_pep.loc["YMLDLQPET", "donor_set_alleles"]

    # Class-only restriction stays class_or_serotype_only (no fabricated 4-digit allele).
    assert by_pep.loc["VVVGADGVGK", "ms_attribution_tier"] == CLASS_OR_SEROTYPE_ONLY
    assert by_pep.loc["VVVGADGVGK", "mhc_resolution"] == "class_only"
    assert by_pep.loc["VVVGADGVGK", "matched_requested_alleles"] == ""

    # Explicit computational prediction is strictly separated as predicted_restriction.
    assert by_pep.loc["KLQCVDLHV", "ms_attribution_tier"] == PREDICTED_RESTRICTION
    assert by_pep.loc["KLQCVDLHV", "mhc_restriction_raw"] == "HLA-B*44:02"

    # Nested peptide from longer ligand KKACDEFGHIKRR does NOT inherit MS support.
    assert by_pep.loc["ACDEFGHIK", "ms_evidence_status"] == "no_match"
    assert by_pep.loc["ACDEFGHIK", "ms_attribution_tier"] == NO_MS_MATCH

    # Structural and unknown-outcome rows are excluded -> no_ms_match.
    assert by_pep.loc["MNPQRSTVW", "ms_evidence_status"] == "no_match"
    assert by_pep.loc["MNPQRSTVW", "ms_attribution_tier"] == NO_MS_MATCH

    # Without predicted_assignments, KLQCVDLHV falls cleanly into any_mhc_peptide.
    no_pred = attribute_ms_evidence(
        cands, obs, requested_alleles=["HLA-A*02:01"], evidence_scope="all"
    ).set_index("peptide")
    assert no_pred.loc["KLQCVDLHV", "ms_attribution_tier"] == ANY_MHC_PEPTIDE


def test_evidence_unavailable_stays_distinct_from_no_match():
    cands = _mixed_family_candidates().head(2)
    unloaded = attribute_ms_evidence(cands, None, evidence_loaded=False)
    assert set(unloaded["ms_evidence_status"]) == {"not_loaded"}
    assert set(unloaded["ms_attribution_tier"]) == {EVIDENCE_UNAVAILABLE}

    empty_loaded = attribute_ms_evidence(
        cands, pd.DataFrame(columns=["peptide"]), evidence_loaded=True
    )
    assert set(empty_loaded["ms_evidence_status"]) == {"no_match"}
    assert set(empty_loaded["ms_attribution_tier"]) == {NO_MS_MATCH}


def test_evidence_scope_filters_requested_mhc_vs_any_mhc():
    cands = _mixed_family_candidates()
    obs = _sample_observations()

    req_only = attribute_ms_evidence(
        cands,
        obs,
        requested_alleles=["HLA-A*02:01"],
        evidence_scope="requested_mhc",
    )
    assert set(req_only["peptide"]) == {"SLLQHLIGL", "YMLDLQPET"}

    any_mhc = attribute_ms_evidence(
        cands,
        obs,
        requested_alleles=["HLA-A*02:01"],
        evidence_scope="any_mhc",
    )
    assert set(any_mhc["peptide"]) == {"SLLQHLIGL", "YMLDLQPET", "VVVGADGVGK", "KLQCVDLHV"}

    with pytest.raises(ValueError, match="requested_alleles"):
        attribute_ms_evidence(cands, obs, evidence_scope="requested_mhc")

    summary = summarize_ms_attribution(
        attribute_ms_evidence(cands, obs, requested_alleles=["HLA-A*02:01"])
    ).set_index("peptide")
    assert summary.loc["SLLQHLIGL", "ms_exact_mhc_hits"] == 1
    assert summary.loc["SLLQHLIGL", "ms_donor_set_hits"] == 0
    assert summary.loc["YMLDLQPET", "ms_exact_mhc_hits"] == 0
    assert summary.loc["YMLDLQPET", "ms_donor_set_hits"] == 1
    assert summary.loc["VVVGADGVGK", "ms_class_or_serotype_labels"] != ""


def test_personalize_viral_and_mutant_share_attribution_contract(monkeypatch):
    obs = _sample_observations()
    custom = _mixed_family_candidates()[["peptide", "category", "source", "source_detail"]].copy()

    def _fake_score(peptides, alleles, predictor="mhcflurry"):
        return pd.DataFrame(
            [
                {
                    "peptide": p,
                    "allele": a,
                    "presentation_percentile": 0.3,
                    "presentation_score": 0.9,
                    "affinity_nm": 50.0,
                }
                for p in peptides
                for a in alleles
            ]
        )

    monkeypatch.setattr("tsarina.scoring.score_presentation", _fake_score)

    df = personalized_targets(
        hla_alleles=["HLA-A*02:01"],
        custom_targets=custom,
        ms_observations=obs,
        show_progress=False,
    ).set_index("peptide")
    assert df.loc["SLLQHLIGL", "ms_attribution_tier"] == EXACT_REQUESTED_MHC
    assert df.loc["YMLDLQPET", "ms_attribution_tier"] == DONOR_SET_AMBIGUOUS
    # Scored best_allele upgrades KLQCVDLHV to predicted_restriction, not exact_requested_mhc!
    assert df.loc["KLQCVDLHV", "ms_attribution_tier"] == PREDICTED_RESTRICTION
    assert df.loc["KLQCVDLHV", "ms_exact_mhc_hits"] == 0
    assert df.loc["ACDEFGHIK", "ms_evidence_status"] == "no_match"

    # Viral overlap with requested_alleles & offline observations
    viral_df = viral_iedb_overlap(
        proteins={"E7": "AAAYMLDLQPETAAA"},
        lengths=(9,),
        requested_alleles=["HLA-A*02:01"],
        observations=obs,
    )
    ym = viral_df.loc[viral_df["peptide"].eq("YMLDLQPET")].iloc[0]
    assert ym["ms_attribution_tier"] == DONOR_SET_AMBIGUOUS
    assert ym["ms_donor_set_hits"] == 1
    assert ym["ms_exact_mhc_hits"] == 0

    # Mutant overlap with offline wt_sequence + observations
    mut_df = mutant_iedb_overlap(
        mutations=[
            {
                "label": "KRAS G12D",
                "gene": "KRAS",
                "gene_id": "ENSG00000133703",
                "protein_position": 6,
                "ref_aa": "G",
                "alt_aa": "D",
                "wt_sequence": "VVVGAGGVGK",
            }
        ],
        lengths=(10,),
        requested_alleles=["HLA-A*02:01"],
        observations=obs,
        host_proteins={"KRAS_WT": "VVVGAGGVGK"},
    )
    row = mut_df.loc[mut_df["peptide"].eq("VVVGADGVGK")].iloc[0]
    assert row["ms_attribution_tier"] == CLASS_OR_SEROTYPE_ONLY
    assert row["has_iedb_hit"]

    # Unified target_peptides with custom_targets + requested_alleles
    tp = target_peptides(
        cta=False,
        mutations=False,
        custom_targets=custom,
        requested_alleles=["HLA-A*02:01"],
        evidence_scope="requested_mhc",
        observations=obs,
    )
    assert set(tp["peptide"]) == {"SLLQHLIGL", "YMLDLQPET"}
