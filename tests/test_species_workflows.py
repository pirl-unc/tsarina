"""Tests for species-scoped target, self-exclusion, and MS-evidence workflows (#217)."""

from __future__ import annotations

import pandas as pd
import pytest

from tsarina import (
    SpeciesContext,
    allele_prediction_capability,
    attribute_ms_evidence,
    personalized_targets,
    resolve_species_context,
    target_peptides,
)
from tsarina.ms_attribution import DONOR_SET_AMBIGUOUS, EXACT_REQUESTED_MHC
from tsarina.mutations import mutant_peptides
from tsarina.species import validate_alleles_for_species
from tsarina.viral import host_exclusive_viral_peptides


def test_resolve_species_context_human_mouse_dog_and_transgenic():
    human = resolve_species_context("human", ensembl_release=112)
    assert isinstance(human, SpeciesContext)
    assert human.host_species == "Homo sapiens"
    assert human.host_taxon == 9606
    assert human.assembly == "GRCh38"
    assert not human.is_cross_species_mhc

    mouse = resolve_species_context("Mus musculus", ensembl_release=112)
    assert mouse.host_species == "Mus musculus"
    assert mouse.host_taxon == 10090
    assert mouse.pyensembl_species == "mus_musculus"
    assert mouse.assembly == "GRCm39"

    dog = resolve_species_context("Canis lupus familiaris", ensembl_release=112)
    assert dog.host_species == "Canis lupus familiaris"
    assert dog.host_taxon == 9615
    assert dog.pyensembl_species == "canis_familiaris"
    assert dog.assembly == "ROS_Cfam_1.0"

    # Transgenic mouse host with human MHC and viral antigen source stays distinct.
    transgenic = resolve_species_context(
        "mouse",
        antigen_source_species="Human papillomavirus type 16",
        antigen_source_taxon=333760,
        mhc_species="Homo sapiens",
        ensembl_release=112,
        reference_source="orfdata",
        reference_key="orfdata.mm39.test",
    )
    assert transgenic.host_species == "Mus musculus"
    assert transgenic.host_taxon == 10090
    assert transgenic.mhc_species == "Homo sapiens"
    assert transgenic.mhc_taxon == 9606
    assert transgenic.antigen_source_species == "Human papillomavirus type 16"
    assert transgenic.is_cross_species_mhc is True
    assert transgenic.reference_key == "orfdata.mm39.test"


def test_species_mismatch_and_human_catalog_transfer_rejected():
    with pytest.raises(ValueError, match="Conflicting host_taxon"):
        resolve_species_context("mouse", host_taxon=9606)

    # Human HLA allele cannot silently run against mouse MHC species.
    with pytest.raises(ValueError, match="disagrees with expected mhc_species"):
        validate_alleles_for_species(["HLA-A*02:01"], "Mus musculus")

    # Mouse H-2 allele cannot silently run against human MHC species.
    with pytest.raises(ValueError, match="disagrees with expected mhc_species"):
        validate_alleles_for_species(["H-2-Kb"], "Homo sapiens")

    # Dog DLA allele matches domestic dog cleanly.
    validated_dog = validate_alleles_for_species(["DLA-88*508:01"], "Canis lupus familiaris")
    assert validated_dog[0]["normalized_allele"] == "DLA-88*508:01"

    # Human CTA and HOTSPOT_MUTATIONS catalogs must never transfer silently to mouse or dog.
    with pytest.raises(ValueError, match="Human HOTSPOT_MUTATIONS"):
        mutant_peptides(host_species="mouse")

    with pytest.raises(ValueError, match="Human CTA catalog"):
        target_peptides(cta=True, mutations=False, host_species="dog")

    with pytest.raises(ValueError, match="Human HOTSPOT_MUTATIONS"):
        target_peptides(cta=False, mutations=True, host_species="mouse")

    with pytest.raises(ValueError, match="Human OncoRef CTA"):
        personalized_targets(
            hla_alleles=["H-2-Kb"],
            cta_expression={"PRAME": 50.0},
            host_species="mouse",
            score_presentation=False,
            skip_ms_evidence=True,
        )


def test_mouse_h2_offline_workflow_across_cta_viral_and_mutation():
    """Offline mouse (Mus musculus, H-2-K*b / H-2-D*b) workflow covering custom
    CTA, viral self-exclusion against mouse proteins, mouse somatic mutation,
    exact vs donor-set MS attribution, and explicit non-human capability audit."""
    mouse_host_proteins = {
        "ENSMUSP_SELF": "MMMMMSIINFEKLSELFSELF",  # self protein contains SELFSELF
        "ENSMUSP_KRAS_WT": "MTEYKLVVVGAGGVGKSALTIQ",
    }

    # 1. Viral self-exclusion against mouse host proteome (not human!)
    viral_df = host_exclusive_viral_peptides(
        proteins={"mCMV_IE1": "AAAYPHFMPTNLAAASELFSELF"},
        lengths=(9,),
        host_species="mouse",
        host_proteins=mouse_host_proteins,
    )
    assert "YPHFMPTNL" in set(viral_df["peptide"])
    assert "SELFSELF" not in set(viral_df["peptide"])
    assert viral_df.attrs["species_context"]["host_taxon"] == 10090

    # 2. Mouse somatic mutation screened against mouse host proteome
    mut_df = mutant_peptides(
        mutations=[
            {
                "label": "mKras G12D",
                "gene": "Kras",
                "gene_id": "ENSMUSG00000030265",
                "protein_position": 12,
                "ref_aa": "G",
                "alt_aa": "D",
                "wt_sequence": mouse_host_proteins["ENSMUSP_KRAS_WT"],
            }
        ],
        lengths=(9,),
        host_species="mouse",
        host_proteins=mouse_host_proteins,
    )
    assert "LVVVGADGV" in set(mut_df["peptide"])
    assert mut_df.attrs["species_context"]["host_species"] == "Mus musculus"

    # 3. Combined mouse targets + H-2 MS observations
    mouse_targets = pd.DataFrame(
        [
            {
                "peptide": "SIINFEKL",
                "category": "cta",
                "source": "mCta1",
                "source_detail": "ENSMUSG1",
            },
            {"peptide": "YPHFMPTNL", "category": "viral", "source": "mCMV", "source_detail": "IE1"},
            {
                "peptide": "LVVVGADGV",
                "category": "mutant",
                "source": "mKras G12D",
                "source_detail": "G12D",
            },
        ]
    )
    mouse_obs = pd.DataFrame(
        [
            {
                "peptide": "SIINFEKL",
                "mhc_restriction": "H-2-Kb",
                "is_monoallelic": True,
                "assay_method": "cellular MHC/mass spectrometry",
                "qualitative_measurement": "Positive",
                "src_cancer": True,
                "src_healthy_tissue": False,
                "pmid": "90000001",
                "source": "iedb",
            },
            {
                "peptide": "YPHFMPTNL",
                "mhc_restriction": "H-2-Kb;H-2-Db",
                "mhc_allele_set": "H-2-Kb;H-2-Db",
                "is_monoallelic": False,
                "mhc_allele_provenance": "sample_allele_match",
                "assay_method": "cellular MHC/mass spectrometry",
                "qualitative_measurement": "Positive",
                "src_cancer": True,
                "src_healthy_tissue": False,
                "pmid": "90000002",
                "source": "iedb",
            },
        ]
    )

    result = personalized_targets(
        hla_alleles=["H-2-Kb", "H-2-Db"],
        custom_targets=mouse_targets,
        ms_observations=mouse_obs,
        host_species="mouse",
        score_presentation=False,
        drop_weak_tier=False,
        show_progress=False,
    )
    by_pep = result.set_index("peptide")
    assert by_pep.loc["SIINFEKL", "ms_attribution_tier"] == EXACT_REQUESTED_MHC
    assert by_pep.loc["YPHFMPTNL", "ms_attribution_tier"] == DONOR_SET_AMBIGUOUS
    assert by_pep.loc["LVVVGADGV", "ms_evidence_status"] == "no_match"
    assert result.attrs["species_context"]["host_taxon"] == 10090
    assert result.attrs["species_context"]["mhc_species"] == "Mus musculus"

    caps = result.attrs["allele_prediction_capability"]
    assert not caps["calibrated"].any()
    assert set(caps["normalized_allele"]) == {"H2-K*b", "H2-D*b"}


def test_dog_dla_capability_and_attribution_preserve_reviewed_tiers():
    reviewed_caps = pd.DataFrame(
        [
            {"allele": "DLA-88*508:01", "tier": "empirically_supported", "reason": "frozen_bundle"},
            {"allele": "DLA-88*001:01", "tier": "sequence_extrapolated", "reason": "frozen_bundle"},
        ]
    )
    caps = allele_prediction_capability(
        ["DLA-88*508:01", "DLA-88*001:01", "DLA-88*999:01"],
        predictor="frozen",
        mhc_species="Canis lupus familiaris",
        capability_table=reviewed_caps,
    ).set_index("allele")
    assert caps.loc["DLA-88*508:01", "capability_tier"] == "empirically_supported"
    assert bool(caps.loc["DLA-88*508:01", "calibrated"]) is True
    assert caps.loc["DLA-88*001:01", "capability_tier"] == "sequence_extrapolated"
    assert bool(caps.loc["DLA-88*001:01", "calibrated"]) is False
    assert caps.loc["DLA-88*999:01", "capability_tier"] == "unverified"

    dog_obs = pd.DataFrame(
        [
            {
                "peptide": "KLMNPQRST",
                "mhc_restriction": "DLA-88*508:01",
                "is_monoallelic": True,
                "assay_method": "cellular MHC/mass spectrometry",
                "qualitative_measurement": "Positive",
                "src_cancer": True,
                "src_healthy_tissue": False,
                "pmid": "88888888",
                "source": "supplement",
            }
        ]
    )
    attributed = attribute_ms_evidence(
        ["KLMNPQRST"],
        dog_obs,
        requested_alleles=["DLA-88*508:01"],
        host_species="dog",
    )
    assert attributed.iloc[0]["ms_attribution_tier"] == EXACT_REQUESTED_MHC
    assert attributed.iloc[0]["host_species"] == "Canis lupus familiaris"
    assert attributed.attrs["species_context"]["host_taxon"] == 9615
