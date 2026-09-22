import pandas as pd
import pytest
from oncoref import cta as oncoref_cta

import tsarina.gene_sets as gene_sets
from tsarina import (
    CTA_by_axes,
    CTA_evidence,
    CTA_excluded_gene_ids,
    CTA_excluded_gene_names,
    CTA_filtered_gene_ids,
    CTA_filtered_gene_names,
    CTA_gene_ids,
    CTA_gene_names,
    CTA_never_expressed_gene_ids,
    CTA_never_expressed_gene_names,
    CTA_placental_restricted_gene_ids,
    CTA_placental_restricted_gene_names,
    CTA_relaxed_reproductive_gene_ids,
    CTA_relaxed_reproductive_gene_names,
    CTA_testis_restricted_gene_ids,
    CTA_testis_restricted_gene_names,
    CTA_unfiltered_gene_ids,
    CTA_unfiltered_gene_names,
    cta_symbol_for_alias,
)
from tsarina.loader import _local_ms_evidence


def test_foundational_cta_helpers_are_direct_oncoref_aliases():
    aliases = {
        "CTA_gene_names": "cta_gene_names",
        "CTA_gene_ids": "cta_gene_ids",
        "CTA_filtered_gene_names": "cta_filtered_gene_names",
        "CTA_filtered_gene_ids": "cta_filtered_gene_ids",
        "CTA_never_expressed_gene_names": "cta_never_expressed_gene_names",
        "CTA_never_expressed_gene_ids": "cta_never_expressed_gene_ids",
        "CTA_unfiltered_gene_names": "cta_unfiltered_gene_names",
        "CTA_unfiltered_gene_ids": "cta_unfiltered_gene_ids",
        "CTA_excluded_gene_names": "cta_excluded_gene_names",
        "CTA_excluded_gene_ids": "cta_excluded_gene_ids",
        "CTA_relaxed_reproductive_gene_names": "cta_relaxed_reproductive_gene_names",
        "CTA_relaxed_reproductive_gene_ids": "cta_relaxed_reproductive_gene_ids",
        "CTA_testis_restricted_gene_names": "cta_testis_restricted_gene_names",
        "CTA_testis_restricted_gene_ids": "cta_testis_restricted_gene_ids",
        "CTA_placental_restricted_gene_names": "cta_placental_restricted_gene_names",
        "CTA_placental_restricted_gene_ids": "cta_placental_restricted_gene_ids",
        "CTA_clinical_target_gene_names": "cta_clinical_target_gene_names",
        "CTA_clinical_target_gene_ids": "cta_clinical_target_gene_ids",
        "cta_symbol_for_alias": "cta_symbol_for_alias",
    }
    for tsarina_name, oncoref_name in aliases.items():
        assert getattr(gene_sets, tsarina_name) is getattr(oncoref_cta, oncoref_name)


def test_ctag2_is_a_clinical_target_excluded_from_the_strict_set():
    """Regression pin: CTAG2/LAGE-1 is oncoref's motivating example for
    exclude_default_keep_clinical -- excluded from the strict CTA set for
    a low-level HPA heart RNA signal, but retained as a known clinical
    target given its NY-ESO-1-family therapeutic history."""
    from tsarina import CTA_clinical_target_gene_names, CTA_gene_names

    assert "CTAG2" in CTA_clinical_target_gene_names()
    assert "CTAG2" not in CTA_gene_names()


def test_cta_set_relationships_are_oncoref_relationships():
    assert CTA_gene_names() == oncoref_cta.cta_gene_names()
    assert CTA_gene_ids() == oncoref_cta.cta_gene_ids()
    assert CTA_filtered_gene_names() == oncoref_cta.cta_filtered_gene_names()
    assert CTA_filtered_gene_ids() == oncoref_cta.cta_filtered_gene_ids()
    assert CTA_never_expressed_gene_names() == CTA_filtered_gene_names() - CTA_gene_names()
    assert CTA_never_expressed_gene_ids() == CTA_filtered_gene_ids() - CTA_gene_ids()
    assert CTA_excluded_gene_names() == CTA_unfiltered_gene_names() - CTA_filtered_gene_names()
    assert CTA_excluded_gene_ids() == CTA_unfiltered_gene_ids() - CTA_filtered_gene_ids()


def test_evidence_row_universe_and_upstream_columns_match_oncoref_exactly():
    actual = CTA_evidence().copy()
    expected = oncoref_cta.cta_evidence().copy()

    assert list(actual["Ensembl_Gene_ID"]) == list(expected["Ensembl_Gene_ID"])
    assert list(actual["Symbol"]) == list(expected["Symbol"])
    pd.testing.assert_frame_equal(
        actual[list(expected.columns)].reset_index(drop=True),
        expected.reset_index(drop=True),
        check_dtype=False,
        check_column_type=False,
    )


def test_tsarina_evidence_adds_only_ms_columns_and_filtered_compatibility_alias():
    actual = CTA_evidence()
    expected_columns = set(oncoref_cta.cta_evidence().columns)
    added = set(actual.columns) - expected_columns
    assert added == {
        "filtered",
        "ms_restriction",
        "ms_healthy_somatic_tissues",
        "ms_pmids",
    }
    assert actual["filtered"].equals(actual["passes_filters"])
    assert not any(column.startswith("ms_") and "count" in column for column in actual)


def test_ms_overlay_is_gene_evidence_not_a_second_cta_definition():
    overlay = _local_ms_evidence()
    assert list(overlay.columns) == [
        "Ensembl_Gene_ID",
        "ms_restriction",
        "ms_healthy_somatic_tissues",
        "ms_pmids",
    ]
    assert overlay["Ensembl_Gene_ID"].is_unique
    assert not {
        "Symbol",
        "passes_filters",
        "never_expressed",
        "specificity_action",
        "specificity_status",
    } & set(overlay.columns)
    assert (
        overlay["ms_restriction"].ne("NO_MS_DATA")
        | overlay["ms_healthy_somatic_tissues"].ne("")
        | overlay["ms_pmids"].ne("")
    ).all()


def test_non_cta_ms_evidence_does_not_expand_oncoref_cta_universe():
    # H1-6 retains generic safety evidence, but oncoref's histone-family rule
    # excludes it from the CTA universe. The left join must not reintroduce it.
    overlay = _local_ms_evidence().set_index("Ensembl_Gene_ID")
    assert overlay.loc["ENSG00000187475", "ms_restriction"] == "RECURRENT_HEALTHY"
    assert "H1-6" not in set(CTA_evidence()["Symbol"])
    assert "H1-6" not in CTA_unfiltered_gene_names()


def test_missing_ms_overlay_rows_receive_explicit_defaults():
    df = CTA_evidence().set_index("Symbol")
    assert df.loc["SUN5", "ms_restriction"] == "NO_MS_DATA"
    assert df.loc["SUN5", "ms_healthy_somatic_tissues"] == ""
    assert df.loc["SUN5", "ms_pmids"] == ""


@pytest.mark.parametrize(
    "symbol, tissues, reliability",
    [
        ("CGB2", "placenta; testis", "Supported"),
        ("CGB3", "placenta; testis", "Enhanced"),
        ("CGB5", "placenta; testis", "Supported"),
        ("CGB7", "placenta; testis", "Supported"),
        ("PSG4", "placenta; testis", "Supported"),
        ("PSG6", "placenta; testis", "Supported"),
        ("PSG7", "placenta; testis", "Supported"),
        ("CT45A5", "testis", "Supported"),
        ("CSH1", "pituitary gland; placenta", "Enhanced"),
        ("GAGE10", "testis", "Supported"),
    ],
)
def test_historical_rna_only_additions_retain_hpa_v23_ihc(symbol, tissues, reliability):
    # #130: RNA-only additions once replaced these real HPA v23 IHC calls
    # with "no data". Verified against normal_tissue, independently of the
    # oncoref evidence frame. This is evidence availability, not CTA admission.
    row = CTA_evidence().set_index("Symbol").loc[symbol]
    assert row["protein_strict_expression"] == tissues
    assert row["protein_reliability"] == reliability
    assert not row["never_expressed"]


@pytest.mark.parametrize(
    "symbol, gene_id, protein_restriction, selected",
    [
        ("SUN5", "ENSG00000167098", "TESTIS", True),
        ("SUN3", "ENSG00000164744", "SOMATIC", False),
    ],
)
def test_sun_candidates_preserve_discordant_rna_and_protein_evidence(
    symbol, gene_id, protein_restriction, selected
):
    row = CTA_evidence().set_index("Symbol").loc[symbol]
    assert row["Ensembl_Gene_ID"] == gene_id
    assert row["rna_restriction"] == "TESTIS"
    assert row["protein_restriction"] == protein_restriction
    assert symbol in CTA_unfiltered_gene_names()
    assert (symbol in CTA_gene_names()) == selected
    assert (symbol in CTA_excluded_gene_names()) != selected


@pytest.mark.parametrize("symbol", ["SUN1", "SUN2", "SPAG4"])
def test_somatic_sun_members_are_not_default_ctas(symbol):
    # A SUN-family name or CT registry designation does not establish
    # reproductive restriction (SPAG4/SUN4 has substantial human pancreas RNA).
    assert symbol not in CTA_gene_names()


def test_csh1_exclusion_is_oncoref_owned():
    row = CTA_evidence().set_index("Symbol").loc["CSH1"]
    upstream = oncoref_cta.cta_evidence().set_index("Symbol").loc["CSH1"]
    assert row["specificity_action"] == upstream["specificity_action"] == "exclude_default"
    assert "CSH1" in CTA_excluded_gene_names()
    assert "CSH1" not in CTA_filtered_gene_names()


def test_hpa_axis_helpers_match_oncoref():
    assert CTA_testis_restricted_gene_names() == oncoref_cta.cta_testis_restricted_gene_names()
    assert CTA_testis_restricted_gene_ids() == oncoref_cta.cta_testis_restricted_gene_ids()
    assert (
        CTA_placental_restricted_gene_names() == oncoref_cta.cta_placental_restricted_gene_names()
    )
    assert CTA_placental_restricted_gene_ids() == oncoref_cta.cta_placental_restricted_gene_ids()
    assert (
        CTA_relaxed_reproductive_gene_names() == oncoref_cta.cta_relaxed_reproductive_gene_names()
    )
    assert CTA_relaxed_reproductive_gene_ids() == oncoref_cta.cta_relaxed_reproductive_gene_ids()


def test_by_axes_matches_oncoref_for_hpa_filters():
    cases = [
        {"restriction": "TESTIS"},
        {"protein_restriction": "TESTIS"},
        {"rna_restriction_level": "STRICT"},
        {"restriction": "TESTIS", "rna_restriction_level": "STRICT"},
        {"restriction_confidence": "HIGH"},
        {"restriction": "TESTIS", "column": "Ensembl_Gene_ID"},
    ]
    for kwargs in cases:
        assert CTA_by_axes(**kwargs) == oncoref_cta.cta_by_axes(**kwargs)


def test_by_axes_adds_explicit_ms_filter_without_changing_membership():
    cancer_only = CTA_by_axes(ms_restriction="CANCER_ONLY")
    recurrent_healthy = CTA_by_axes(ms_restriction="RECURRENT_HEALTHY")
    assert cancer_only
    assert recurrent_healthy
    assert cancer_only <= CTA_filtered_gene_names()
    assert recurrent_healthy <= CTA_filtered_gene_names()


def test_by_axes_unfiltered_still_uses_oncoref_universe():
    all_restrictions = {"TESTIS", "PLACENTAL", "REPRODUCTIVE", "SOMATIC", "NO_DATA"}
    assert (
        CTA_by_axes(restriction=all_restrictions, filtered_only=False)
        == CTA_unfiltered_gene_names()
    )


def test_alias_resolution_is_oncoref_owned():
    assert cta_symbol_for_alias("NY-ESO-1") == "CTAG1B"
    assert cta_symbol_for_alias("CT12.2") == "XAGE2"
    assert cta_symbol_for_alias("not-a-gene") is None
