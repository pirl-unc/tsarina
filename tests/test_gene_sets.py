from tsarina import (
    CTA_by_axes,
    CTA_evidence,
    CTA_excluded_gene_ids,
    CTA_excluded_gene_names,
    CTA_filtered_gene_ids,
    CTA_filtered_gene_names,
    CTA_gene_ids,
    CTA_gene_names,
    CTA_never_expressed_gene_names,
    CTA_placental_restricted_gene_names,
    CTA_testis_restricted_gene_ids,
    CTA_testis_restricted_gene_names,
    CTA_unfiltered_gene_ids,
    CTA_unfiltered_gene_names,
)


def test_gene_names_nonempty():
    assert len(CTA_gene_names()) == 257


def test_gene_ids_nonempty():
    assert len(CTA_gene_ids()) == 257


def test_expressed_is_strict_subset_of_filtered():
    assert CTA_gene_names() < CTA_filtered_gene_names()
    assert CTA_gene_ids() < CTA_filtered_gene_ids()


def test_filtered_is_strict_subset_of_unfiltered():
    assert CTA_filtered_gene_names() < CTA_unfiltered_gene_names()
    assert CTA_filtered_gene_ids() < CTA_unfiltered_gene_ids()


def test_expressed_and_never_expressed_partition_filtered():
    expressed = CTA_gene_names()
    never_expr = CTA_never_expressed_gene_names()
    filtered = CTA_filtered_gene_names()
    assert expressed & never_expr == set()
    assert expressed | never_expr == filtered


def test_filtered_and_excluded_partition_unfiltered():
    filtered = CTA_filtered_gene_names()
    excluded = CTA_excluded_gene_names()
    unfiltered = CTA_unfiltered_gene_names()
    assert filtered & excluded == set()
    assert filtered | excluded == unfiltered


def test_excluded_ids_consistent():
    assert CTA_excluded_gene_ids() == CTA_unfiltered_gene_ids() - CTA_filtered_gene_ids()


def test_evidence_row_count_matches_unfiltered():
    df = CTA_evidence()
    assert len(df) == len(CTA_unfiltered_gene_names())


def test_evidence_has_expected_columns():
    df = CTA_evidence()
    expected = [
        "Symbol",
        "Ensembl_Gene_ID",
        "source_databases",
        "protein_reproductive",
        "protein_reliability",
        "rna_reproductive",
        "rna_deflated_reproductive_frac",
        "filtered",
        "never_expressed",
        "rna_max_ntpm",
        "protein_restriction",
        "rna_restriction",
        "rna_restriction_level",
        "rna_testis_ntpm",
        "restriction",
        "restriction_confidence",
        "ms_restriction",
    ]
    for col in expected:
        assert col in df.columns, f"Missing column: {col}"


def test_magea4_is_expressed_cta():
    assert "MAGEA4" in CTA_gene_names()


def test_magea4_id_is_expressed_cta():
    assert "ENSG00000147381" in CTA_gene_ids()


# ── Restriction accessor tests ────────────────────────────────────────────


def test_testis_restricted_nonempty():
    assert len(CTA_testis_restricted_gene_names()) >= 200


def test_testis_restricted_is_subset_of_filtered():
    assert CTA_testis_restricted_gene_names() <= CTA_filtered_gene_names()


def test_testis_restricted_ids_match_names():
    assert len(CTA_testis_restricted_gene_ids()) == len(CTA_testis_restricted_gene_names())


def test_placental_restricted_nonempty():
    assert len(CTA_placental_restricted_gene_names()) >= 3


# ── CTA_by_axes tests ────────────────────────────────────────────────────


def test_by_axes_restriction_testis():
    assert CTA_by_axes(restriction="TESTIS") == CTA_testis_restricted_gene_names()


def test_by_axes_all_restrictions_covers_filtered():
    all_r = CTA_by_axes(restriction={"TESTIS", "PLACENTAL", "REPRODUCTIVE", "SOMATIC", "NO_DATA"})
    assert all_r == CTA_filtered_gene_names()


def test_by_axes_unfiltered_covers_all():
    all_r = CTA_by_axes(
        restriction={"TESTIS", "PLACENTAL", "REPRODUCTIVE", "SOMATIC", "NO_DATA"},
        filtered_only=False,
    )
    assert all_r == CTA_unfiltered_gene_names()


def test_by_axes_protein_restriction():
    pr = CTA_by_axes(protein_restriction="TESTIS")
    assert len(pr) >= 150


def test_by_axes_rna_restriction_level():
    strict = CTA_by_axes(rna_restriction_level="STRICT")
    assert len(strict) >= 150


def test_by_axes_combined():
    testis_strict = CTA_by_axes(restriction="TESTIS", rna_restriction_level="STRICT")
    assert len(testis_strict) > 0
    assert testis_strict <= CTA_testis_restricted_gene_names()


def test_by_axes_confidence():
    high = CTA_by_axes(restriction_confidence="HIGH")
    assert len(high) >= 150


def test_by_axes_gene_id_column():
    ids = CTA_by_axes(restriction="TESTIS", column="Ensembl_Gene_ID")
    assert len(ids) == len(CTA_testis_restricted_gene_names())
