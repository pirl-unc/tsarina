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
    assert len(CTA_gene_names()) == 262


def test_gene_ids_nonempty():
    assert len(CTA_gene_ids()) == 262


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
        "rna_98_pct_filter",
        "passes_filters",
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


def test_evidence_uses_passes_filters_column_name():
    """``passes_filters`` is the canonical name. ``filtered`` is now
    surfaced as a backward-compat alias (tsarina#61) — both must be
    present and identical."""
    df = CTA_evidence()
    assert "passes_filters" in df.columns
    assert "filtered" in df.columns
    assert df["filtered"].equals(df["passes_filters"])


def test_xage1b_passes_relaxed_no_protein_threshold():
    df = CTA_evidence()
    row = df[df["Symbol"] == "XAGE1B"].iloc[0]
    assert bool(row["passes_filters"])
    assert bool(row["rna_98_pct_filter"])
    assert not bool(row["rna_99_pct_filter"])
    assert "XAGE1B" in CTA_gene_names()


def test_xage2_added_with_lung_safety_flag():
    """XAGE2 (CTpedia CT12.2) is added; it passes the relaxed 0.97 gate but
    carries a visible lung safety flag (real somatic expression). See
    tsarina#79."""
    df = CTA_evidence()
    rows = df[df["Symbol"] == "XAGE2"]
    assert len(rows) == 1
    row = rows.iloc[0]
    assert row["Ensembl_Gene_ID"] == "ENSG00000155622"
    assert bool(row["passes_filters"])
    assert not bool(row["never_expressed"])
    assert not bool(row["rna_98_pct_filter"])  # 0.977 < 0.98, passes only at 0.97
    assert "lung" in str(row["safety_flags"])
    assert "XAGE2" in CTA_gene_names()


def test_ct83_and_prm3_pass_at_097_threshold():
    """Lowering the no-protein gate 0.98 -> 0.97 admits CT83 (KK-LC-1) and
    PRM3, whose deflated reproductive fraction is in [0.97, 0.98)."""
    expressed = CTA_gene_names()
    assert {"CT83", "PRM3"} <= expressed


def test_xage5_rescued_into_expressed_set():
    """XAGE5 stays flagged never_expressed (HPA truth) but is rescued into the
    expressed set via MANUALLY_EXPRESSED_CTA. See tsarina#78."""
    df = CTA_evidence()
    row = df[df["Symbol"] == "XAGE5"].iloc[0]
    assert bool(row["never_expressed"])  # HPA-derived flag preserved
    assert "XAGE5" in CTA_gene_names()  # but rescued into expressed set
    assert "XAGE5" not in CTA_never_expressed_gene_names()


def test_evidence_does_not_bundle_runtime_ms_count_columns():
    df = CTA_evidence()
    count_columns = [
        column for column in df.columns if column.startswith("ms_") and "count" in column
    ]
    assert count_columns == []


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
    assert len(high) >= 120


def test_by_axes_gene_id_column():
    ids = CTA_by_axes(restriction="TESTIS", column="Ensembl_Gene_ID")
    assert len(ids) == len(CTA_testis_restricted_gene_names())


def test_no_protein_threshold_is_0_97():
    """The no-protein/Uncertain adaptive gate is the parameterized 0.97."""
    from tsarina.tissues import HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS

    assert HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS["Missing"] == 0.97
    assert HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS["Uncertain"] == 0.97


def test_never_expressed_column_matches_parameterized_floor():
    """The bundled ``never_expressed`` column equals the documented rule:
    no protein IHC data AND max RNA nTPM < HPA_EXPRESSION_FLOOR_NTPM."""
    from tsarina.tissues import HPA_EXPRESSION_FLOOR_NTPM

    df = CTA_evidence()
    no_protein = df["protein_reliability"].astype(str).str.lower().isin(["no data", "nan", ""])
    rule = no_protein & (df["rna_max_ntpm"] < HPA_EXPRESSION_FLOOR_NTPM)
    shipped = df["never_expressed"].astype(str).str.lower() == "true"
    assert (rule == shipped).all()
