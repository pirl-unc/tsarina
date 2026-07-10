import oncoref

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
    CTA_relaxed_reproductive_gene_ids,
    CTA_relaxed_reproductive_gene_names,
    CTA_testis_restricted_gene_ids,
    CTA_testis_restricted_gene_names,
    CTA_unfiltered_gene_ids,
    CTA_unfiltered_gene_names,
)


def test_gene_names_nonempty():
    assert len(CTA_gene_names()) == 293


def test_gene_ids_nonempty():
    assert len(CTA_gene_ids()) == 293


def test_default_cta_set_matches_oncoref():
    assert CTA_gene_names() == oncoref.cta_gene_names()
    assert CTA_gene_ids() == oncoref.cta_gene_ids()


def test_canonical_filtered_cta_set_matches_oncoref():
    assert CTA_filtered_gene_names() == oncoref.cta_filtered_gene_names()
    assert CTA_filtered_gene_ids() == oncoref.cta_filtered_gene_ids()


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


def test_relaxed_reproductive_tier_is_opt_in_and_disjoint():
    """tsarina#119: the relaxed tier is RNA-only failing genes that remain
    reproductive-dominant; disjoint from the default sets."""
    relaxed = CTA_relaxed_reproductive_gene_names()  # default min_deflated_frac=0.80
    # Opt-in: never overlaps the default expressed/filtered sets.
    assert relaxed & CTA_filtered_gene_names() == set()
    assert relaxed & CTA_gene_names() == set()
    # Every member is an excluded (filter-failing) candidate.
    assert relaxed <= CTA_excluded_gene_names()
    # Placental ERV envelopes (syncytin-1/2) are the motivating inclusions.
    assert {"ERVW-1", "ERVFRD-1"} <= relaxed
    # A higher floor is a strict subset (fewer, more reproductive-dominant).
    assert CTA_relaxed_reproductive_gene_names(0.95) <= relaxed
    # ERVV-1 (deflated frac 0.954) is in at 0.90 but the syncytins (~0.85) are not.
    at_90 = CTA_relaxed_reproductive_gene_names(0.90)
    assert "ERVV-1" in at_90
    assert "ERVW-1" not in at_90
    # IDs accessor agrees by row.
    assert len(CTA_relaxed_reproductive_gene_ids()) == len(relaxed)


def test_relaxed_reproductive_excludes_somatic_protein_genes():
    """Genes that fail on somatic *protein* (not RNA fraction) must be excluded
    -- the tier is restricted to RNA-only genes so the relaxed RNA threshold is
    meaningful."""
    df = CTA_evidence()
    relaxed = CTA_relaxed_reproductive_gene_names()
    for sym in relaxed:
        rel = str(df.loc[df["Symbol"] == sym, "protein_reliability"].iloc[0]).strip().lower()
        assert rel in ("no data", "nan", ""), f"{sym} has protein data but is in relaxed tier"


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
        "specificity_status",
        "specificity_action",
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


def test_csh1_and_h1_6_are_excluded_candidates_with_supporting_evidence():
    df = CTA_evidence().set_index("Symbol")
    default = CTA_gene_names()
    filtered = CTA_filtered_gene_names()
    excluded = CTA_excluded_gene_names()

    assert {"CSH1", "H1-6"}.isdisjoint(default)
    assert {"CSH1", "H1-6"}.isdisjoint(filtered)
    assert {"CSH1", "H1-6"} <= excluded

    csh1 = df.loc["CSH1"]
    assert csh1["Ensembl_Gene_ID"] == "ENSG00000136488"
    assert not bool(csh1["passes_filters"])
    assert not bool(csh1["rna_99_pct_filter"])
    assert csh1["rna_deflated_reproductive_frac"] == 0.9898
    assert csh1["rna_max_somatic_tissue"] == "smooth muscle"
    assert csh1["rna_max_somatic_ntpm"] == 38.3
    assert csh1["rna_lung_max_ntpm"] == 10.4
    assert csh1["restriction"] == "SOMATIC"
    assert "lung" in str(csh1["safety_flags"])
    assert csh1["specificity_status"] == "excluded_normal_expression"
    assert csh1["specificity_action"] == "exclude_default"

    h1_6 = df.loc["H1-6"]
    assert h1_6["Ensembl_Gene_ID"] == "ENSG00000187475"
    assert bool(h1_6["passes_filters"])  # raw HPA row retained for audit.
    assert h1_6["specificity_status"] == "excluded_structural_gene"
    assert h1_6["specificity_action"] == "exclude_default"
    assert h1_6["ms_restriction"] == "RECURRENT_HEALTHY"
    tissues = {t.strip() for t in str(h1_6["ms_healthy_somatic_tissues"]).split(";") if t}
    assert len(tissues) == 15
    assert {"lung", "liver", "blood", "bone marrow"} <= tissues
    pmids = {p for p in str(h1_6["ms_pmids"]).split(";") if p}
    assert len(pmids) == 48


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


def test_mageb6_added_clean_testis():
    """MAGEB6 (CT3.2) — dual-corroborated CTA from #79, the only one of the 9
    candidates that passes the HPA reproductive-restriction filter (testis-only
    3.0 nTPM, no protein)."""
    df = CTA_evidence()
    rows = df[df["Symbol"] == "MAGEB6"]
    assert len(rows) == 1
    row = rows.iloc[0]
    assert row["Ensembl_Gene_ID"] == "ENSG00000176746"
    assert bool(row["passes_filters"])
    assert not bool(row["never_expressed"])
    assert row["restriction"] == "TESTIS"
    assert "MAGEB6" in CTA_gene_names()


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


def test_never_expressed_column_matches_oncoref_evidence():
    """oncoref owns raw HPA CTA evidence, including ``never_expressed``."""
    df = CTA_evidence().copy()
    df["_gid"] = df["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
    ref = oncoref.cta_evidence().copy()
    ref["_gid"] = ref["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
    merged = df.merge(ref[["_gid", "never_expressed"]], on="_gid", suffixes=("", "_oncoref"))
    left = merged["never_expressed"].astype(str).str.lower()
    right = merged["never_expressed_oncoref"].astype(str).str.lower()
    assert left.equals(right)


def test_cta_symbol_for_alias_resolves_common_names():
    """Common antigen names + NCBI synonyms resolve to official symbols (#77)."""
    from tsarina import cta_symbol_for_alias

    # The issue's headline example: NY-ESO-1 / ESO1 -> CTAG1B (not the paralog).
    assert cta_symbol_for_alias("NY-ESO-1") == "CTAG1B"
    assert cta_symbol_for_alias("ESO1") == "CTAG1B"
    # Punctuation/case-insensitive.
    assert cta_symbol_for_alias("ny eso 1") == "CTAG1B"
    # NCBI + curated synonyms across the table.
    assert cta_symbol_for_alias("CT12.2") == "XAGE2"
    assert cta_symbol_for_alias("XAGE1C") == "XAGE1B"
    assert cta_symbol_for_alias("MAGE4") == "MAGEA4"
    # Official symbols resolve to themselves.
    assert cta_symbol_for_alias("MAGEA4") == "MAGEA4"
    assert cta_symbol_for_alias("CTAG1A") == "CTAG1A"
    # Unknown names return None.
    assert cta_symbol_for_alias("not-a-gene") is None


def test_aliases_backfilled_for_most_genes():
    """The Aliases column is populated broadly, not ~96% empty (#77)."""
    df = CTA_evidence()
    with_aliases = (df["Aliases"].fillna("").astype(str).str.len() > 0).sum()
    assert with_aliases >= 300


def test_non_cta_conserved_genes_excluded_from_universe():
    """Conserved/multicopy non-CTAs (core histones, alpha-tubulins) are dropped
    from the CTA universe (tsarina#92). H1-6 is kept only as a tsarina evidence
    candidate and is excluded from canonical oncoref CTA defaults
    (tsarina#141). hCG-beta (CGB8) is no longer excluded by hand — it now enters
    the candidate universe with the rest of the placental-antigen families and
    is judged by the reproductive-restriction filter (tsarina#111)."""
    from tsarina import CTA_evidence, CTA_unfiltered_gene_names

    universe = CTA_unfiltered_gene_names()
    evidence_symbols = set(CTA_evidence()["Symbol"])
    for gene in ("H4C6", "H2BC1", "H2BC3", "H1-1", "TUBA3C", "TUBA3E"):
        assert gene not in universe, f"{gene} should be excluded from the universe"
        assert gene not in evidence_symbols, f"{gene} should not be a CTA_evidence row"
    # The testis-specific linker histone H1t (H1-6) is retained only for
    # evidence auditability, not canonical CTA membership.
    assert "H1-6" in universe
    assert "H1-6" not in CTA_filtered_gene_names()
    assert "H1-6" not in CTA_gene_names()
    # hCG-beta now enters the universe as a placental candidate (tsarina#111).
    assert "CGB8" in universe


def test_non_cta_exclusion_preserves_evidence_unfiltered_invariant():
    from tsarina import CTA_evidence, CTA_unfiltered_gene_names

    assert len(CTA_evidence()) == len(CTA_unfiltered_gene_names())


def test_paralog_copies_added_to_universe():
    """Near-identical CTA paralog copies are in the universe so they don't
    pollute downstream non-CTA negative sets (tsarina#93). Tagged paralog:<sibling>."""
    from tsarina import CTA_evidence, CTA_unfiltered_gene_names

    universe = CTA_unfiltered_gene_names()
    for gene in ("CT47A8", "CT47A9", "CT47A10", "CT45A8", "CT45A9", "GAGE12D"):
        assert gene in universe, f"{gene} should be in the universe"
    df = CTA_evidence().set_index("Symbol")
    assert df.loc["GAGE12D", "source_databases"].startswith("paralog:")
    # MAGEA2B/SSX4B/CT45A5 previously carried their siblings' Ensembl gene IDs,
    # so they collided with the existing rows and were silently dropped. With the
    # corrected IDs (tsarina#111) they are now distinct universe rows of their own
    # — identical-protein paralogs that must be present so their peptides don't
    # leak into the non-CTA negative set.
    for gene, sibling in (("MAGEA2B", "MAGEA2"), ("SSX4B", "SSX4"), ("CT45A5", "CT45A1")):
        assert gene in universe, f"{gene} should be in the universe"
        assert df.loc[gene, "source_databases"] == f"paralog:{sibling}"


def test_rna_97_pct_filter_column_matches_active_threshold():
    """The bundled table carries rna_97_pct_filter for the active no-protein
    gate (0.97), consistent with the deflated reproductive fraction."""
    df = CTA_evidence()
    assert "rna_97_pct_filter" in df.columns
    rule = df["rna_deflated_reproductive_frac"] >= 0.97
    col = df["rna_97_pct_filter"].astype(str).str.lower() == "true"
    assert (rule == col).all()
    # Ordered between the 95 and 98 columns.
    cols = list(df.columns)
    assert (
        cols.index("rna_95_pct_filter")
        < cols.index("rna_97_pct_filter")
        < cols.index("rna_98_pct_filter")
    )
