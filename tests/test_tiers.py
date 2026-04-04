import pandas as pd

from tsarina import CTA_evidence, CTA_filtered_gene_names, CTA_gene_names
from tsarina.tiers import (
    CONFIDENCE_VALUES,
    MS_RESTRICTION_VALUES,
    RESTRICTION_VALUES,
    RNA_RESTRICTION_LEVELS,
    aggregate_gene_ms_safety,
    assign_all_axes,
    confidence_rank,
    ms_restriction_rank,
    restriction_rank,
)


def _axis_sets(axis: str):
    """Return dict mapping axis value -> set of gene symbols."""
    df = CTA_evidence()
    result = {}
    for val in df[axis].dropna().unique():
        if val:
            result[val] = set(df.loc[df[axis] == val, "Symbol"])
    return result


# ── Constants ──────────────────────────────────────────────────────────────


def test_restriction_values():
    assert RESTRICTION_VALUES == ["TESTIS", "PLACENTAL", "REPRODUCTIVE", "SOMATIC"]


def test_rna_restriction_levels():
    assert RNA_RESTRICTION_LEVELS == ["STRICT", "MODERATE", "PERMISSIVE", "LEAKY"]


def test_ms_restriction_values():
    assert len(MS_RESTRICTION_VALUES) == 5


def test_confidence_values():
    assert CONFIDENCE_VALUES == ["HIGH", "MODERATE", "LOW"]


# ── Rank functions ─────────────────────────────────────────────────────────


def test_restriction_rank_monotonic():
    ranks = [restriction_rank(v) for v in RESTRICTION_VALUES]
    assert ranks == sorted(ranks)


def test_ms_restriction_rank_monotonic():
    ranks = [ms_restriction_rank(v) for v in MS_RESTRICTION_VALUES]
    assert ranks == sorted(ranks)


def test_confidence_rank_monotonic():
    ranks = [confidence_rank(v) for v in CONFIDENCE_VALUES]
    assert ranks == sorted(ranks)


def test_empty_rank_is_worst():
    assert restriction_rank("") > restriction_rank(RESTRICTION_VALUES[-1])
    assert ms_restriction_rank("") > ms_restriction_rank(MS_RESTRICTION_VALUES[-1])


# ── Protein restriction ──────────────────────────────────────────────────


def test_protein_restriction_counts():
    sets = _axis_sets("protein_restriction")
    assert len(sets.get("TESTIS", set())) >= 150
    assert len(sets.get("PLACENTAL", set())) >= 3
    assert len(sets.get("REPRODUCTIVE", set())) >= 15


def test_protein_restriction_mutually_exclusive():
    sets = _axis_sets("protein_restriction")
    vals = list(sets.values())
    for i in range(len(vals)):
        for j in range(i + 1, len(vals)):
            assert vals[i] & vals[j] == set()


# ── RNA restriction ──────────────────────────────────────────────────────


def test_rna_restriction_counts():
    sets = _axis_sets("rna_restriction")
    assert len(sets.get("TESTIS", set())) >= 150
    assert len(sets.get("SOMATIC", set())) >= 100
    assert len(sets.get("PLACENTAL", set())) >= 3


def test_rna_restriction_level_counts():
    sets = _axis_sets("rna_restriction_level")
    assert len(sets.get("STRICT", set())) >= 180
    assert len(sets.get("MODERATE", set())) >= 50
    assert len(sets.get("PERMISSIVE", set())) >= 20
    assert len(sets.get("LEAKY", set())) >= 40


# ── Per-tissue RNA nTPMs ─────────────────────────────────────────────────


def test_rna_testis_ntpm_populated():
    df = CTA_evidence()
    assert "rna_testis_ntpm" in df.columns
    # Most CTA genes should have testis expression
    assert (df["rna_testis_ntpm"].astype(float) > 0).sum() > 200


def test_rna_somatic_columns_exist():
    df = CTA_evidence()
    for col in ["rna_max_somatic_tissue", "rna_max_somatic_ntpm", "rna_somatic_detected_count"]:
        assert col in df.columns


# ── Synthesized restriction ──────────────────────────────────────────────


def test_synthesized_restriction_covers_all_genes():
    """Every gene gets a synthesized restriction (no empties)."""
    df = CTA_evidence()
    assert (df["restriction"].fillna("") != "").all()


def test_synthesized_restriction_expressed_have_tissue():
    """Expressed (filtered + not never_expressed) genes all have a tissue."""
    df = CTA_evidence()
    filt = df["filtered"].astype(str).str.lower() == "true"
    ne = df["never_expressed"].astype(str).str.lower() == "true"
    expressed = df[filt & ~ne]
    assert (expressed["restriction"] != "NO_DATA").all()


def test_restriction_confidence_exists():
    df = CTA_evidence()
    assert "restriction_confidence" in df.columns
    sets = _axis_sets("restriction_confidence")
    assert "HIGH" in sets


# ── Known gene spot checks ───────────────────────────────────────────────


def test_magea1_testis():
    df = CTA_evidence()
    row = df[df["Symbol"] == "MAGEA1"].iloc[0]
    assert row["protein_restriction"] == "TESTIS"
    assert row["rna_restriction"] == "TESTIS"
    assert row["restriction"] == "TESTIS"
    assert row["restriction_confidence"] == "HIGH"


def test_prame_testis_but_permissive_rna():
    df = CTA_evidence()
    row = df[df["Symbol"] == "PRAME"].iloc[0]
    assert row["protein_restriction"] == "TESTIS"
    assert row["rna_restriction_level"] == "PERMISSIVE"
    assert row["restriction"] == "TESTIS"


def test_magea4_reproductive():
    """MAGEA4 has protein in placenta+testis → REPRODUCTIVE."""
    df = CTA_evidence()
    row = df[df["Symbol"] == "MAGEA4"].iloc[0]
    assert row["protein_restriction"] == "REPRODUCTIVE"


# ── Per-tissue protein flags ─────────────────────────────────────────────


def test_protein_flags_for_testis_gene():
    df = CTA_evidence()
    row = df[df["Symbol"] == "MAGEA1"].iloc[0]
    assert str(row["protein_testis"]) == "True"
    assert str(row["protein_ovary"]) == "False"
    assert str(row["protein_placenta"]) == "False"


# ── MS defaults ──────────────────────────────────────────────────────────


def test_ms_restriction_defaults_to_no_ms_data():
    df = CTA_evidence()
    assert (df["ms_restriction"] == "NO_MS_DATA").all()


# ── CSV columns ──────────────────────────────────────────────────────────


def test_csv_has_all_new_columns():
    df = CTA_evidence()
    expected = [
        "protein_restriction",
        "protein_testis",
        "protein_ovary",
        "protein_placenta",
        "rna_restriction",
        "rna_restriction_level",
        "rna_testis_ntpm",
        "rna_ovary_ntpm",
        "rna_placenta_ntpm",
        "rna_max_somatic_tissue",
        "rna_max_somatic_ntpm",
        "rna_somatic_detected_count",
        "ms_restriction",
        "restriction",
        "restriction_confidence",
    ]
    for col in expected:
        assert col in df.columns, f"Missing column: {col}"


def test_csv_has_no_old_columns():
    df = CTA_evidence()
    for col in ["cta_tier", "evidence", "ms_safety"]:
        assert col not in df.columns, f"Stale column present: {col}"


# ── Backward compatibility ────────────────────────────────────────────────


def test_gene_names_count_unchanged():
    assert len(CTA_gene_names()) == 257


def test_filtered_count_unchanged():
    assert len(CTA_filtered_gene_names()) == 278


# ── assign_all_axes runtime consistency ──────────────────────────────────


def test_assign_all_axes_matches_csv():
    df = CTA_evidence()
    recomputed = assign_all_axes(df)
    assert (
        recomputed["protein_restriction"].fillna("") == df["protein_restriction"].fillna("")
    ).all()
    assert (recomputed["rna_restriction"].fillna("") == df["rna_restriction"].fillna("")).all()
    assert (recomputed["restriction"].fillna("") == df["restriction"].fillna("")).all()


# ── MS safety aggregation (unit tests with mock data) ────────────────────


def test_ms_cancer_only():
    hits = pd.DataFrame(
        {
            "peptide": ["ABCDEFGH", "ABCDEFGH"],
            "src_cancer": [True, True],
            "src_healthy_tissue": [False, False],
            "src_healthy_reproductive": [False, False],
            "src_healthy_thymus": [False, False],
            "source_tissue": ["melanoma", "melanoma"],
        }
    )
    gene_map = pd.DataFrame({"peptide": ["ABCDEFGH"], "gene_name": ["GENE1"]})
    result = aggregate_gene_ms_safety(hits, gene_map)
    assert result.iloc[0]["ms_restriction"] == "CANCER_ONLY"


def test_ms_expected_tissue():
    hits = pd.DataFrame(
        {
            "peptide": ["ABCDEFGH"],
            "src_cancer": [False],
            "src_healthy_tissue": [False],
            "src_healthy_reproductive": [True],
            "src_healthy_thymus": [False],
            "source_tissue": ["testis"],
        }
    )
    gene_map = pd.DataFrame({"peptide": ["ABCDEFGH"], "gene_name": ["GENE1"]})
    result = aggregate_gene_ms_safety(hits, gene_map)
    assert result.iloc[0]["ms_restriction"] == "EXPECTED_TISSUE"


def test_ms_singleton_healthy():
    hits = pd.DataFrame(
        {
            "peptide": ["ABCDEFGH"],
            "src_cancer": [True],
            "src_healthy_tissue": [True],
            "src_healthy_reproductive": [False],
            "src_healthy_thymus": [False],
            "source_tissue": ["liver"],
        }
    )
    gene_map = pd.DataFrame({"peptide": ["ABCDEFGH"], "gene_name": ["GENE1"]})
    result = aggregate_gene_ms_safety(hits, gene_map)
    assert result.iloc[0]["ms_restriction"] == "SINGLETON_HEALTHY"


def test_ms_recurrent_healthy():
    hits = pd.DataFrame(
        {
            "peptide": ["PEP1", "PEP2", "PEP2"],
            "src_cancer": [True, True, False],
            "src_healthy_tissue": [True, False, True],
            "src_healthy_reproductive": [False, False, False],
            "src_healthy_thymus": [False, False, False],
            "source_tissue": ["liver", "melanoma", "kidney"],
        }
    )
    gene_map = pd.DataFrame({"peptide": ["PEP1", "PEP2"], "gene_name": ["GENE1", "GENE1"]})
    result = aggregate_gene_ms_safety(hits, gene_map)
    assert result.iloc[0]["ms_restriction"] == "RECURRENT_HEALTHY"


def test_ms_empty_inputs():
    result = aggregate_gene_ms_safety(pd.DataFrame(), pd.DataFrame())
    assert len(result) == 0
    assert "ms_restriction" in result.columns
    assert "ms_ebv_lcl_peptide_count" in result.columns
