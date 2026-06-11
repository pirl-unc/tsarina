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
    synthesize_restriction,
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
    assert len(MS_RESTRICTION_VALUES) == 6


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
    assert len(sets.get("MODERATE", set())) >= 30
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
    """Expressed (passes_filters + not never_expressed) genes all have a tissue."""
    df = CTA_evidence()
    filt = df["passes_filters"].astype(str).str.lower() == "true"
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
    assert row["restriction_confidence"] == "MODERATE"


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


def test_ms_restriction_has_real_data():
    """Bundled CSV includes MS restriction from IEDB/CEDAR scan."""
    df = CTA_evidence()
    values = df["ms_restriction"].value_counts()
    assert "CANCER_ONLY" in values.index
    assert "NO_MS_DATA" in values.index
    assert values["CANCER_ONLY"] >= 100


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
    """Guard against pre-rename column names creeping back in.

    ``filtered`` is deliberately retained as a backward-compat alias of
    ``passes_filters`` (tsarina#61), so it's exempted here even though
    the rest of the list still represents stale columns."""
    df = CTA_evidence()
    for col in ["cta_tier", "evidence", "ms_safety"]:
        assert col not in df.columns, f"Stale column present: {col}"


def test_csv_has_no_runtime_ms_count_columns():
    df = CTA_evidence()
    count_columns = [
        column for column in df.columns if column.startswith("ms_") and "count" in column
    ]
    assert count_columns == []


# ── Backward compatibility ────────────────────────────────────────────────


def test_gene_names_count_unchanged():
    assert len(CTA_gene_names()) == 272


def test_filtered_count():
    # 297 = 286 + 11 reproductive-restriction-filter passers from the tsarina#111
    # batch (placental-antigen families plus the corrected MAGEA2B/SSX4B
    # identical-protein paralogs, which previously carried their siblings' gene
    # IDs and were silently dropped). Somatically-leaky members (most CGB, PSG4/7,
    # GH2, and the distinct-protein CT45A5) land in the universe as excluded
    # candidates, not counted here.
    assert len(CTA_filtered_gene_names()) == 297


def test_gage10_added_gage12b_excluded():
    """Pin the two GAGE-audit decisions by identity (not just count; tsarina#108).

    GAGE10 is a full-length GAGE paralog added to the universe; it passes the
    reproductive-restriction filter but is never_expressed (1.6 < 2.0 floor), so
    it is in the filtered set but not the expressed set. GAGE12B is deliberately
    excluded -- its Ensembl ID is a degenerate 117-bp fragment, not a GAGE
    protein -- so it must be in neither.
    """
    filtered = CTA_filtered_gene_names()
    expressed = CTA_gene_names()
    assert "GAGE10" in filtered
    assert "GAGE10" not in expressed
    assert "GAGE12B" not in filtered
    assert "GAGE12B" not in expressed


# ── assign_all_axes runtime consistency ──────────────────────────────────


def test_assign_all_axes_matches_csv():
    df = CTA_evidence()
    recomputed = assign_all_axes(df)
    assert (
        recomputed["protein_restriction"].fillna("") == df["protein_restriction"].fillna("")
    ).all()
    assert (recomputed["rna_restriction"].fillna("") == df["rna_restriction"].fillna("")).all()
    assert (recomputed["restriction"].fillna("") == df["restriction"].fillna("")).all()
    assert (
        recomputed["restriction_confidence"].fillna("") == df["restriction_confidence"].fillna("")
    ).all()


def test_confidence_capped_for_subfloor_rna_only():
    """tsarina#114: a restriction resting only on RNA below the expression floor
    (no protein, no MS) is capped at MODERATE — the scorer otherwise grants the
    same STRICT credit to near-noise RNA as to robust expression."""
    base = {
        "protein_restriction": "NO_DATA",
        "protein_reliability": "no data",
        "rna_restriction": "TESTIS",
        "rna_restriction_level": "STRICT",
        "ms_restriction": "NO_MS_DATA",
    }
    # Sub-floor RNA only -> capped to MODERATE (would be HIGH without the cap).
    assert synthesize_restriction(pd.Series({**base, "rna_max_ntpm": 1.6})) == (
        "TESTIS",
        "MODERATE",
    )
    # Same gene expressed above the 2.0 floor keeps HIGH.
    assert synthesize_restriction(pd.Series({**base, "rna_max_ntpm": 5.0})) == ("TESTIS", "HIGH")
    # MS evidence means real expression despite low RNA -> not capped.
    ms_row = {**base, "rna_max_ntpm": 1.6, "ms_restriction": "CANCER_ONLY"}
    assert synthesize_restriction(pd.Series(ms_row))[1] == "HIGH"
    # Protein evidence (IHC) likewise means real expression -> not capped even
    # at sub-floor RNA.
    prot_row = {
        **base,
        "rna_max_ntpm": 1.0,
        "protein_restriction": "TESTIS",
        "protein_reliability": "Supported",
    }
    assert synthesize_restriction(pd.Series(prot_row))[1] == "HIGH"
    # Boundary: exactly at the floor is "expressed" (cap is strict <), keeps HIGH.
    assert synthesize_restriction(pd.Series({**base, "rna_max_ntpm": 2.0})) == ("TESTIS", "HIGH")


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


def test_ms_unclassified_for_source_free_evidence():
    # A peptide with MS evidence but no classifiable source flag (e.g. cell-line
    # only, which the classifier ignores) must be UNCLASSIFIED_MS, *not*
    # NO_MS_DATA -- the gene does have MS data, just none we can place.
    hits = pd.DataFrame(
        {
            "peptide": ["ABCDEFGH"],
            "src_cancer": [False],
            "src_healthy_tissue": [False],
            "src_healthy_reproductive": [False],
            "src_healthy_thymus": [False],
            "source_tissue": ["cell_line"],
        }
    )
    gene_map = pd.DataFrame({"peptide": ["ABCDEFGH"], "gene_name": ["GENE1"]})
    result = aggregate_gene_ms_safety(hits, gene_map)
    assert result.iloc[0]["ms_peptide_count"] > 0
    assert result.iloc[0]["ms_restriction"] == "UNCLASSIFIED_MS"


def test_ms_cta_exclusive_counts_are_separate_from_all_cta_counts():
    hits = pd.DataFrame(
        {
            "peptide": ["EXCLUSIVE", "SHARED", "EXCLUSIVE", "REPRO", "THYMUS"],
            "src_cancer": [True, True, False, False, False],
            "src_healthy_tissue": [False, False, True, False, False],
            "src_healthy_reproductive": [False, False, False, True, False],
            "src_healthy_thymus": [False, False, False, False, True],
            "source_tissue": ["melanoma", "melanoma", "liver", "testis", "thymus"],
        }
    )
    gene_map = pd.DataFrame(
        {
            "peptide": ["EXCLUSIVE", "SHARED", "REPRO", "THYMUS"],
            "gene_name": ["GENE1", "GENE1", "GENE1", "GENE1"],
        }
    )
    exclusive_gene_map = pd.DataFrame(
        {
            "peptide": ["EXCLUSIVE", "REPRO", "THYMUS"],
            "gene_name": ["GENE1", "GENE1", "GENE1"],
        }
    )

    result = aggregate_gene_ms_safety(
        hits,
        gene_map,
        exclusive_peptide_gene_map=exclusive_gene_map,
    )
    row = result.iloc[0]

    assert row["ms_peptide_count"] == 4
    assert row["ms_cancer_peptide_count"] == 2
    assert row["ms_cta_exclusive_peptide_count"] == 3
    assert row["ms_cta_exclusive_cancer_peptide_count"] == 1
    assert row["ms_cta_exclusive_healthy_somatic_peptide_count"] == 1
    assert row["ms_cta_exclusive_healthy_reproductive_peptide_count"] == 1
    assert row["ms_cta_exclusive_healthy_thymus_peptide_count"] == 1


def test_ms_empty_inputs():
    result = aggregate_gene_ms_safety(pd.DataFrame(), pd.DataFrame())
    assert len(result) == 0
    assert "ms_restriction" in result.columns
    assert "ms_ebv_lcl_peptide_count" in result.columns
    assert "ms_cta_exclusive_cancer_peptide_count" in result.columns


def test_vital_organ_vocabulary_is_publicly_exported():
    """The vital-organ vocabulary is consumable from the package top level so
    downstream tools (e.g. vaxrank safety scoring) share it (vaxrank#303)."""
    import tsarina
    from tsarina.tiers import VITAL_TISSUE_MS_NAMES

    assert {"brain", "heart", "lung", "liver", "pancreas"} <= set(tsarina.SAFETY_TISSUE_GROUPS)
    assert "heart" in tsarina.VITAL_TISSUE_MS_NAMES
    assert tsarina.VITAL_TISSUE_MS_NAMES is VITAL_TISSUE_MS_NAMES
    # spanning's internal alias points at the same shared vocabulary
    from tsarina.spanning import _VITAL_TISSUE_MS_NAMES

    assert _VITAL_TISSUE_MS_NAMES is VITAL_TISSUE_MS_NAMES
