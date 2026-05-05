import pandas as pd

from tsarina.regions import (
    GLOBAL_ALLELE_FREQUENCY_ROWS,
    REGION_POPULATIONS,
    REGION_PRIORITY_ROWS,
    allele_frequency_audit,
    allele_frequency_rows,
    global_allele_frequencies,
    region_allele_frequencies,
    region_names,
)


def test_region_priority_rows_nonempty():
    assert len(REGION_PRIORITY_ROWS) > 50


def test_global_frequency_rows_nonempty():
    assert len(GLOBAL_ALLELE_FREQUENCY_ROWS) > 10


def test_all_rows_have_required_keys():
    required = {
        "region",
        "proxy",
        "locus",
        "allele",
        "frequency",
        "resolution",
        "source_label",
        "source_url",
        "note",
    }
    for i, row in enumerate(REGION_PRIORITY_ROWS):
        assert required.issubset(row.keys()), f"Row {i} missing keys: {required - row.keys()}"


def test_global_rows_have_required_keys():
    required = {
        "region",
        "proxy",
        "locus",
        "allele",
        "frequency",
        "resolution",
        "source_label",
        "source_url",
        "note",
    }
    for i, row in enumerate(GLOBAL_ALLELE_FREQUENCY_ROWS):
        assert required.issubset(row.keys()), f"Row {i} missing keys: {required - row.keys()}"


def test_region_populations_has_seven_regions():
    assert len(REGION_POPULATIONS) == 7


def test_region_names_matches_populations():
    names = set(region_names())
    pops = set(REGION_POPULATIONS.keys())
    assert names == pops, f"Mismatch: {names.symmetric_difference(pops)}"


def test_region_allele_frequencies_returns_dataframe():
    df = region_allele_frequencies()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == len(REGION_PRIORITY_ROWS)
    assert "allele" in df.columns
    assert "frequency" in df.columns
    assert "region" in df.columns


def test_global_allele_frequencies_returns_dataframe():
    df = global_allele_frequencies()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == len(GLOBAL_ALLELE_FREQUENCY_ROWS)
    assert "allele" in df.columns
    assert "frequency" in df.columns
    assert set(df["region"]) == {"Global"}


def test_allele_frequency_rows_returns_compatible_schema():
    df = allele_frequency_rows()
    assert isinstance(df, pd.DataFrame)
    assert {"regional_proxy", "published_global"} <= set(df["frequency_scope"])
    expected = {
        "frequency_scope",
        "region",
        "proxy",
        "locus",
        "allele",
        "frequency",
        "resolution",
        "source_label",
        "source_url",
        "note",
    }
    assert expected <= set(df.columns)


def test_all_alleles_start_with_hla():
    for row in REGION_PRIORITY_ROWS:
        assert row["allele"].startswith("HLA-"), f"{row['allele']} doesn't start with HLA-"
    for row in GLOBAL_ALLELE_FREQUENCY_ROWS:
        assert row["allele"].startswith("HLA-"), f"{row['allele']} doesn't start with HLA-"


def test_frequencies_are_valid():
    for row in REGION_PRIORITY_ROWS:
        freq = row["frequency"]
        if freq is not None:
            assert 0 < freq < 1, f"Invalid frequency {freq} for {row['allele']} in {row['region']}"
    for row in GLOBAL_ALLELE_FREQUENCY_ROWS:
        freq = row["frequency"]
        assert 0 < freq < 1, f"Invalid frequency {freq} for {row['allele']}"


def test_all_loci_are_abc():
    loci = {row["locus"] for row in REGION_PRIORITY_ROWS}
    assert loci == {"A", "B", "C"}
    global_loci = {row["locus"] for row in GLOBAL_ALLELE_FREQUENCY_ROWS}
    assert global_loci == {"A", "B", "C"}


def test_global53_panel_has_frequency_support():
    from tsarina.alleles import get_panel
    from tsarina.spanning import _weighted_allele_frequencies

    alleles = get_panel("global53_abc")
    weighted = _weighted_allele_frequencies(alleles)
    missing = [allele for allele in alleles if weighted.get(allele, 0.0) <= 0.0]
    assert missing == []


def test_global53_panel_has_auditable_subpop_and_global_frequencies():
    from tsarina.alleles import get_panel

    alleles = get_panel("global53_abc")
    audit = allele_frequency_audit(alleles)

    assert audit["allele"].tolist() == alleles
    assert audit["published_global_frequency"].notna().all()
    assert (audit["published_global_frequency"] > 0.0).all()
    assert (audit["coverage_frequency"] > 0.0).all()
    assert set(audit["coverage_frequency_source"]) == {"regional_weighted", "published_global"}
    assert (
        audit.loc[audit["coverage_frequency_source"] == "published_global", "numeric_region_count"]
        == 0
    ).all()
    assert (
        audit.loc[audit["coverage_frequency_source"] == "regional_weighted", "numeric_region_count"]
        > 0
    ).all()


def test_frequency_audit_keeps_global_average_separate_from_coverage_frequency():
    audit = allele_frequency_audit(["HLA-A*02:01", "HLA-A*23:01", "HLA-C*04:03"])
    by_allele = audit.set_index("allele")

    assert by_allele.loc["HLA-A*02:01", "coverage_frequency_source"] == "regional_weighted"
    assert (
        by_allele.loc["HLA-A*02:01", "published_global_frequency"]
        != by_allele.loc["HLA-A*02:01", "coverage_frequency"]
    )
    assert by_allele.loc["HLA-A*23:01", "coverage_frequency_source"] == "published_global"
    assert by_allele.loc["HLA-C*04:03", "global_source_label"] == (
        "Sarkizova HLA-C global allele frequency"
    )
