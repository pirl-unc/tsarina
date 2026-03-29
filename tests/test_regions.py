import pandas as pd

from tsarina.regions import (
    REGION_POPULATIONS,
    REGION_PRIORITY_ROWS,
    region_allele_frequencies,
    region_names,
)


def test_region_priority_rows_nonempty():
    assert len(REGION_PRIORITY_ROWS) > 50


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


def test_all_alleles_start_with_hla():
    for row in REGION_PRIORITY_ROWS:
        assert row["allele"].startswith("HLA-"), f"{row['allele']} doesn't start with HLA-"


def test_frequencies_are_valid():
    for row in REGION_PRIORITY_ROWS:
        freq = row["frequency"]
        if freq is not None:
            assert 0 < freq < 1, f"Invalid frequency {freq} for {row['allele']} in {row['region']}"


def test_all_loci_are_abc():
    loci = {row["locus"] for row in REGION_PRIORITY_ROWS}
    assert loci == {"A", "B", "C"}
