import pytest

from tsarina.mhc import (
    mhc_restriction_matches_any,
    normalize_mhc_restriction,
    normalize_mhc_restriction_set,
    parse_mhc,
    serotype_key,
    serotype_keys,
    split_mhc_restrictions,
)


def test_normalize_mhc_restriction_adds_hla_prefix():
    assert normalize_mhc_restriction("A*02:01") == "HLA-A*02:01"


def test_split_mhc_restrictions_normalizes_semicolon_joined_cell():
    assert split_mhc_restrictions("A*02:01; HLA-B*07:02") == (
        "HLA-A*02:01",
        "HLA-B*07:02",
    )


def test_mhc_restriction_matches_any_uses_normalized_values():
    wanted = normalize_mhc_restriction_set(["A*02:01"])
    assert mhc_restriction_matches_any("HLA-A*02:01;HLA-B*07:02", wanted)
    assert not mhc_restriction_matches_any("HLA-B*07:02", wanted)


def test_parse_mhc_requires_the_expected_reading():
    """A stated expectation wins over whatever mhcgnomes would prefer."""
    assert parse_mhc("A2", expect="serotype") is not None
    assert parse_mhc("A2", expect="allele") is None
    assert parse_mhc("A*02:01", expect="allele") is not None
    assert parse_mhc("A*02:01", expect="serotype") is None


def test_parse_mhc_is_cached():
    parse_mhc.cache_clear()
    parse_mhc("HLA-A*02:01")
    parse_mhc("HLA-A*02:01")
    assert parse_mhc.cache_info().hits == 1


def test_serotype_key_ignores_case_and_prefix():
    assert serotype_key("A2") == serotype_key("hla-a2") == serotype_key("HLA-A2")
    assert serotype_key("Bw4") == serotype_key("bw4")


@pytest.mark.parametrize("name", ["DR1B", "DR3A", "DR7A"])
def test_serotype_key_keeps_legacy_curated_names(name):
    """Existing indexed labels remain queryable across spelling variants."""
    for token in (name, name.lower(), f"HLA-{name}", f" hla-{name.lower()} "):
        assert serotype_key(token) == name


def test_serotype_key_rejects_what_is_not_a_serotype():
    for value in (
        "A*02:01",
        "HLA-A*02:01",
        "A0201",
        "hla-a0201",
        "DRB10401",
        "ABC123",
        "A999",
        "DR1C",
        "HLA class I",
        "nonsense",
        "",
        None,
    ):
        assert serotype_key(value) is None


def test_serotype_keys_splits_a_stored_cell():
    assert serotype_keys("HLA-A24;HLA-Bw4") == {"A24", "BW4"}
    assert serotype_keys("") == frozenset()
    assert serotype_keys(None) == frozenset()
