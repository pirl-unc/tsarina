import pytest

from tsarina.mhc import (
    mhc_restriction_matches_any,
    normalize_mhc_restriction,
    normalize_mhc_restriction_set,
    parse_mhc,
    serotype_key,
    serotype_keys,
    split_mhc_restrictions,
    strip_hla_prefix,
)


def test_normalize_mhc_restriction_adds_hla_prefix():
    assert normalize_mhc_restriction("A*02:01") == "HLA-A*02:01"


def test_normalize_mhc_restriction_ignores_case():
    """A regression guard for the fix, not just the old passing case above:

    mhcgnomes' own case handling covers ``a2``, so no hand-rolled prefix
    trick is needed to accept it (the old hand-rolled ``HLA-`` force-prefix
    this replaced had regressed lowercase queries)."""
    assert normalize_mhc_restriction("a2") == normalize_mhc_restriction("A2")


@pytest.mark.parametrize(
    "value",
    ["SLA1*01:01", "BoLA1*001:01", "DLA88*001:01", "RT1A*01:01"],
)
def test_normalize_mhc_restriction_passes_through_non_human_designations(value):
    """Non-human MHC designations are left unchanged, not silently reparsed.

    mhcgnomes recognizes these directly (no forced HLA- prefix needed to
    reject them), so without a species constraint they would otherwise
    reformat into that species' own canonical string -- a pig SLA becoming
    e.g. 'SLA-1*01:01' -- which breaks the documented "canonical HLA
    restriction string" contract and any exact-match comparison against a
    corpus value in the original format."""
    assert normalize_mhc_restriction(value) == value


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


@pytest.mark.parametrize("value", [None, 123, 1.5])
def test_parse_mhc_returns_none_for_non_string_input(value):
    """A raw non-string value (e.g. a pandas NaN/None cell) must degrade to
    None like any other unparseable input, not raise AttributeError from
    calling .strip() on it."""
    assert parse_mhc(value) is None


def test_parse_mhc_returns_none_for_an_unrecognized_expect():
    """An ``expect`` outside {"", "serotype", "allele"} -- a typo, or a
    future caller passing the wrong casing -- must degrade to None like any
    other unsatisfiable request, not raise KeyError."""
    assert parse_mhc("A2", expect="Serotype") is None
    assert parse_mhc("A2", expect="nonsense") is None


def test_parse_mhc_species_constrains_to_that_species():
    assert parse_mhc("A*02:01", species="HLA") is not None
    assert parse_mhc("SLA1*01:01", species="HLA") is None
    assert parse_mhc("SLA1*01:01", species=None) is not None


def test_strip_hla_prefix_is_case_insensitive():
    assert strip_hla_prefix("HLA-A24") == "A24"
    assert strip_hla_prefix("hla-a24") == "a24"
    assert strip_hla_prefix("A24") == "A24"


def test_serotype_key_does_not_itself_resolve_a_split_to_its_broad_parent():
    """serotype_key keys one token to itself; it is serotype_keys (plural,
    reading a stored ``serotypes`` cell that lists both names) that makes a
    split serotype match its broad parent -- see test_filter_by_serotype_
    finds_split_serotype_members_of_a_broad_query for that behavior."""
    assert serotype_key("A2403") == "A2403"
    assert serotype_key("A2403") != serotype_key("A24")


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
