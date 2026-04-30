from tsarina.mhc import (
    mhc_restriction_matches_any,
    normalize_mhc_restriction,
    normalize_mhc_restriction_set,
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
