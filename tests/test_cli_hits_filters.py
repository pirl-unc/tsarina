"""Unit tests for tsarina.cli_hits filter helpers."""

import pandas as pd

from tsarina.cli_hits import (
    _apply_min_resolution,
    _filter_by_allele,
    _filter_by_serotype,
)


def _hits(mhcs: list[str]) -> pd.DataFrame:
    return pd.DataFrame({"peptide": ["P"] * len(mhcs), "mhc_restriction": mhcs})


def test_filter_by_allele_exact_match():
    df = _hits(["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02"])
    out = _filter_by_allele(df, ["HLA-A*24:02"])
    assert out["mhc_restriction"].tolist() == ["HLA-A*24:02"]


def test_filter_by_allele_empty_passthrough():
    df = _hits(["HLA-A*02:01"])
    assert _filter_by_allele(df, []).equals(df)


def test_filter_by_serotype_accepts_A24_for_A_star_24_02():
    """hitlist#44: canonical serotype for A*24:02 is mis-reported as Bw4.
    Our filter must still find A*24:02 when the user types A24."""
    df = _hits(["HLA-A*24:02", "HLA-A*02:01", "HLA-B*07:02"])
    out = _filter_by_serotype(df, ["A24"])
    assert out["mhc_restriction"].tolist() == ["HLA-A*24:02"]


def test_filter_by_serotype_accepts_HLA_prefix():
    df = _hits(["HLA-A*24:02", "HLA-B*07:02"])
    out = _filter_by_serotype(df, ["HLA-A24"])
    assert out["mhc_restriction"].tolist() == ["HLA-A*24:02"]


def test_filter_by_serotype_matches_serological_restriction():
    df = _hits(["HLA-A24", "HLA-B*07:02"])
    out = _filter_by_serotype(df, ["A24"])
    assert out["mhc_restriction"].tolist() == ["HLA-A24"]


def test_filter_by_serotype_a2_via_canonical_mapping():
    df = _hits(["HLA-A*02:01", "HLA-B*07:02"])
    out = _filter_by_serotype(df, ["A2"])
    assert out["mhc_restriction"].tolist() == ["HLA-A*02:01"]


def test_apply_min_resolution_drops_coarser_alleles():
    df = _hits(["HLA-A*02:01", "HLA-A2", "HLA class I"])
    out = _apply_min_resolution(df, "four_digit")
    assert out["mhc_restriction"].tolist() == ["HLA-A*02:01"]


def test_apply_min_resolution_passthrough_when_none():
    df = _hits(["HLA-A*02:01", "HLA-A2"])
    assert _apply_min_resolution(df, None).equals(df)
