"""Unit tests for tsarina.cli_hits filter helpers."""

import pandas as pd
import pytest

from tsarina.cli_hits import (
    _apply_min_resolution,
    _filter_by_allele,
    _filter_by_serotype,
)


def _hits(mhcs: list[str]) -> pd.DataFrame:
    """Frames annotated exactly as hitlist annotates an observations row.

    ``--serotype`` and ``--min-resolution`` read the ``serotypes`` and
    ``allele_resolution`` columns that hitlist writes alongside
    ``mhc_restriction``, so the fixture derives all three from hitlist's own
    resolver rather than restating them.  A change in how hitlist annotates
    these restrictions therefore shows up here.
    """
    from hitlist.curation import resolve_mhc_annotation

    rows = [resolve_mhc_annotation(mhc).as_record_fields() for mhc in mhcs]
    frame = pd.DataFrame(rows)
    frame.insert(0, "peptide", ["P"] * len(mhcs))
    return frame


def test_filter_by_allele_exact_match():
    df = _hits(["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02"])
    out = _filter_by_allele(df, ["HLA-A*24:02"])
    assert out["mhc_restriction"].tolist() == ["HLA-A*24:02"]


def test_filter_by_allele_normalizes_user_input():
    df = _hits(["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02"])
    out = _filter_by_allele(df, ["A*02:01"])
    assert out["mhc_restriction"].tolist() == ["HLA-A*02:01"]


def test_filter_by_allele_matches_semicolon_joined_restrictions():
    df = _hits(["HLA-A*02:01;HLA-B*07:02", "HLA-A*24:02"])
    out = _filter_by_allele(df, ["B*07:02"])
    assert out["mhc_restriction"].tolist() == ["HLA-A*02:01;HLA-B*07:02"]


def test_filter_by_allele_empty_passthrough():
    df = _hits(["HLA-A*02:01"])
    assert _filter_by_allele(df, []).equals(df)


def test_filter_by_serotype_accepts_A24_for_A_star_24_02():
    """A*24:02 carries both HLA-A24 and the public epitope HLA-Bw4.

    hitlist#44 (canonical serotype mis-reported as Bw4) is fixed upstream and
    ``serotypes`` lists every membership, so the locus query still selects it.
    """
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


def test_filter_by_serotype_matches_public_epitope_across_loci():
    """A public epitope is an orthogonal serotype axis, not a locus label.

    Bw4 is carried by subsets of both A- and B-locus alleles, so it selects
    rows whose locus serotype is A23 or A24.
    """
    df = _hits(["HLA-A*23:01", "HLA-A*24:02", "HLA-A*02:01", "HLA-B*07:02"])
    out = _filter_by_serotype(df, ["Bw4"])
    assert out["mhc_restriction"].tolist() == ["HLA-A*23:01", "HLA-A*24:02"]


def test_filter_by_serotype_excludes_donor_sets():
    """A donor bag makes the serotype a candidate, not the restriction.

    The sample carried an A2 molecule, but which of its alleles presented the
    peptide is unknown, so this is not an A2-restricted observation. tsarina
    only credits a donor bag to one allele after deconvolution, and
    ``--min-resolution donor_set`` is how you look at these rows.
    """
    df = _hits(["HLA-A*01:01;HLA-A*02:01;HLA-B*44:03", "HLA-A*02:01"])
    out = _filter_by_serotype(df, ["A2"])
    assert out["mhc_restriction"].tolist() == ["HLA-A*02:01"]


def test_filter_by_serotype_ignores_query_case():
    df = _hits(["HLA-A*24:02", "HLA-A24", "HLA-B*07:02"])
    for query in ("A24", "a24", "HLA-A24", "hla-a24"):
        out = _filter_by_serotype(df, [query])
        assert out["mhc_restriction"].tolist() == ["HLA-A*24:02", "HLA-A24"], query


def test_filter_by_serotype_finds_split_serotype_members_of_a_broad_query():
    """A*24:03 is listed under the split A2403, and a split implies its parent.

    mhcgnomes' own ``A24`` member list omits it, so expanding the query instead
    of reading the allele's serotypes used to miss this row.
    """
    df = _hits(["HLA-A*24:03", "HLA-A*02:01"])
    assert _filter_by_serotype(df, ["A24"])["mhc_restriction"].tolist() == ["HLA-A*24:03"]
    assert _filter_by_serotype(df, ["A2403"])["mhc_restriction"].tolist() == ["HLA-A*24:03"]


def test_filter_by_serotype_rejects_a_query_it_cannot_read():
    """Silently returning every row, or none, would both be wrong answers."""
    df = _hits(["HLA-A*02:01"])
    with pytest.raises(ValueError, match="could not read"):
        _filter_by_serotype(df, ["HLA-A*02:01"])
    with pytest.raises(ValueError, match="could not read"):
        _filter_by_serotype(df, ["nonsense"])


def test_filter_by_serotype_accepts_multiple_queries():
    df = _hits(["HLA-A*24:02", "HLA-A*02:01", "HLA-B*07:02"])
    out = _filter_by_serotype(df, ["A24", "A2"])
    assert out["mhc_restriction"].tolist() == ["HLA-A*24:02", "HLA-A*02:01"]


def test_filter_by_serotype_empty_passthrough():
    df = _hits(["HLA-A*02:01"])
    assert _filter_by_serotype(df, []).equals(df)
    assert _filter_by_serotype(df, ["  "]).equals(df)


def test_filter_by_serotype_requires_annotated_index():
    legacy = pd.DataFrame({"peptide": ["P"], "mhc_restriction": ["HLA-A*24:02"]})
    with pytest.raises(ValueError, match="serotypes"):
        _filter_by_serotype(legacy, ["A24"])


def test_apply_min_resolution_drops_coarser_alleles():
    df = _hits(["HLA-A*02:01", "HLA-A2", "HLA class I"])
    out = _apply_min_resolution(df, "four_digit")
    assert out["mhc_restriction"].tolist() == ["HLA-A*02:01"]


def test_apply_min_resolution_keeps_donor_sets_above_two_digit():
    """A promoted donor set is more specific than a two-digit allele.

    hitlist stores ``donor_set`` for these rows; reclassifying the joined
    restriction string is what used to decide this, and disagreed with the
    stored annotation on pre-1.55.7 indexes.
    """
    df = _hits(["HLA-A*02:01", "HLA-A*01:01;HLA-A*23:01", "HLA-A*02", "HLA class I"])
    out = _apply_min_resolution(df, "donor_set")
    assert out["mhc_restriction"].tolist() == ["HLA-A*02:01", "HLA-A*01:01;HLA-A*23:01"]


def test_apply_min_resolution_reads_the_stored_annotation():
    """The stored label decides, not a re-derivation from the restriction."""
    df = _hits(["HLA-A*02:01"])
    df.loc[0, "allele_resolution"] = "class_only"
    assert _apply_min_resolution(df, "four_digit").empty


def test_apply_min_resolution_passthrough_when_none():
    df = _hits(["HLA-A*02:01", "HLA-A2"])
    assert _apply_min_resolution(df, None).equals(df)


def test_apply_min_resolution_requires_annotated_index():
    legacy = pd.DataFrame({"peptide": ["P"], "mhc_restriction": ["HLA-A*02:01"]})
    with pytest.raises(ValueError, match="allele_resolution"):
        _apply_min_resolution(legacy, "four_digit")
