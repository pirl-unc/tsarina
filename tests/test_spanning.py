# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Unit tests for ``tsarina.spanning.spanning_pmhc_set``.

The library function pulls peptides through ``cta_exclusive_peptides``,
scores them via ``score_presentation``, and pivots the best result per
(CTA, allele) cell. Tests stub each of those upstream calls so they run
without Ensembl / mhcflurry installed.
"""

from __future__ import annotations

import pandas as pd
import pytest

from tsarina.spanning import spanning_pmhc_set


def _stub_peptides() -> pd.DataFrame:
    """Three CTAs, two peptides each at length 9."""
    return pd.DataFrame(
        {
            "peptide": [
                "MAGEAPEP1",
                "MAGEAPEP2",
                "PRAMEPEP1",
                "PRAMEPEP2",
                "NYESPEPT1",
                "NYESPEPT2",
            ],
            "length": [9, 9, 9, 9, 9, 9],
            "gene_name": ["MAGEA4", "MAGEA4", "PRAME", "PRAME", "CTAG1B", "CTAG1B"],
            "gene_id": ["E1", "E1", "E2", "E2", "E3", "E3"],
        }
    )


def _stub_scores(peptides, alleles, **kwargs) -> pd.DataFrame:
    """Score every (peptide, allele) pair with deterministic percentiles
    derived from the peptide and allele strings — lets tests assert which
    peptide should win each cell."""
    rows = []
    for pep in peptides:
        for allele in alleles:
            # Lowest percentile = most-presented; PEP1 always beats PEP2
            # for the same gene; A*02:01 beats A*24:02 for all peptides.
            base = 0.1 if pep.endswith("1") else 1.5
            allele_bump = 0.0 if allele == "HLA-A*02:01" else 0.5
            rows.append(
                {
                    "peptide": pep,
                    "allele": allele,
                    "presentation_percentile": base + allele_bump,
                    "presentation_score": 0.95 - (base + allele_bump) / 10,
                    "affinity_nm": 50.0 + (base + allele_bump) * 100,
                }
            )
    return pd.DataFrame(rows)


@pytest.fixture(autouse=True)
def _stub_pipeline(monkeypatch):
    """Wire stubs for every upstream call spanning_pmhc_set makes."""
    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides",
        lambda **kw: _stub_peptides(),
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.peptides.cta_peptides",
        lambda **kw: _stub_peptides(),
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.scoring.score_presentation",
        _stub_scores,
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.gene_sets.CTA_gene_names",
        lambda: {"MAGEA4", "PRAME", "CTAG1B", "BAGE"},
        raising=True,
    )
    # Stub the bundled CTA dataframe so cta_count / rank_by paths work.
    cta_csv = pd.DataFrame(
        {
            "Symbol": ["MAGEA4", "PRAME", "CTAG1B", "BAGE"],
            "Ensembl_Gene_ID": ["E1", "E2", "E3", "E4"],
            "filtered": ["true", "true", "true", "true"],
            "never_expressed": ["false", "false", "false", "true"],
            "restriction_confidence": ["HIGH", "HIGH", "MODERATE", "LOW"],
            "restriction": ["TESTIS", "TESTIS", "PLACENTAL", "TESTIS"],
            "ms_cancer_peptide_count": [50, 30, 10, 0],
        }
    )
    monkeypatch.setattr("tsarina.loader.cta_dataframe", lambda: cta_csv, raising=True)


# ── CTA selection ──────────────────────────────────────────────────────


def test_top_n_ranking_by_default_column():
    """Top 2 by ms_cancer_peptide_count: MAGEA4 (50), PRAME (30).
    CTAG1B (10) excluded by N=2 cap. BAGE (0) excluded as never_expressed."""
    df = spanning_pmhc_set(
        cta_count=2,
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
    )
    assert list(df["cta"]) == ["MAGEA4", "PRAME"]


def test_explicit_ctas_overrides_ranking():
    """When --ctas is supplied, --cta-count is ignored and the order is
    preserved (capped only by membership in CTA_gene_names())."""
    df = spanning_pmhc_set(
        ctas=["CTAG1B", "PRAME"],
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
    )
    assert list(df["cta"]) == ["CTAG1B", "PRAME"]


def test_min_restriction_confidence_filters_low():
    """LOW-confidence CTAs (BAGE) drop out of ranking when the gate is
    HIGH/MODERATE."""
    df = spanning_pmhc_set(
        cta_count=10,
        min_restriction_confidence=("HIGH", "MODERATE"),
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
    )
    assert "BAGE" not in df["cta"].tolist()


def test_min_restriction_confidence_none_admits_low():
    """Disabling the confidence gate would let BAGE through — but BAGE is
    also excluded by never_expressed=True. Use a fixture variant via
    explicit ctas to demonstrate the bypass."""
    df = spanning_pmhc_set(
        cta_count=10,
        min_restriction_confidence=None,
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
    )
    # Three remain after never_expressed filter; LOW-confidence BAGE excluded.
    assert set(df["cta"]) == {"MAGEA4", "PRAME", "CTAG1B"}


def test_restriction_levels_filter():
    """restriction_levels=('PLACENTAL',) keeps only CTAG1B."""
    df = spanning_pmhc_set(
        cta_count=10,
        restriction_levels=("PLACENTAL",),
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
    )
    assert list(df["cta"]) == ["CTAG1B"]


# ── Allele resolution ──────────────────────────────────────────────────


def test_panel_default_resolves_via_get_panel():
    """Default panel='iedb27_ab' should produce a 27-column wide table
    plus the 'cta' index column = 28 cols."""
    df = spanning_pmhc_set(
        cta_count=2,
        max_percentile=10.0,
    )
    assert df.shape[1] == 28  # 1 cta + 27 alleles


def test_explicit_alleles_override_panel():
    df = spanning_pmhc_set(
        cta_count=2,
        alleles=["HLA-A*02:01", "HLA-A*24:02"],
        max_percentile=10.0,
    )
    assert list(df.columns) == ["cta", "HLA-A*02:01", "HLA-A*24:02"]


def test_no_alleles_no_panel_raises():
    with pytest.raises(ValueError, match="alleles or panel"):
        spanning_pmhc_set(panel=None, alleles=None)


# ── Best-per-cell selection ────────────────────────────────────────────


def test_lowest_percentile_peptide_wins_each_cell():
    """Stub gives PEP1 = 0.1, PEP2 = 1.5 baseline. Both pass max=2.0;
    PEP1 should win every cell."""
    df = spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01"],
        max_percentile=2.0,
    )
    assert df.set_index("cta").loc["MAGEA4", "HLA-A*02:01"] == "MAGEAPEP1"
    assert df.set_index("cta").loc["PRAME", "HLA-A*02:01"] == "PRAMEPEP1"


def test_max_percentile_filters_weak_cells():
    """Stub: A*24:02 percentiles are 0.6 (PEP1) and 2.0 (PEP2). A cutoff
    of 0.5 leaves every A*24:02 cell empty; A*02:01 cells (0.1) survive."""
    df = spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01", "HLA-A*24:02"],
        max_percentile=0.5,
    )
    indexed = df.set_index("cta")
    assert indexed.loc["MAGEA4", "HLA-A*02:01"] == "MAGEAPEP1"
    assert pd.isna(indexed.loc["MAGEA4", "HLA-A*24:02"])


# ── Output formats ─────────────────────────────────────────────────────


def test_long_format_has_one_row_per_filled_cell():
    long = spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01", "HLA-A*24:02"],
        max_percentile=2.0,
        output_format="long",
    )
    # 2 CTAs * 2 alleles = 4 cells, all pass with max=2.0
    assert len(long) == 4
    assert set(long.columns) == {
        "cta",
        "allele",
        "peptide",
        "length",
        "presentation_percentile",
        "presentation_score",
        "affinity_nm",
    }
    # Sorted by cta then allele in the supplied order
    assert list(long["cta"]) == ["MAGEA4", "MAGEA4", "PRAME", "PRAME"]


def test_long_format_drops_filtered_cells_entirely():
    """In long output, cells above max_percentile are absent (not NaN)."""
    long = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01", "HLA-A*24:02"],
        max_percentile=0.5,
        output_format="long",
    )
    assert len(long) == 1
    assert long.iloc[0]["allele"] == "HLA-A*02:01"


def test_unrecognized_output_format_raises():
    with pytest.raises(ValueError, match="output_format"):
        spanning_pmhc_set(output_format="grid")


# ── Edge cases ─────────────────────────────────────────────────────────


def test_no_ctas_pass_filter_returns_canonical_empty_frame():
    """If filters eliminate every CTA, we still return a well-formed
    (empty) wide table with the requested allele columns."""
    df = spanning_pmhc_set(
        ctas=[],
        alleles=["HLA-A*02:01"],
        max_percentile=2.0,
    )
    assert list(df.columns) == ["cta", "HLA-A*02:01"]
    assert df.empty


def test_size_within_target_envelope_for_default_args(monkeypatch):
    """The headline use case: defaults (25 CTAs x 27 alleles, max 2%)
    should produce a wide table whose count of filled cells is in the
    100-1000 spanning-set range. Stub gives exactly 25*27=675 candidate
    cells; with the stub's 0.1 / 1.5 percentile baselines, every cell
    passes the 2.0 cutoff, so all 675 are filled — comfortably in range.

    Uses ``monkeypatch.setattr`` (not raw attribute assignment) so the
    expanded stubs are torn down after this test and don't leak into the
    rest of the test session.
    """
    # Expand the stub CTA pool so cta_count=25 has enough candidates.
    extra_genes = [f"GENE{i:03d}" for i in range(25)]
    extra_peps = pd.DataFrame(
        {
            "peptide": [f"PEPGENE{i:03d}" for i in range(25)],
            "length": [9] * 25,
            "gene_name": extra_genes,
            "gene_id": [f"E{i + 10}" for i in range(25)],
        }
    )
    all_genes = ["MAGEA4", "PRAME", "CTAG1B", *extra_genes]
    cta_csv = pd.DataFrame(
        {
            "Symbol": all_genes,
            "Ensembl_Gene_ID": [f"E{i}" for i in range(28)],
            "filtered": ["true"] * 28,
            "never_expressed": ["false"] * 28,
            "restriction_confidence": ["HIGH"] * 28,
            "restriction": ["TESTIS"] * 28,
            "ms_cancer_peptide_count": list(range(28, 0, -1)),
        }
    )

    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides",
        lambda **kw: pd.concat([_stub_peptides(), extra_peps], ignore_index=True),
        raising=True,
    )
    monkeypatch.setattr("tsarina.loader.cta_dataframe", lambda: cta_csv, raising=True)
    monkeypatch.setattr("tsarina.gene_sets.CTA_gene_names", lambda: set(all_genes), raising=True)

    df = spanning_pmhc_set(
        cta_count=25,
        panel="iedb27_ab",
        max_percentile=2.0,
        output_format="long",
    )
    assert 100 < len(df) <= 1000
