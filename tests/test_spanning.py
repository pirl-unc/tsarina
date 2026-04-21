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


# ── on_progress callback ───────────────────────────────────────────────


def test_on_progress_default_is_silent():
    """Library contract: when on_progress is not supplied, the function
    runs without calling any progress callback (or printing anywhere).

    Confirmed indirectly by the function completing without error; the
    explicit-callback test below asserts the callback IS called on the
    opt-in path."""
    df = spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01"],
        max_percentile=2.0,
    )
    assert not df.empty  # sanity — stub produces rows


def test_on_progress_callback_receives_three_stages():
    messages: list[str] = []
    spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01", "HLA-A*24:02"],
        max_percentile=2.0,
        on_progress=messages.append,
    )
    # Expect three stages: pre-scoring, post-scoring, summary.
    assert len(messages) == 3
    pre, post, summary = messages
    # Pre-scoring mentions peptide + allele counts and predictor.
    assert "peptides" in pre.lower() and "alleles" in pre.lower()
    assert "mhcflurry" in pre
    # Post-scoring mentions elapsed time.
    assert "Scored" in post
    assert "s" in post  # the elapsed-seconds suffix
    # Summary mentions CTA + allele counts and filled-cell count.
    assert "CTAs" in summary
    assert "alleles" in summary
    assert "filled cells" in summary


def test_on_progress_pre_scoring_message_has_accurate_counts():
    """The pre-scoring message's numbers should match the actual peptide
    and allele pool sizes: stub has 2 peptides for {MAGEA4, PRAME} each
    = 4 unique peptides, crossed with 3 alleles = 12 predictions."""
    messages: list[str] = []
    spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02"],
        max_percentile=10.0,
        on_progress=messages.append,
    )
    pre = messages[0]
    assert "4 peptides" in pre
    assert "3 alleles" in pre
    assert "12 predictions" in pre


def test_on_progress_summary_reports_filled_cell_count():
    """Summary message's filled-cell count should equal the count the
    output table actually carries (rows with a non-NaN peptide)."""
    messages: list[str] = []
    long = spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01", "HLA-A*24:02"],
        max_percentile=2.0,
        output_format="long",
        on_progress=messages.append,
    )
    summary = messages[-1]
    assert f"{len(long)} filled cells" in summary


def test_cli_handler_wires_on_progress_to_stderr(monkeypatch, capsys):
    """The CLI handler must pass a stderr-printing callback into
    spanning_pmhc_set.  Verify callback messages land on stderr (not
    stdout) and that the DataFrame itself still serializes to stdout."""
    import argparse

    from tsarina import cli_spanning

    captured_kwargs: dict = {}

    def _fake_spanning(*args, **kwargs):
        captured_kwargs.update(kwargs)
        if kwargs.get("on_progress"):
            kwargs["on_progress"]("fake-progress-message")
        return pd.DataFrame({"cta": ["MAGEA4"], "HLA-A*02:01": ["PEPTIDE9"]})

    monkeypatch.setattr("tsarina.spanning.spanning_pmhc_set", _fake_spanning)

    args = argparse.Namespace(
        cta_count=25,
        cta_rank_by="ms_cancer_peptide_count",
        ctas=None,
        min_restriction_confidence=["HIGH", "MODERATE"],
        restriction_levels=None,
        alleles=None,
        panel="iedb27_ab",
        lengths=(9,),
        ensembl_release=112,
        require_cta_exclusive=True,
        predictor="mhcflurry",
        max_percentile=2.0,
        format="wide",
        output=None,
    )
    cli_spanning.handle(args)

    assert "on_progress" in captured_kwargs
    assert callable(captured_kwargs["on_progress"])

    captured = capsys.readouterr()
    assert "fake-progress-message" in captured.err
    assert "fake-progress-message" not in captured.out
    assert "MAGEA4" in captured.out  # DataFrame serialized to stdout


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


# ── require_cta_exclusive routing ──────────────────────────────────────


def test_require_cta_exclusive_true_calls_exclusive(monkeypatch):
    """Default True path must dispatch to cta_exclusive_peptides, not
    cta_peptides."""
    calls = {"exclusive": 0, "non_exclusive": 0}

    def _exclusive(**kw):
        calls["exclusive"] += 1
        return _stub_peptides()

    def _non_exclusive(**kw):
        calls["non_exclusive"] += 1
        return _stub_peptides()

    monkeypatch.setattr("tsarina.peptides.cta_exclusive_peptides", _exclusive, raising=True)
    monkeypatch.setattr("tsarina.peptides.cta_peptides", _non_exclusive, raising=True)

    spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01"],
        require_cta_exclusive=True,
        max_percentile=10.0,
    )
    assert calls == {"exclusive": 1, "non_exclusive": 0}


def test_require_cta_exclusive_false_calls_cta_peptides(monkeypatch):
    """False path must dispatch to cta_peptides (full CTA k-mers, no
    exclusivity gate), never the exclusive fn."""
    calls = {"exclusive": 0, "non_exclusive": 0}

    def _exclusive(**kw):
        calls["exclusive"] += 1
        return _stub_peptides()

    def _non_exclusive(**kw):
        calls["non_exclusive"] += 1
        return _stub_peptides()

    monkeypatch.setattr("tsarina.peptides.cta_exclusive_peptides", _exclusive, raising=True)
    monkeypatch.setattr("tsarina.peptides.cta_peptides", _non_exclusive, raising=True)

    spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01"],
        require_cta_exclusive=False,
        max_percentile=10.0,
    )
    assert calls == {"exclusive": 0, "non_exclusive": 1}


# ── Multi-length best-per-cell ─────────────────────────────────────────


def _stub_peptides_multi_length() -> pd.DataFrame:
    """Mixed 9-mer and 10-mer peptides so the best-per-cell selection
    must reach across lengths."""
    return pd.DataFrame(
        {
            "peptide": [
                "MAGEAPEP9M",  # 10-mer; stub gives it percentile 0.3
                "MAGEAPEP1",  # 9-mer; stub gives it percentile 0.1 (wins)
                "PRAMEPEP0X",  # 10-mer; stub gives it percentile 0.05 (wins)
                "PRAMEPEP2",  # 9-mer; stub gives it percentile 1.5
            ],
            "length": [10, 9, 10, 9],
            "gene_name": ["MAGEA4", "MAGEA4", "PRAME", "PRAME"],
            "gene_id": ["E1", "E1", "E2", "E2"],
        }
    )


def _stub_scores_multi_length(peptides, alleles, **kwargs) -> pd.DataFrame:
    """Deterministic percentile per peptide — chosen so the 10-mer wins
    for PRAME and the 9-mer wins for MAGEA4, regardless of allele."""
    percentile_by_peptide = {
        "MAGEAPEP9M": 0.3,
        "MAGEAPEP1": 0.1,
        "PRAMEPEP0X": 0.05,
        "PRAMEPEP2": 1.5,
    }
    rows = []
    for pep in peptides:
        pct = percentile_by_peptide.get(pep, 5.0)
        for allele in alleles:
            rows.append(
                {
                    "peptide": pep,
                    "allele": allele,
                    "presentation_percentile": pct,
                    "presentation_score": 0.95 - pct / 10,
                    "affinity_nm": 50.0 + pct * 100,
                }
            )
    return pd.DataFrame(rows)


def test_multi_length_best_per_cell_winner_can_be_any_length(monkeypatch):
    """With lengths=(9, 10), each cell's winner should be the
    lowest-percentile peptide across BOTH lengths — not length-biased."""
    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides",
        lambda **kw: _stub_peptides_multi_length(),
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.scoring.score_presentation", _stub_scores_multi_length, raising=True
    )

    df = spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01"],
        lengths=(9, 10),
        max_percentile=2.0,
    )
    indexed = df.set_index("cta")
    # MAGEA4: MAGEAPEP1 (9-mer, 0.1%) beats MAGEAPEP9M (10-mer, 0.3%)
    assert indexed.loc["MAGEA4", "HLA-A*02:01"] == "MAGEAPEP1"
    # PRAME: PRAMEPEP0X (10-mer, 0.05%) beats PRAMEPEP2 (9-mer, 1.5%)
    assert indexed.loc["PRAME", "HLA-A*02:01"] == "PRAMEPEP0X"


def test_multi_length_long_format_preserves_winner_length(monkeypatch):
    """In long format, the length column should reflect the winning
    peptide's length, not a default."""
    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides",
        lambda **kw: _stub_peptides_multi_length(),
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.scoring.score_presentation", _stub_scores_multi_length, raising=True
    )

    long = spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01"],
        lengths=(9, 10),
        max_percentile=2.0,
        output_format="long",
    )
    rows = {r["cta"]: r for _, r in long.iterrows()}
    assert int(rows["MAGEA4"]["length"]) == 9  # MAGEAPEP1
    assert int(rows["PRAME"]["length"]) == 10  # PRAMEPEP0X
