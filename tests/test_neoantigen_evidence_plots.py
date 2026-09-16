"""Unit tests for tsarina.neoantigen_evidence_plots.

Split in two tiers, matching tests/test_hitlist_integration.py's convention:

- Pure-Python logic (parsing, window/epitope overlap math, plot rendering)
  is tested against small synthetic data, always runs, no hitlist needed.
- The one thing that genuinely needs real data -- does hitlist's public MS
  index actually contain a hit for this fixture's peptides -- is gated on
  ``hitlist.observations.is_built()`` and skips cleanly without it.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pandas as pd
import pytest

from tsarina.neoantigen_evidence_plots import (
    _gene_hit_layout,
    _minimal_epitope_for_gene,
    load_vaccine_peptide_table,
    plot_ms_hit_ranking,
    plot_sequence_overlay,
    plot_tissue_provenance,
)

_FIXTURE = Path(__file__).parent / "data" / "osteosarc" / "vaccine_overlap_summary.tsv"


# ── load_vaccine_peptide_table ──────────────────────────────────────────


def test_load_vaccine_peptide_table_parses_the_real_fixture():
    df = load_vaccine_peptide_table(_FIXTURE)
    assert not df.empty
    assert set(df.columns) == {"gene", "mutation", "construct", "peptide"}
    assert "SMC5" in set(df["gene"])
    row = df[(df["gene"] == "SMC5") & (df["construct"] == "mRNA_full")]
    assert row["peptide"].iloc[0] == "EELQQALIVKQNEELDRQKRIGNTRKMIEDL"


def test_load_vaccine_peptide_table_splits_comma_separated_elispot_variants():
    """MAP2's elispot cell holds two peptides separated by a comma -- both
    must become their own row, not one malformed concatenated string."""
    df = load_vaccine_peptide_table(_FIXTURE)
    sub = df[(df["gene"] == "MAP2") & (df["construct"] == "JLF_elispot")]
    assert set(sub["peptide"]) == {
        "KTVRIYQGRVVPFTKALMIKFEE",
        "KTVRIYQGRVVPFTKALMIKFEEI",
    }


def test_load_vaccine_peptide_table_rejects_placeholder_values():
    """'NA', '?', and empty cells must never become peptide rows."""
    df = load_vaccine_peptide_table(_FIXTURE)
    assert not (df["peptide"] == "NA").any()
    assert not (df["peptide"] == "?").any()
    assert not (df["peptide"] == "").any()
    # KRT18's CeGaT cell is literally "?" in the source -- confirm it produced no row.
    assert "KRT18" not in set(df["gene"])


def test_load_vaccine_peptide_table_dedupes_identical_sequences_across_constructs():
    """BTD's mRNA_minimal/JLF_elispot/JLF_V3/JLF_V2/JLF_V1 constructs share
    the identical 23-residue peptide -- one row, not five."""
    df = load_vaccine_peptide_table(_FIXTURE)
    sub = df[(df["gene"] == "BTD") & (df["peptide"] == "RPTLSKELYALGVVDGLHTVHGT")]
    assert len(sub) == 1


# ── _minimal_epitope_for_gene ────────────────────────────────────────────


def test_minimal_epitope_for_gene_reads_the_mrna_minimal_column():
    df = load_vaccine_peptide_table(_FIXTURE)
    epitope, is_minimal = _minimal_epitope_for_gene(df, "SMC5")
    assert epitope == "RQKRIGNTR"
    assert is_minimal is True


def test_minimal_epitope_for_gene_prefers_the_shorter_elispot_variant():
    """EPG5 has no true mRNA_minimal construct, only two comma-separated
    ELISPOT testing peptides of different lengths; the shorter one is used,
    but it is flagged as not a true minimal epitope."""
    df = load_vaccine_peptide_table(_FIXTURE)
    epitope, is_minimal = _minimal_epitope_for_gene(df, "EPG5")
    assert epitope == "KELPLYLWQPSTSEIAVIRD"
    assert is_minimal is False


def test_minimal_epitope_for_gene_returns_none_when_undocumented():
    df = load_vaccine_peptide_table(_FIXTURE)
    assert _minimal_epitope_for_gene(df, "CD109") is None


# ── epitope/hit overlap math (the safety-relevant core logic) ───────────


def _synthetic_pipeline(gene_peptide, epitope, hit_windows):
    """Build a minimal (df, window_source, window_info) triple by hand,
    bypassing hitlist entirely, to test the overlap arithmetic in isolation."""
    gene, peptide = gene_peptide
    df = pd.DataFrame(
        [
            {"gene": gene, "mutation": "p.X1Y", "construct": "mRNA_full", "peptide": peptide},
            {"gene": gene, "mutation": "p.X1Y", "construct": "mRNA_minimal", "peptide": epitope},
        ]
    )
    window_source = {}
    window_info = {}
    for window in hit_windows:
        pos = peptide.find(window)
        assert pos >= 0, f"test setup error: {window!r} not in {peptide!r}"
        window_source[window] = [(gene, peptide, pos, len(window))]
        window_info[window] = {
            "n": 1,
            "categories": {"cancer": 1},
            "tissues": ["Blood"],
            "pmids": ["1"],
        }
    return df, window_source, window_info


def test_gene_hit_layout_flags_a_hit_window_that_overlaps_the_epitope():
    # peptide:  A B C D E F G H I J K
    # epitope:      C D E F G           (positions 2-7)
    # hit:              E F G H I       (positions 4-9) -- overlaps epitope at E,F,G
    peptide = "ABCDEFGHIJK"
    df, window_source, window_info = _synthetic_pipeline(
        ("GENE1", peptide), epitope="CDEFG", hit_windows=["EFGHI"]
    )
    layout = _gene_hit_layout(df, window_source, window_info)
    assert len(layout) == 1
    row = layout[0]
    assert row["epitope_span"] == (2, 7)
    (window, start, end, _category, overlaps) = row["hit_spans"][0]
    assert (window, start, end) == ("EFGHI", 4, 9)
    assert overlaps is True


def test_gene_hit_layout_does_not_flag_a_purely_flanking_hit():
    # epitope at 2-7 ("CDEFG"); hit at 7-11 ("HIJK") starts exactly where
    # the epitope ends -- adjacent, not overlapping.
    peptide = "ABCDEFGHIJK"
    df, window_source, window_info = _synthetic_pipeline(
        ("GENE2", peptide), epitope="CDEFG", hit_windows=["HIJK"]
    )
    layout = _gene_hit_layout(df, window_source, window_info)
    (_window, _start, _end, _category, overlaps) = layout[0]["hit_spans"][0]
    assert overlaps is False


def test_gene_hit_layout_leaves_overlap_undetermined_without_a_known_epitope():
    peptide = "ABCDEFGHIJK"
    df = pd.DataFrame(
        [{"gene": "GENE3", "mutation": "p.X1Y", "construct": "mRNA_full", "peptide": peptide}]
    )
    window_source = {"DEF": [("GENE3", peptide, 3, 3)]}
    window_info = {
        "DEF": {"n": 1, "categories": {"cancer": 1}, "tissues": ["Blood"], "pmids": ["1"]}
    }
    layout = _gene_hit_layout(df, window_source, window_info)
    assert layout[0]["epitope_span"] is None
    assert layout[0]["hit_spans"][0][-1] is False  # no epitope known -> cannot claim overlap


def test_gene_hit_layout_flags_true_minimal_epitope_as_minimal():
    peptide = "ABCDEFGHIJK"
    df, window_source, window_info = _synthetic_pipeline(
        ("GENE1", peptide), epitope="CDEFG", hit_windows=["EFGHI"]
    )
    layout = _gene_hit_layout(df, window_source, window_info)
    assert layout[0]["epitope_is_minimal"] is True


def test_gene_hit_layout_flags_elispot_fallback_as_not_minimal():
    """Regression pin for the EPG5 bug: a gene with no mRNA_minimal
    construct falls back to its (much longer) JLF_elispot testing peptide,
    and the layout must flag that region as not a true minimal epitope so
    the plot doesn't imply the whole span is the localized mutation."""
    peptide = "ABCDEFGHIJKLMNOPQRSTUV"
    gene = "GENE4"
    df = pd.DataFrame(
        [
            {"gene": gene, "mutation": "p.X1Y", "construct": "mRNA_full", "peptide": peptide},
            {
                "gene": gene,
                "mutation": "p.X1Y",
                "construct": "JLF_elispot",
                "peptide": "CDEFGHIJKLMNOPQRST",
            },
        ]
    )
    window_source = {"EFGHI": [(gene, peptide, 4, 5)]}
    window_info = {
        "EFGHI": {"n": 1, "categories": {"cancer": 1}, "tissues": ["Blood"], "pmids": ["1"]}
    }
    layout = _gene_hit_layout(df, window_source, window_info)
    assert layout[0]["epitope_is_minimal"] is False
    assert layout[0]["epitope_span"] == (2, 20)


# ── plotting: renders without error on synthetic data ───────────────────


def test_plot_sequence_overlay_renders():
    peptide = "ABCDEFGHIJK"
    df, window_source, window_info = _synthetic_pipeline(
        ("GENE1", peptide), epitope="CDEFG", hit_windows=["EFGHI"]
    )
    layout = _gene_hit_layout(df, window_source, window_info)
    fig, ax = plot_sequence_overlay(layout, window_info)
    assert fig is not None
    assert len(ax.texts) > 0


def test_plot_ms_hit_ranking_includes_zero_evidence_genes():
    peptide = "ABCDEFGHIJK"
    df, window_source, window_info = _synthetic_pipeline(
        ("GENE1", peptide), epitope="CDEFG", hit_windows=["EFGHI"]
    )
    layout = _gene_hit_layout(df, window_source, window_info)
    _fig, ax = plot_ms_hit_ranking(["GENE1", "GENE_WITH_NOTHING"], layout, window_info)
    labels = [t.get_text() for t in ax.get_yticklabels()]
    assert "GENE1" in labels
    assert "GENE_WITH_NOTHING" in labels


def test_plot_tissue_provenance_renders():
    hits = pd.DataFrame(
        [
            {
                "peptide": "EFGHI",
                "source_tissue": "Blood",
                "disease": "",
                "src_cancer": True,
                "src_healthy_tissue": False,
                "pmid": "1",
            },
            {
                "peptide": "EFGHI",
                "source_tissue": "Skin",
                "disease": "healthy",
                "src_cancer": False,
                "src_healthy_tissue": True,
                "pmid": "2",
            },
        ]
    )
    fig, axes = plot_tissue_provenance(hits, ["EFGHI"])
    assert fig is not None
    assert len(axes) == 1


# ── real-data integration: gated on a built local hitlist index ─────────


def _hitlist_index_built() -> bool:
    try:
        from hitlist.observations import is_built
    except ImportError:
        return False
    return is_built()


@pytest.mark.skipif(
    not _hitlist_index_built(), reason="hitlist observations index not built locally"
)
def test_render_all_finds_the_known_smc5_hit_against_real_hitlist_data(tmp_path):
    """Regression pin: SMC5's flanking window has real public MS support
    (verified manually against the built index -- 71 observations, 15
    tissues, including a bone osteosarcoma cell line, PMID 36010968).
    If this ever returns nothing, either hitlist's corpus changed
    materially or the windowing/lookup logic broke."""
    from tsarina.neoantigen_evidence_plots import render_all

    written = render_all(_FIXTURE, out_dir=tmp_path)
    assert set(written) == {"sequence-overlay", "tissue-provenance", "ms-hit-ranking"}
    for path in written.values():
        assert path.exists() and path.stat().st_size > 0
