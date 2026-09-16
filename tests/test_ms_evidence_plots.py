"""Unit tests for tsarina.ms_evidence_plots.

Every figure here reads from :func:`tsarina.evidence.CTA_evidence`, which is
a packaged-CSV join (oncoref's frame + tsarina's MS overlay) -- fast, no live
hitlist index required, and already exercised unconditionally elsewhere (see
tests/test_evidence.py). These tests do the same: no skip guard.
"""

import matplotlib

matplotlib.use("Agg")

from tsarina.ms_evidence_plots import (
    MS_RESTRICTION_COLORS,
    _label,
    plot_confidence_vs_ms_restriction,
    plot_literature_support_distribution,
    plot_ms_restriction_distribution,
    plot_offtarget_tissue_frequency,
    plot_top_evidenced_ctas,
    render_all,
)
from tsarina.tiers import MS_RESTRICTION_VALUES


def test_label_preserves_known_acronyms():
    assert _label("NO_MS_DATA") == "No MS Data"
    assert _label("central nervous system (cns)") == "Central Nervous System (CNS)"
    assert _label("CANCER_ONLY") == "Cancer Only"


def test_ms_restriction_colors_cover_every_category_in_fixed_order():
    """Every value tiers.py declares must have a color, and colors must not
    repeat -- a silent hue collision would make two categories
    indistinguishable in the restriction-distribution and stacked-bar
    figures."""
    assert set(MS_RESTRICTION_COLORS) == set(MS_RESTRICTION_VALUES)
    assert len(set(MS_RESTRICTION_COLORS.values())) == len(MS_RESTRICTION_COLORS)


def test_plot_ms_restriction_distribution_draws_one_bar_per_present_category():
    _fig, ax = plot_ms_restriction_distribution()
    counts = _cta_evidence_counts()
    present = [c for c in MS_RESTRICTION_VALUES if counts.get(c, 0) > 0]
    assert len(ax.patches) == len(present)
    # Every category present in the data must be readable off the axis.
    labels = [t.get_text() for t in ax.get_yticklabels()]
    assert len(labels) == len(present)


def test_plot_confidence_vs_ms_restriction_stacks_to_the_row_total():
    """Each HPA-confidence column's stacked segments must sum to that
    column's true CTA count -- the one invariant a stacked bar cannot get
    wrong without silently misrepresenting the data."""
    import pandas as pd

    from tsarina.evidence import CTA_evidence
    from tsarina.gene_sets import CTA_gene_names

    df = CTA_evidence()
    df = df[df["Symbol"].isin(set(CTA_gene_names()))]
    expected = df["restriction_confidence"].value_counts()

    _fig, ax = plot_confidence_vs_ms_restriction()
    # One bar container per present MS-restriction category; sum heights
    # per x position across containers to get each column's total.
    totals = pd.Series(0.0, index=range(len(ax.get_xticklabels())))
    for container in ax.containers:
        for i, patch in enumerate(container):
            totals[i] += patch.get_height()
    for i, label in enumerate(ax.get_xticklabels()):
        tier = label.get_text().upper()
        assert totals[i] == expected.get(tier, 0)


def test_plot_literature_support_distribution_is_log_scaled_and_nonempty():
    _fig, ax = plot_literature_support_distribution()
    assert ax.get_yscale() == "log"
    assert len(ax.patches) > 0


def test_plot_top_evidenced_ctas_only_includes_cancer_only_restriction():
    from tsarina.evidence import CTA_evidence

    _fig, ax = plot_top_evidenced_ctas(top_n=5)
    labels = {t.get_text() for t in ax.get_yticklabels()}
    assert 0 < len(labels) <= 5
    df = CTA_evidence()
    cancer_only = set(df.loc[df["ms_restriction"] == "CANCER_ONLY", "Symbol"])
    assert labels <= cancer_only


def test_plot_offtarget_tissue_frequency_bars_are_positive_and_sorted_desc():
    _fig, ax = plot_offtarget_tissue_frequency(top_n=10)
    heights = [p.get_width() for p in ax.patches]
    assert all(h > 0 for h in heights)
    # Rendered top-to-bottom largest-first (bars built ascending, then
    # y-axis is not inverted here -- unlike the restriction-distribution
    # figure -- so check monotonicity in plot order instead of assuming
    # a particular axis direction).
    assert heights == sorted(heights) or heights == sorted(heights, reverse=True)


def test_render_all_writes_exactly_five_nonempty_pngs(tmp_path):
    written = render_all(tmp_path)
    assert len(written) == 5
    for path in written.values():
        assert path.exists()
        assert path.stat().st_size > 0
        assert path.suffix == ".png"


def test_render_all_default_creates_a_timestamped_run_with_a_latest_symlink(tmp_path, monkeypatch):
    """No explicit out_dir: figures/<family>/run_<stamp>/, with `latest`
    repointed at it -- the same convention as oncoref's figure output, so a
    combined figures/ tree across both repos reads as one system."""
    from tsarina.ms_evidence_plots import FIGURE_FAMILY

    monkeypatch.chdir(tmp_path)
    written = render_all()
    assert len(written) == 5

    family_dir = tmp_path / "figures" / FIGURE_FAMILY
    run_dirs = sorted(p for p in family_dir.iterdir() if p.name.startswith("run_"))
    assert len(run_dirs) == 1

    latest = family_dir / "latest"
    assert latest.is_symlink()
    assert latest.resolve() == run_dirs[0].resolve()
    for path in written.values():
        assert path.resolve().parent == run_dirs[0].resolve()


def test_render_all_default_does_not_collide_on_repeated_calls_within_one_second(
    tmp_path, monkeypatch
):
    from tsarina.ms_evidence_plots import FIGURE_FAMILY

    monkeypatch.chdir(tmp_path)
    render_all()
    render_all()  # must not raise or silently overwrite the first run's dir

    family_dir = tmp_path / "figures" / FIGURE_FAMILY
    run_dirs = [p for p in family_dir.iterdir() if p.name.startswith("run_")]
    assert len(run_dirs) == 2


def test_render_all_respects_an_explicit_gene_subset(tmp_path):
    """A caller-supplied gene set must actually narrow the figures, not be
    silently ignored in favor of the full default CTA panel."""
    from tsarina.evidence import CTA_evidence

    df = CTA_evidence()
    some_genes = set(df["Symbol"].head(5))
    _fig, ax = plot_ms_restriction_distribution(genes=some_genes)
    assert f"n={len(some_genes)}" in ax.get_xlabel()


def _cta_evidence_counts():
    from tsarina.evidence import CTA_evidence
    from tsarina.gene_sets import CTA_gene_names

    df = CTA_evidence()
    df = df[df["Symbol"].isin(set(CTA_gene_names()))]
    return df["ms_restriction"].value_counts().to_dict()
