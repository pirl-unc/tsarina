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

"""Publication figures for tsarina's mass-spec evidence layer.

oncoref decides CTA membership from HPA normal-tissue restriction alone --
no mass-spec evidence enters that call (see ``oncoref.cta`` and
``oncoref.cta_curation_plots``, which figure the *definition* funnel).  This
module figures the layer tsarina adds on top: does public immunopeptidomics
mass spec (IEDB/CEDAR, via hitlist) actually corroborate a CTA's restriction,
and what does the safety picture look like across the panel.

Two candidate tumor-antigen categories currently have a mass-spec evidence
surface worth figuring:

- **CTAs** -- :func:`tsarina.evidence.CTA_evidence`'s ``ms_restriction`` /
  ``ms_healthy_somatic_tissues`` / ``ms_pmids`` columns, a versioned overlay
  (``data/gene-ms-safety-evidence.csv``) computed once against IEDB/CEDAR and
  packaged with the release. Fast: reads the packaged CSV, no live hitlist
  index required.
- **Viral oncoantigens** -- :func:`tsarina.viral.viral_iedb_overlap`, run
  per virus. Slower: enumerates each virus's proteome and queries hitlist
  live, so figure generation here is opt-in (``include_viral=True``).

Mutation-derived neoantigens are not figured here: they are patient-private
by construction, so a public MS evidence corpus has essentially nothing to
say about any specific one -- there is no general-population category to
plot the way there is for CTAs or a fixed set of oncoviruses.

Every figure function takes an optional ``ax``; :func:`render_all` builds the
full set and saves one PNG per figure, to a fresh timestamped run directory
under ``figures/ms-evidence-cta/`` by default (with a ``latest`` symlink) --
or to an explicit path when one is given. Design intent is print-ready:
white background, no baked-in title (captions belong in the manuscript),
minimal in-figure text -- selective direct labels, no per-bar value labels
except where they are the point of the chart.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from .tiers import MS_RESTRICTION_VALUES

#: Colorblind-safe categorical palette (validated: adjacent-pair CVD Delta E
#: >= 8 in both light and dark surfaces). Assigned in fixed order to
#: :data:`tsarina.tiers.MS_RESTRICTION_VALUES`, never re-cycled per plot.
_CATEGORICAL_HUES = (
    "#2a78d6",  # blue
    "#eb6834",  # orange
    "#1baf7a",  # aqua
    "#eda100",  # yellow
    "#e87ba4",  # magenta
    "#4a3aa7",  # violet
)

#: MS restriction category -> fixed hue. Built once from the shared
#: vocabulary so a plot never silently drifts from ``tiers.py``.
MS_RESTRICTION_COLORS: dict[str, str] = dict(zip(MS_RESTRICTION_VALUES, _CATEGORICAL_HUES))

#: Muted ink for axis text -- never the series color (text carries no
#: identity here, the mark next to it does).
_TEXT_COLOR = "#3a3a38"
_GRID_COLOR = "#e4e3df"

#: Acronyms `str.title()` would otherwise mangle ("No Ms Data", "(Cns)").
_ACRONYMS = ("MS", "CNS")


def _label(text: str) -> str:
    """Title-case a snake/space-separated label, preserving acronyms."""
    words = text.replace("_", " ").split(" ")
    out = []
    for word in words:
        upper = word.upper().strip("()")
        if upper in _ACRONYMS:
            out.append(word.upper())
        else:
            out.append(word.title())
    return " ".join(out)


def _apply_publication_style() -> None:
    """Set matplotlib rcParams for a white-background, low-chrome figure.

    Deliberately minimal: opaque white surfaces (never transparent -- a
    transparent PNG placed on a dark manuscript background is unreadable),
    thin hairline spines with the top/right pair removed, no title/legend
    frame, and muted-ink text so a colored mark -- not the text -- carries
    category identity.
    """
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.transparent": False,
            "savefig.dpi": 300,
            "figure.dpi": 150,
            "font.size": 9,
            "font.family": "sans-serif",
            "text.color": _TEXT_COLOR,
            "axes.edgecolor": _TEXT_COLOR,
            "axes.labelcolor": _TEXT_COLOR,
            "xtick.color": _TEXT_COLOR,
            "ytick.color": _TEXT_COLOR,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "axes.grid": False,
            "grid.color": _GRID_COLOR,
            "grid.linewidth": 0.6,
            "legend.frameon": False,
            "pdf.fonttype": 42,  # editable text in vector output, not paths
            "ps.fonttype": 42,
        }
    )


def _new_ax(ax, figsize: tuple[float, float]):
    _apply_publication_style()
    if ax is not None:
        return ax.figure, ax
    import matplotlib.pyplot as plt

    return plt.subplots(figsize=figsize)


def _cta_evidence_frame(genes: set[str] | list[str] | None = None) -> pd.DataFrame:
    """oncoref's default CTA set, joined with tsarina's MS overlay.

    ``genes=None`` uses oncoref's default expressed CTA set
    (:func:`tsarina.gene_sets.CTA_gene_names`) -- the set an evidence figure
    should describe, not the pre-filter candidate universe (that funnel is
    oncoref's to figure).
    """
    from .evidence import CTA_evidence
    from .gene_sets import CTA_gene_names

    df = CTA_evidence()
    wanted = set(genes) if genes is not None else set(CTA_gene_names())
    return df[df["Symbol"].isin(wanted)].copy()


def _pmid_count(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).apply(lambda s: len([p for p in s.split(";") if p]))


def _somatic_tissues(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str)


# ── Figure 1: MS restriction category distribution ─────────────────────────


def plot_ms_restriction_distribution(genes=None, ax=None, figsize=(4.2, 3.2)):
    """Horizontal bar of CTAs by mass-spec restriction category.

    Answers: of the CTA panel, how many have direct public MS evidence at
    all, and of those, how many are cancer-specific versus seen on healthy
    tissue by mass spec.
    """
    fig, ax = _new_ax(ax, figsize)
    df = _cta_evidence_frame(genes)
    counts = df["ms_restriction"].value_counts().reindex(MS_RESTRICTION_VALUES, fill_value=0)
    counts = counts[counts > 0]

    colors = [MS_RESTRICTION_COLORS[k] for k in counts.index]
    y = np.arange(len(counts))
    ax.barh(y, counts.values, color=colors, height=0.62)
    ax.set_yticks(y)
    ax.set_yticklabels([_label(c) for c in counts.index])
    ax.invert_yaxis()
    ax.set_xlabel(f"CTAs (n={len(df)})")
    for yi, v in zip(y, counts.values):
        ax.text(v + max(counts.values) * 0.015, yi, str(int(v)), va="center", fontsize=8)
    fig.tight_layout()
    return fig, ax


# ── Figure 2: does MS evidence corroborate HPA restriction confidence? ─────


def plot_confidence_vs_ms_restriction(genes=None, ax=None, figsize=(5.0, 3.2)):
    """Stacked bar: HPA restriction confidence x MS restriction category.

    HPA confidence (HIGH/MODERATE/LOW) comes entirely from oncoref's
    normal-tissue RNA/protein data. This asks a genuinely independent
    question mass spec can answer: within each confidence tier, how much of
    the panel actually has cancer-specific MS support versus a healthy-tissue
    hit HPA's normal-tissue call would not have caught (a different modality,
    same organism, same peptide-MHC axis).
    """
    fig, ax = _new_ax(ax, figsize)
    df = _cta_evidence_frame(genes)
    order = ["HIGH", "MODERATE", "LOW"]
    present = [c for c in order if c in df["restriction_confidence"].unique()]
    cross = pd.crosstab(df["restriction_confidence"], df["ms_restriction"]).reindex(
        index=present, columns=MS_RESTRICTION_VALUES, fill_value=0
    )

    x = np.arange(len(present))
    bottom = np.zeros(len(present))
    for cat in MS_RESTRICTION_VALUES:
        if cat not in cross.columns or cross[cat].sum() == 0:
            continue
        vals = cross[cat].values.astype(float)
        ax.bar(
            x,
            vals,
            bottom=bottom,
            width=0.6,
            color=MS_RESTRICTION_COLORS[cat],
            label=_label(cat),
            edgecolor="white",
            linewidth=1.2,
        )
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels([_label(c) for c in present])
    ax.set_ylabel("CTAs")
    ax.set_xlabel("HPA restriction confidence")
    ax.legend(loc="upper right", fontsize=7, ncol=1)
    fig.tight_layout()
    return fig, ax


# ── Figure 3: literature support depth ──────────────────────────────────


def plot_literature_support_distribution(genes=None, ax=None, figsize=(4.2, 3.0)):
    """Histogram of independent-study count (distinct MS PMIDs) per CTA.

    Answers: how deep is the evidence base per antigen. Long-tailed by
    construction -- a handful of well-studied antigens (NY-ESO-1, MAGE
    family) accumulate dozens of supporting studies while most of the panel
    has a handful or none.
    """
    fig, ax = _new_ax(ax, figsize)
    df = _cta_evidence_frame(genes)
    counts = _pmid_count(df["ms_pmids"])

    max_n = int(counts.max()) if len(counts) else 0
    bins = list(range(0, max_n + 2))
    ax.hist(counts, bins=bins, color=_CATEGORICAL_HUES[0], edgecolor="white", linewidth=0.6)
    ax.set_yscale("log")
    ax.set_xlabel("Distinct supporting studies (PMIDs)")
    ax.set_ylabel("CTAs")
    fig.tight_layout()
    return fig, ax


# ── Figure 4: best-evidenced cancer-specific CTAs ───────────────────────


def plot_top_evidenced_ctas(genes=None, top_n=20, ax=None, figsize=(4.2, 5.0)):
    """Ranked bar of the most literature-supported cancer-only CTAs.

    Restricted to ``CANCER_ONLY`` restriction -- the direct answer to "which
    antigens have the deepest, cleanest public mass-spec support" for
    target-selection purposes.
    """
    fig, ax = _new_ax(ax, figsize)
    df = _cta_evidence_frame(genes)
    df = df[df["ms_restriction"] == "CANCER_ONLY"].copy()
    df["pmid_count"] = _pmid_count(df["ms_pmids"])
    df = df.sort_values("pmid_count", ascending=False).head(top_n)
    df = df.iloc[::-1]  # largest at top

    y = np.arange(len(df))
    ax.barh(y, df["pmid_count"].values, color=MS_RESTRICTION_COLORS["CANCER_ONLY"], height=0.65)
    ax.set_yticks(y)
    ax.set_yticklabels(df["Symbol"].values, fontsize=8)
    ax.set_xlabel("Distinct supporting studies (PMIDs)")
    fig.tight_layout()
    return fig, ax


# ── Figure 5: off-target somatic tissue frequency ───────────────────────


def plot_offtarget_tissue_frequency(genes=None, top_n=15, ax=None, figsize=(4.2, 3.6)):
    """Bar of which healthy somatic tissues most often carry an MS hit.

    Built from every CTA with >=1 healthy-somatic MS observation
    (``SINGLETON_HEALTHY`` or ``RECURRENT_HEALTHY``). Answers a safety
    question the restriction category alone does not: when a CTA does show
    up on normal tissue, where does that risk concentrate across the panel.
    """
    fig, ax = _new_ax(ax, figsize)
    df = _cta_evidence_frame(genes)
    tissues = _somatic_tissues(df["ms_healthy_somatic_tissues"])
    counter: Counter[str] = Counter()
    for row in tissues:
        for t in row.split(";"):
            t = t.strip()
            if t:
                counter[t] += 1
    top = counter.most_common(top_n)
    top = top[::-1]
    labels = [_label(t) for t, _ in top]
    values = [c for _, c in top]

    y = np.arange(len(labels))
    ax.barh(y, values, color=_CATEGORICAL_HUES[1], height=0.65)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("CTAs with an MS hit on this tissue")
    fig.tight_layout()
    return fig, ax


# ── Driver ───────────────────────────────────────────────────────────────

_FIGURE_BUILDERS = {
    "ms-restriction-distribution": plot_ms_restriction_distribution,
    "confidence-vs-ms-restriction": plot_confidence_vs_ms_restriction,
    "literature-support-distribution": plot_literature_support_distribution,
    "top-evidenced-ctas": plot_top_evidenced_ctas,
    "offtarget-tissue-frequency": plot_offtarget_tissue_frequency,
}


#: Figure family name -- namespaces this module's output under
#: ``figures/<name>/`` alongside any other figure family tsarina grows.
FIGURE_FAMILY = "ms-evidence-cta"


def _default_run_dir() -> Path:
    """``figures/<family>/run_<YYYYMMDD-HHMMSS>/``, with a ``latest`` symlink.

    Matches the timestamped-run convention used for oncoref's figure output
    (``oncoref.figure_output``) so a mixed multi-repo figures/ tree reads as
    one system rather than two ad-hoc layouts. Not a dependency on that
    module -- this is a small, self-contained equivalent, since it isn't
    released yet.
    """
    import time

    root = Path("figures") / FIGURE_FAMILY
    stamp = time.strftime("run_%Y%m%d-%H%M%S")
    run_dir = root / stamp
    suffix = 1
    while run_dir.exists():
        suffix += 1
        run_dir = root / f"{stamp}-{suffix}"
    run_dir.mkdir(parents=True)

    latest = root / "latest"
    if latest.is_symlink() or latest.exists():
        latest.unlink()
    latest.symlink_to(run_dir.name)
    return run_dir


def render_all(out_dir: str | Path | None = None, genes=None) -> dict[str, Path]:
    """Render every CTA MS-evidence figure as 300 dpi PNGs.

    Parameters
    ----------
    out_dir
        Where to write the PNGs. ``None`` (the default) creates a fresh
        timestamped run directory under ``figures/ms-evidence-cta/`` and
        repoints the ``latest`` symlink at it. An explicit path is used
        verbatim -- no timestamping, no symlink.
    genes
        Restrict the figures to this gene set; ``None`` uses oncoref's
        default expressed CTA set.

    Returns
    -------
    dict[str, Path]
        Figure name -> path written, so a caller (a CLI, a script, a test)
        can assert on exactly what was produced.
    """
    out_dir = _default_run_dir() if out_dir is None else Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for name, builder in _FIGURE_BUILDERS.items():
        fig, _ = builder(genes=genes)
        path = out_dir / f"{name}.png"
        fig.savefig(path)
        import matplotlib.pyplot as plt

        plt.close(fig)
        written[name] = path
    return written
