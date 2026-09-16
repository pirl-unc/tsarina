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

"""Public MS evidence overlap for mutation-derived (neoantigen) vaccine peptides.

CTAs (:mod:`tsarina.ms_evidence_plots`) are shared, population-level tumor
antigens: the same peptide is meaningful across patients, so its public MS
evidence answers "how well is this target characterized in general". A
neoantigen vaccine peptide is the opposite -- it is built around one
patient's private somatic mutation, so the peptide itself is essentially
never in a public corpus. What *can* be in a public corpus is a window
overlapping the peptide's wild-type-side flanking sequence, or -- less
expected, and worth flagging when it happens -- a window that overlaps the
mutation site itself.

This module takes a "vaccine overlap summary" table (see
``tests/data/osteosarc/README.md`` for the real fixture this was built
against and its provenance/licensing), enumerates every 9-11mer window of
every candidate peptide, checks each against hitlist's public IEDB/CEDAR
observations, and produces three figures:

- :func:`plot_sequence_overlay` -- one row per (gene, peptide) with any
  public hit, the peptide drawn as text with the minimal epitope highlighted
  and MS-evidence spans colored by tissue category and outlined in red
  wherever they overlap the epitope itself.
- :func:`plot_tissue_provenance` -- a per-tissue breakdown for whichever
  windows have enough hits that a bare count doesn't tell the real story.
- :func:`plot_ms_hit_ranking` -- every gene, including the ones with zero
  public evidence, ranked by total public MS observations.

Deliberately not library code for *production* use the way
``ms_evidence_plots`` is: there is no general-population neoantigen
category to characterize the way there is for CTAs. This exists to make
the fixture data in ``tests/data/osteosarc/`` exercisable and visualizable,
for a specific patient's specific vaccine constructs.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

#: The "vaccine overlap summary" wide-format columns, mapped to a short
#: construct name. Order matters only for readability of output frames.
_CONSTRUCT_COLUMNS = {
    "mRNA minimal epitope": "mRNA_minimal",
    "Peptide from JLF for elispot": "JLF_elispot",
    "mRNA full peptide sequence": "mRNA_full",
    "JLF V3 vaccine Full Peptide": "JLF_V3",
    "JLF V2 vaccine Full Peptide": "JLF_V2",
    "JLF V1 vaccine Full Peptide": "JLF_V1",
    "CeGaT vaccine Full Peptide": "CeGaT",
}

#: Construct names treated as candidates for "the minimal epitope" when
#: looking for the mutation-bearing core of a longer flanking sequence.
_MINIMAL_CONSTRUCTS = frozenset({"mRNA_minimal", "JLF_elispot"})

_NOT_A_SEQUENCE = frozenset({"", "NA", "NAN", "?"})


def load_vaccine_peptide_table(path: str | Path) -> pd.DataFrame:
    """Parse a wide-format vaccine-overlap-summary TSV into a tidy frame.

    Returns
    -------
    pd.DataFrame
        Columns ``gene``, ``mutation``, ``construct``, ``peptide`` -- one
        row per non-empty peptide-sequence cell. The elispot column's
        comma-separated multi-peptide cells are split into one row each.
        Deduplicated on (gene, peptide): several constructs commonly carry
        the identical sequence.
    """
    raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    rows = []
    for _, r in raw.iterrows():
        gene = str(r.get("Gene Target", "")).strip()
        if not gene:
            continue
        mutation = str(r.get("Mutation: HGVSp", "")).strip()
        for column, construct in _CONSTRUCT_COLUMNS.items():
            if column not in raw.columns:
                continue
            cell = str(r.get(column, "")).strip()
            for candidate in cell.split(","):
                pep = candidate.strip().upper()
                if pep and pep not in _NOT_A_SEQUENCE and pep.isalpha():
                    rows.append((gene, mutation, construct, pep))
    df = pd.DataFrame(rows, columns=["gene", "mutation", "construct", "peptide"])
    return df.drop_duplicates(subset=["gene", "peptide"]).reset_index(drop=True)


def _minimal_epitope_for_gene(df: pd.DataFrame, gene: str) -> str | None:
    """The gene's shortest minimal/elispot-construct peptide, if any."""
    sub = df[(df["gene"] == gene) & (df["construct"].isin(_MINIMAL_CONSTRUCTS))]
    if sub.empty:
        return None
    return sub.loc[sub["peptide"].str.len().idxmin(), "peptide"]


def ms_evidence_for_peptides(
    df: pd.DataFrame,
    window_lengths: tuple[int, ...] = (9, 10, 11),
) -> tuple[pd.DataFrame, dict[str, list[tuple[str, str, int, int]]]]:
    """Check every window of every peptide against hitlist's public MS index.

    Parameters
    ----------
    df
        A tidy peptide table, as returned by :func:`load_vaccine_peptide_table`.
    window_lengths
        Contiguous window lengths to enumerate. Defaults to the dominant
        MHC class I lengths.

    Returns
    -------
    hits : pd.DataFrame
        One row per (window, public observation) match, with tissue/source
        annotation columns.
    window_source : dict[str, list[(gene, peptide, start, length)]]
        Every window that appears in ``hits``, mapped to every (gene,
        peptide, position) it was generated from -- a window can come from
        more than one peptide if two sequences share a sub-region.
    """
    from hitlist.observations import load_observations

    windows: set[str] = set()
    window_source: dict[str, list[tuple[str, str, int, int]]] = {}
    for _, row in df.iterrows():
        seq = row["peptide"]
        for length in window_lengths:
            if len(seq) < length:
                continue
            for start in range(len(seq) - length + 1):
                w = seq[start : start + length]
                windows.add(w)
                window_source.setdefault(w, []).append((row["gene"], seq, start, length))

    if not windows:
        return pd.DataFrame(), {}

    hits = load_observations(
        peptide=sorted(windows),
        columns=[
            "peptide",
            "mhc_restriction",
            "source_tissue",
            "disease",
            "src_cancer",
            "src_healthy_tissue",
            "src_healthy_reproductive",
            "src_healthy_thymus",
            "cell_name",
            "pmid",
        ],
    )
    matched = {w for w in window_source if w in set(hits["peptide"])}
    return hits, {w: window_source[w] for w in matched}


def _categorize_hit(row) -> str:
    """cancer / healthy / mixed / other, from hitlist's per-row curation."""
    disease = str(row.get("disease") or "").lower()
    cancer = bool(row.get("src_cancer"))
    healthy = bool(row.get("src_healthy_tissue")) or disease == "healthy"
    if cancer and not healthy:
        return "cancer"
    if healthy and not cancer:
        return "healthy"
    if cancer and healthy:
        return "mixed"
    return "other"


CATEGORY_COLOR = {
    "cancer": "#eb6834",
    "healthy": "#2a78d6",
    "mixed": "#4a3aa7",
    "other": "#9a9a94",
}
CATEGORY_ORDER = ("cancer", "healthy", "mixed", "other")
_TEXT_COLOR = "#3a3a38"


def _apply_publication_style() -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.transparent": False,
            "savefig.dpi": 300,
            "font.size": 9,
            "text.color": _TEXT_COLOR,
            "axes.edgecolor": _TEXT_COLOR,
            "axes.labelcolor": _TEXT_COLOR,
            "xtick.color": _TEXT_COLOR,
            "ytick.color": _TEXT_COLOR,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "pdf.fonttype": 42,
        }
    )


def _window_summary(hits: pd.DataFrame) -> dict[str, dict]:
    """Per matched window: observation count, category counts, tissues."""
    if hits.empty:
        return {}
    hits = hits.copy()
    hits["category"] = hits.apply(_categorize_hit, axis=1)
    out = {}
    for pep, grp in hits.groupby("peptide"):
        tissue_values = grp["source_tissue"].astype("string").fillna("")
        tissues = sorted(set(tissue_values.replace("", "(unspecified)")))
        out[pep] = {
            "n": len(grp),
            "categories": grp["category"].value_counts().to_dict(),
            "tissues": tissues,
            "pmids": sorted(set(grp["pmid"].dropna().astype(str))),
        }
    return out


def _gene_hit_layout(df: pd.DataFrame, window_source: dict, window_info: dict) -> list[dict]:
    """One entry per (gene, peptide) that carries >=1 matched window.

    Each entry: gene, peptide, epitope_span (or None), hit_spans (list of
    (window, start, end, dominant_category, overlaps_epitope)).
    """
    # window -> (gene, peptide) -> (start, length) it occurs at in that
    # specific peptide (a window can be a substring of more than one
    # peptide; each occurrence gets its own row).
    by_gene_peptide: dict[tuple[str, str], list[tuple[str, int, int]]] = {}
    for window, occurrences in window_source.items():
        for gene, peptide, start, length in occurrences:
            by_gene_peptide.setdefault((gene, peptide), []).append((window, start, start + length))

    rows = []
    for (gene, peptide), spans in sorted(by_gene_peptide.items()):
        epitope = _minimal_epitope_for_gene(df, gene)
        epitope_span = None
        if epitope and epitope != peptide:
            pos = peptide.find(epitope)
            if pos >= 0:
                epitope_span = (pos, pos + len(epitope))
        elif epitope == peptide:
            # The peptide *is* the minimal epitope -- nothing shorter to bracket.
            epitope_span = None

        hit_spans = []
        for window, s, e in spans:
            info = window_info[window]
            dominant = max(info["categories"], key=info["categories"].get)
            overlaps = False
            if epitope_span is not None:
                e0, e1 = epitope_span
                overlaps = not (e1 <= s or e <= e0)
            hit_spans.append((window, s, e, dominant, overlaps))

        rows.append(
            {
                "gene": gene,
                "peptide": peptide,
                "epitope_span": epitope_span,
                "hit_spans": hit_spans,
            }
        )
    return rows


#: Pale, low-chroma fill for the epitope highlight -- a highlighter-style
#: wash behind the letters, not a bordered box competing with the
#: MS-evidence bars for visual weight.
_EPITOPE_HIGHLIGHT = "#fbe6a8"


def plot_sequence_overlay(
    layout: list[dict],
    window_info: dict,
    ax=None,
    figsize: tuple[float, float] | None = None,
):
    """One row per (gene, peptide) with public MS evidence.

    A pale highlight behind the letters marks the vaccinated (mutant)
    minimal epitope, when known. A colored bar beneath the sequence marks
    each public MS-evidence span, colored by dominant tissue category; a
    span outlined in red, rather than plain, overlaps the epitope itself
    instead of sitting purely in flanking sequence -- the distinction that
    actually matters for interpreting the hit.
    """
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    _apply_publication_style()
    order = sorted(layout, key=lambda r: -sum(window_info[w]["n"] for w, *_ in r["hit_spans"]))
    n = len(order)
    max_len = max((len(r["peptide"]) for r in order), default=1)
    row_h = 1.5
    if figsize is None:
        figsize = (max(10, max_len * 0.32), n * row_h * 0.62 + 1.2)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    for i, row in enumerate(order):
        y = (n - i) * row_h
        gene, seq = row["gene"], row["peptide"]
        ax.text(
            -2, y, gene, ha="right", va="center", fontsize=9, fontweight="bold", family="sans-serif"
        )
        if row["epitope_span"]:
            e0, e1 = row["epitope_span"]
            ax.add_patch(
                plt.Rectangle(
                    (e0 - 0.5, y - 0.34),
                    (e1 - e0),
                    0.68,
                    facecolor=_EPITOPE_HIGHLIGHT,
                    edgecolor="none",
                    zorder=1,
                )
            )
        for pos, aa in enumerate(seq):
            ax.text(
                pos, y, aa, ha="center", va="center", fontsize=8.5, family="monospace", zorder=2
            )
        notes = []
        for window, s, e, category, overlaps in row["hit_spans"]:
            info = window_info[window]
            ax.add_patch(
                plt.Rectangle(
                    (s - 0.5, y - 0.62),
                    (e - s),
                    0.16,
                    facecolor=CATEGORY_COLOR[category],
                    edgecolor="#b3261e" if overlaps else "none",
                    linewidth=1.1 if overlaps else 0,
                    zorder=2,
                )
            )
            flag = " [OVERLAPS MUTATION]" if overlaps else ""
            n_tissue = len(info["tissues"])
            tissue_note = (
                f"n={info['n']}, {n_tissue}t" if n_tissue > 2 else ", ".join(info["tissues"])
            )
            notes.append(f"{window} ({tissue_note}){flag}")
        ax.text(
            len(seq) + 1,
            y,
            "  |  ".join(notes),
            ha="left",
            va="center",
            fontsize=6.8,
            family="sans-serif",
            color=_TEXT_COLOR,
        )

    ax.set_xlim(-9, max_len + 26)
    ax.set_ylim(0, (n + 1) * row_h)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    legend_handles = [
        mpatches.Patch(color=_EPITOPE_HIGHLIGHT, label="vaccinated (mutant) epitope"),
        mpatches.Patch(color=CATEGORY_COLOR["cancer"], label="cancer / tumor cell line"),
        mpatches.Patch(color=CATEGORY_COLOR["healthy"], label="healthy tissue"),
        mpatches.Patch(color=CATEGORY_COLOR["mixed"], label="both, across studies"),
        mpatches.Patch(color=CATEGORY_COLOR["other"], label="ambiguous curation"),
        mpatches.Patch(
            facecolor="white", edgecolor="#b3261e", linewidth=1.4, label="hit overlaps mutation"
        ),
    ]
    ax.legend(
        handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, 1.05), ncol=3, fontsize=7
    )
    fig.tight_layout()
    return fig, ax


def plot_tissue_provenance(hits: pd.DataFrame, windows: list[str], ax=None, figsize=(11, 4.6)):
    """Per-tissue, per-category breakdown for the given windows.

    Meant for whichever windows have enough observations that
    :func:`plot_sequence_overlay`'s compact text note doesn't tell the
    real story -- pass the busiest one or two windows, not all of them.
    """
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    _apply_publication_style()
    hits = hits.copy()
    hits["category"] = hits.apply(_categorize_hit, axis=1)
    hits["source_tissue"] = (
        hits["source_tissue"].astype("string").fillna("").replace("", "(unspecified)")
    )

    fig, axes = plt.subplots(1, len(windows), figsize=figsize, squeeze=False)
    axes = axes[0]
    for ax, window in zip(axes, windows):
        sub = hits[hits["peptide"] == window]
        tt = sub.groupby(["source_tissue", "category"], observed=True).size().unstack(fill_value=0)
        tt = tt.reindex(columns=list(CATEGORY_ORDER), fill_value=0)
        tt = tt.loc[tt.sum(axis=1).sort_values(ascending=True).index]
        y = range(len(tt))
        left = np.zeros(len(tt))
        for cat in CATEGORY_ORDER:
            vals = tt[cat].to_numpy()
            if vals.sum() == 0:
                continue
            ax.barh(list(y), vals, left=left, color=CATEGORY_COLOR[cat], height=0.65)
            left = left + vals
        ax.set_yticks(list(y))
        ax.set_yticklabels(tt.index, fontsize=8)
        ax.set_xlabel("public MS observations")
        n_pmid = sub["pmid"].dropna().nunique()
        ax.set_title(f"{window}  (n={len(sub)}, {n_pmid} PMIDs)", fontsize=9)

    legend_handles = [mpatches.Patch(color=CATEGORY_COLOR[c], label=c) for c in CATEGORY_ORDER]
    fig.legend(
        handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, 1.06), ncol=4, fontsize=8
    )
    fig.tight_layout()
    return fig, axes


def plot_ms_hit_ranking(
    all_genes: list[str],
    layout: list[dict],
    window_info: dict,
    ax=None,
    figsize=(7, 9),
):
    """Every candidate gene, ranked by total public MS observations.

    Includes genes with zero hits, so the contrast is visible. Each bar is
    colored by the strongest mutation-overlap finding for that gene:
    a public hit overlapping the vaccinated epitope itself outranks a
    flanking-only hit for coloring purposes, since it is the more
    consequential finding to see at a glance.
    """
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    _apply_publication_style()
    STATUS_COLOR = {
        "overlaps_mutation": "#e34948",  # red -- the one finding worth double-checking
        "flanking_only": CATEGORY_COLOR["cancer"],
        "unknown_position": "#9a9a94",
        "no_evidence": "#d8d7d0",
    }
    STATUS_LABEL = {
        "overlaps_mutation": "public hit overlaps the mutation",
        "flanking_only": "public hit, flanking sequence only",
        "unknown_position": "public hit, epitope position unknown",
        "no_evidence": "no public MS evidence",
    }

    by_gene = {}
    for row in layout:
        by_gene.setdefault(row["gene"], []).append(row)

    ranking = []
    for gene in all_genes:
        rows = by_gene.get(gene, [])
        total = sum(window_info[w]["n"] for r in rows for w, *_ in r["hit_spans"])
        if not rows:
            status = "no_evidence"
        else:
            has_epitope_info = any(r["epitope_span"] is not None for r in rows)
            overlaps = any(ov for r in rows for *_, ov in r["hit_spans"])
            if overlaps:
                status = "overlaps_mutation"
            elif has_epitope_info:
                status = "flanking_only"
            else:
                status = "unknown_position"
        ranking.append((gene, total, status))
    ranking.sort(key=lambda r: -r[1])

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    y = range(len(ranking))
    colors = [STATUS_COLOR[status] for _, _, status in ranking]
    totals = [t for _, t, _ in ranking]
    ax.barh(list(y), totals, color=colors, height=0.7)
    ax.set_yticks(list(y))
    ax.set_yticklabels([g for g, _, _ in ranking], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("total public MS observations (all windows, all constructs)")
    for yi, total in zip(y, totals):
        if total > 0:
            ax.text(total + max(totals) * 0.01, yi, str(total), va="center", fontsize=7)

    legend_handles = [
        mpatches.Patch(color=STATUS_COLOR[s], label=STATUS_LABEL[s]) for s in STATUS_COLOR
    ]
    ax.legend(handles=legend_handles, loc="lower right", fontsize=7)
    fig.tight_layout()
    return fig, ax


def render_all(
    peptide_table_path: str | Path,
    all_gene_names: list[str] | None = None,
    out_dir: str | Path | None = None,
) -> dict[str, Path]:
    """Full pipeline: load table -> public MS lookup -> all three figures.

    Parameters
    ----------
    peptide_table_path
        Path to a wide-format vaccine-overlap-summary TSV.
    all_gene_names
        Every candidate gene name, including ones with no public MS
        evidence, for the ranking figure. Defaults to every gene in the
        table.
    out_dir
        Where to write PNGs. ``None`` creates a fresh timestamped run
        directory under ``figures/ms-evidence-neoantigen-vaccine/`` with a
        ``latest`` symlink, matching :mod:`tsarina.ms_evidence_plots`. An
        explicit path is used verbatim.
    """
    import matplotlib.pyplot as plt

    df = load_vaccine_peptide_table(peptide_table_path)
    if all_gene_names is None:
        all_gene_names = sorted(df["gene"].unique())

    hits, window_source = ms_evidence_for_peptides(df)
    window_info = _window_summary(hits)
    layout = _gene_hit_layout(df, window_source, window_info)

    if out_dir is None:
        out_dir = _default_run_dir()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    written = {}

    fig, _ = plot_sequence_overlay(layout, window_info)
    path = out_dir / "sequence-overlay.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    written["sequence-overlay"] = path

    if window_info:
        busiest = sorted(window_info, key=lambda w: -window_info[w]["n"])[:2]
        fig, _ = plot_tissue_provenance(hits, busiest)
        path = out_dir / "tissue-provenance.png"
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
        written["tissue-provenance"] = path

    fig, _ = plot_ms_hit_ranking(all_gene_names, layout, window_info)
    path = out_dir / "ms-hit-ranking.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    written["ms-hit-ranking"] = path

    return written


FIGURE_FAMILY = "ms-evidence-neoantigen-vaccine"


def _default_run_dir() -> Path:
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
