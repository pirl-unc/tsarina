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

"""Population-spanning pMHC sets for off-the-shelf vaccine / TCR programs.

A *spanning set* is a small panel of peptide-HLA pairs designed to
maximize patient coverage across a population: the most common cancer-
testis antigens (CTAs) crossed with the most frequent HLA class I
alleles, with each (CTA, HLA) cell filled by the single best-presented
peptide that the CTA's protein offers for that allele.

The output is a wide pivot table — CTA rows, HLA columns, peptide at the
intersection — plus an optional long format that preserves percentile,
presentation score, and affinity for downstream filtering.

Typical usage::

    from tsarina.spanning import spanning_pmhc_set

    table = spanning_pmhc_set(
        cta_count=25,
        panel="iedb27_ab",
        max_percentile=2.0,
    )
"""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

_DEFAULT_RANK_COLUMN = "ms_cancer_peptide_count"

_LONG_OUTPUT_COLUMNS: tuple[str, ...] = (
    "cta",
    "allele",
    "peptide",
    "length",
    "presentation_percentile",
    "presentation_score",
    "affinity_nm",
)


def spanning_pmhc_set(
    cta_count: int = 25,
    cta_rank_by: str = _DEFAULT_RANK_COLUMN,
    ctas: Iterable[str] | None = None,
    min_restriction_confidence: Iterable[str] | None = ("HIGH", "MODERATE"),
    restriction_levels: Iterable[str] | None = None,
    alleles: Iterable[str] | None = None,
    panel: str | None = "iedb27_ab",
    lengths: tuple[int, ...] = (9,),
    ensembl_release: int = 112,
    require_cta_exclusive: bool = True,
    predictor: str = "mhcflurry",
    max_percentile: float = 2.0,
    output_format: str = "wide",
) -> pd.DataFrame:
    """Build a CTA x HLA spanning pMHC table.

    Parameters
    ----------
    cta_count
        Number of top CTAs to include when ``ctas`` is not supplied.
    cta_rank_by
        Column in the bundled CTA CSV used to rank candidates.  Default
        ``"ms_cancer_peptide_count"`` (most-observed-in-cancer first).
        Falls back to alphabetical ``Symbol`` order when the column is
        missing or all-NaN.
    ctas
        Explicit list of gene symbols.  Overrides ``cta_count`` /
        ``cta_rank_by``.  Genes not present in the bundled CTA CSV are
        silently dropped.
    min_restriction_confidence
        Allowed ``restriction_confidence`` bins (e.g. ``("HIGH",
        "MODERATE")``).  Pass ``None`` to disable.  Applied before
        ranking.
    restriction_levels
        Optional set of ``restriction`` values to keep (e.g.
        ``("TESTIS", "PLACENTAL")``).  ``None`` (default) keeps all.
    alleles
        Explicit allele list.  Overrides ``panel``.
    panel
        Named panel from :mod:`tsarina.alleles`.  Default
        ``"iedb27_ab"`` (27 globally-common HLA-A/B alleles).
    lengths
        Peptide lengths (default 9-mers — the dominant class I length and
        the most therapeutically relevant for vaccine / TCR design).
    ensembl_release
        Ensembl release for peptide enumeration (default 112).
    require_cta_exclusive
        If True (default), use :func:`tsarina.peptides.cta_exclusive_peptides`
        — peptides that do NOT appear in any non-CTA protein.  Set False
        to fall back to :func:`tsarina.peptides.cta_peptides` (all CTA
        k-mers, no exclusivity gate).
    predictor
        mhctools predictor key.  ``"mhcflurry"`` (default), ``"netmhcpan"``,
        or ``"netmhcpan_el"``.
    max_percentile
        Maximum presentation percentile for a cell to be filled.  Cells
        whose best peptide exceeds this threshold are left blank in wide
        output (or absent in long output).  Default 2.0% — the standard
        MHC class I "presented" threshold.
    output_format
        ``"wide"`` (default) — pivot with CTA rows, allele columns, and
        the chosen peptide as the cell value.
        ``"long"`` — one row per filled cell with peptide, length,
        percentile, score, and affinity columns.

    Returns
    -------
    pd.DataFrame
        Wide pivot or long table per ``output_format``.

    Raises
    ------
    ValueError
        If neither ``alleles`` nor ``panel`` is supplied, or if
        ``output_format`` is unrecognized.
    ImportError
        If the scoring backend (topiary + mhctools + chosen predictor)
        is not installed.
    """
    if output_format not in ("wide", "long"):
        raise ValueError(f"output_format must be 'wide' or 'long', got {output_format!r}")

    allele_list = _resolve_alleles(alleles, panel)
    cta_list = _resolve_ctas(
        ctas=ctas,
        cta_count=cta_count,
        cta_rank_by=cta_rank_by,
        min_restriction_confidence=min_restriction_confidence,
        restriction_levels=restriction_levels,
    )

    if not cta_list or not allele_list:
        return _empty_output(cta_list, allele_list, output_format)

    pep_df = _resolve_peptides(
        cta_list=cta_list,
        lengths=lengths,
        ensembl_release=ensembl_release,
        require_cta_exclusive=require_cta_exclusive,
    )
    if pep_df.empty:
        return _empty_output(cta_list, allele_list, output_format)

    from .scoring import score_presentation

    unique_peptides = pep_df["peptide"].unique().tolist()
    scores = score_presentation(peptides=unique_peptides, alleles=allele_list, predictor=predictor)

    best = _best_per_cell(scores, pep_df, max_percentile)

    if output_format == "wide":
        return _to_wide(best, cta_list, allele_list)
    return _to_long(best, cta_list, allele_list)


# ── Helpers ────────────────────────────────────────────────────────────


def _resolve_alleles(
    alleles: Iterable[str] | None,
    panel: str | None,
) -> list[str]:
    if alleles is not None:
        out = [a for a in alleles if a]
        if not out:
            raise ValueError("alleles is empty")
        return out
    if panel is None:
        raise ValueError("Specify either alleles or panel")
    from .alleles import get_panel

    return list(get_panel(panel))


def _resolve_ctas(
    ctas: Iterable[str] | None,
    cta_count: int,
    cta_rank_by: str,
    min_restriction_confidence: Iterable[str] | None,
    restriction_levels: Iterable[str] | None,
) -> list[str]:
    """Pick the CTA gene set per the resolution order: explicit override
    wins; else filter the bundled CTA CSV and take top-N by rank."""
    from .gene_sets import CTA_gene_names
    from .loader import cta_dataframe

    valid = CTA_gene_names()

    if ctas is not None:
        return [g for g in ctas if g in valid]

    df = cta_dataframe()
    if "filtered" in df.columns:
        df = df[df["filtered"].astype(str).str.lower() == "true"]
    if "never_expressed" in df.columns:
        df = df[~(df["never_expressed"].astype(str).str.lower() == "true")]

    if min_restriction_confidence is not None and "restriction_confidence" in df.columns:
        wanted = {c.upper() for c in min_restriction_confidence}
        df = df[df["restriction_confidence"].astype(str).str.upper().isin(wanted)]

    if restriction_levels is not None and "restriction" in df.columns:
        wanted_levels = set(restriction_levels)
        df = df[df["restriction"].isin(wanted_levels)]

    if cta_rank_by in df.columns and df[cta_rank_by].notna().any():
        df = df.sort_values(cta_rank_by, ascending=False, na_position="last")
    else:
        df = df.sort_values("Symbol")

    df = df[df["Symbol"].isin(valid)]
    return df["Symbol"].head(cta_count).tolist()


def _resolve_peptides(
    cta_list: list[str],
    lengths: tuple[int, ...],
    ensembl_release: int,
    require_cta_exclusive: bool,
) -> pd.DataFrame:
    if require_cta_exclusive:
        from .peptides import cta_exclusive_peptides

        all_peps = cta_exclusive_peptides(ensembl_release=ensembl_release, lengths=tuple(lengths))
    else:
        from .peptides import cta_peptides

        all_peps = cta_peptides(ensembl_release=ensembl_release, lengths=tuple(lengths))
    if all_peps.empty:
        return all_peps
    return all_peps[all_peps["gene_name"].isin(cta_list)].copy()


def _best_per_cell(
    scores: pd.DataFrame,
    pep_df: pd.DataFrame,
    max_percentile: float,
) -> pd.DataFrame:
    """For each (cta, allele), pick the lowest-percentile peptide above
    the max_percentile cutoff."""
    if scores.empty or "presentation_percentile" not in scores.columns:
        return pd.DataFrame(
            columns=[
                "cta",
                "allele",
                "peptide",
                "length",
                "presentation_percentile",
                "presentation_score",
                "affinity_nm",
            ]
        )

    pep_meta = (
        pep_df[["peptide", "gene_name", "length"]]
        .drop_duplicates(subset="peptide")
        .rename(columns={"gene_name": "cta"})
    )
    scored = scores.merge(pep_meta, on="peptide", how="inner")

    sort_cols = ["presentation_percentile", "peptide", "allele"]
    best = (
        scored.sort_values(sort_cols, kind="mergesort")
        .groupby(["cta", "allele"], as_index=False)
        .first()
    )
    best = best[best["presentation_percentile"] <= max_percentile].copy()
    for col in ("presentation_score", "affinity_nm"):
        if col not in best.columns:
            best[col] = pd.NA
    return best[
        [
            "cta",
            "allele",
            "peptide",
            "length",
            "presentation_percentile",
            "presentation_score",
            "affinity_nm",
        ]
    ]


def _to_wide(
    best: pd.DataFrame,
    cta_list: list[str],
    allele_list: list[str],
) -> pd.DataFrame:
    if best.empty:
        wide = pd.DataFrame(index=cta_list, columns=allele_list, dtype=object)
    else:
        wide = best.pivot(index="cta", columns="allele", values="peptide")
        wide = wide.reindex(index=cta_list, columns=allele_list)
    wide.index.name = "cta"
    wide = wide.reset_index()
    return wide


def _to_long(
    best: pd.DataFrame,
    cta_list: list[str],
    allele_list: list[str],
) -> pd.DataFrame:
    if best.empty:
        return pd.DataFrame(columns=list(_LONG_OUTPUT_COLUMNS))
    cta_order = {gene: i for i, gene in enumerate(cta_list)}
    allele_order = {a: i for i, a in enumerate(allele_list)}
    out = best.copy()
    out["_cta_order"] = out["cta"].map(cta_order)
    out["_allele_order"] = out["allele"].map(allele_order)
    out = out.sort_values(["_cta_order", "_allele_order"]).drop(
        columns=["_cta_order", "_allele_order"]
    )
    return out[list(_LONG_OUTPUT_COLUMNS)].reset_index(drop=True)


def _empty_output(cta_list: list[str], allele_list: list[str], output_format: str) -> pd.DataFrame:
    if output_format == "wide":
        wide = pd.DataFrame(index=cta_list, columns=allele_list, dtype=object)
        wide.index.name = "cta"
        return wide.reset_index()
    return pd.DataFrame(columns=list(_LONG_OUTPUT_COLUMNS))
