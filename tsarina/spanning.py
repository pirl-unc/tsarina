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

"""CTA x HLA pMHC panel matrices for off-the-shelf vaccine / TCR programs.

A panel is a matrix of CTA source proteins crossed with a population HLA
allele panel. Each (CTA, HLA) cell is filled by the best peptide-HLA candidate
that satisfies the configured public-MS evidence tier and prediction cutoff.

The default is MS-evidence-first: mono-allelic MS evidence is allowed the
least stringent prediction cutoff, multi-allelic sample-genotype evidence is
stricter, unrestricted MS evidence is stricter still, and prediction-only
candidates are excluded unless explicitly requested.

The output is a wide pivot table — CTA rows, HLA columns, peptide at the
intersection — plus an optional long format that preserves evidence tier,
MS provenance, percentile, presentation score, and affinity for downstream
filtering.

Typical usage::

    from tsarina.spanning import spanning_pmhc_set

    table = spanning_pmhc_set(
        cta_count=25,
        panel="global51_abc_ssa",
        lengths=(8, 9, 10, 11),
    )
"""

from __future__ import annotations

import time
from collections.abc import Callable, Iterable
from pathlib import Path

import pandas as pd

_DEFAULT_RANK_COLUMN = "ms_cancer_peptide_count"
_DEFAULT_PANEL = "global51_abc_ssa"
_DEFAULT_LENGTHS = (8, 9, 10, 11)
_DEFAULT_MONOALLELIC_MS_MAX_PERCENTILE = 2.0
_DEFAULT_SAMPLE_ALLELE_MS_MAX_PERCENTILE = 1.0
_DEFAULT_UNRESTRICTED_MS_MAX_PERCENTILE = 0.5
_DEFAULT_PREDICTED_ONLY_MAX_PERCENTILE = 0.1

_EVIDENCE_TIER_RANKS: dict[str, int] = {
    "monoallelic_ms": 0,
    "sample_allele_ms": 1,
    "unrestricted_ms": 2,
    "predicted_only": 3,
}

_LONG_OUTPUT_COLUMNS: tuple[str, ...] = (
    "cta",
    "allele",
    "peptide",
    "length",
    "evidence_tier",
    "ms_hit_count",
    "ms_alleles",
    "ms_pmids",
    "ms_samples",
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
    panel: str | None = _DEFAULT_PANEL,
    lengths: tuple[int, ...] = _DEFAULT_LENGTHS,
    ensembl_release: int = 112,
    require_cta_exclusive: bool = True,
    predictor: str = "mhcflurry",
    monoallelic_ms_max_percentile: float = _DEFAULT_MONOALLELIC_MS_MAX_PERCENTILE,
    sample_allele_ms_max_percentile: float = _DEFAULT_SAMPLE_ALLELE_MS_MAX_PERCENTILE,
    unrestricted_ms_max_percentile: float = _DEFAULT_UNRESTRICTED_MS_MAX_PERCENTILE,
    include_predicted_only: bool = False,
    predicted_only_max_percentile: float = _DEFAULT_PREDICTED_ONLY_MAX_PERCENTILE,
    max_percentile: float | None = None,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    output_format: str = "wide",
    on_progress: Callable[[str], None] | None = None,
) -> pd.DataFrame:
    """Build a CTA x HLA pMHC panel matrix.

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
        silently dropped.  Note: ``min_restriction_confidence`` and
        ``restriction_levels`` gates do **not** apply on this path —
        the explicit list wins verbatim.
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
        ``"global51_abc_ssa"`` (51 globally broad HLA-A/B/C alleles).
    lengths
        Peptide lengths (default 8-11-mers for MHC class I).
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
    monoallelic_ms_max_percentile
        Maximum percentile for a peptide-HLA candidate with mono-allelic
        MS support for that allele. Default 2.0.
    sample_allele_ms_max_percentile
        Maximum percentile for a candidate observed in a multi-allelic sample
        whose genotype contains the HLA allele and whose predicted percentile
        is best among that sample's alleles. Default 1.0.
    unrestricted_ms_max_percentile
        Maximum percentile for a candidate whose peptide has MS evidence but
        no usable allele assignment. Default 0.5.
    include_predicted_only
        If True, allow prediction-only candidates. They are ranked below all
        MS-supported candidates and therefore only fill otherwise-empty cells.
        Default False.
    predicted_only_max_percentile
        Maximum percentile for prediction-only candidates when enabled.
        Default 0.1.
    max_percentile
        Deprecated compatibility shortcut. When supplied, sets all three
        MS-supported tier cutoffs to the same value. It does not enable or
        relax prediction-only candidates.
    iedb_path, cedar_path
        Optional explicit raw public-MS files. Defaults use hitlist's cached
        observations path.
    output_format
        ``"wide"`` (default) — pivot with CTA rows, allele columns, and
        the chosen peptide as the cell value.
        ``"long"`` — one row per filled cell with peptide, length,
        percentile, score, and affinity columns.
    on_progress
        Optional callable invoked with a human-readable status string at
        each pipeline stage (pre-scoring, post-scoring with elapsed
        seconds, post-pivot summary).  Used by the CLI handler to emit
        stderr progress lines on multi-minute runs; library callers
        leave it ``None`` (default) for silent operation.  The library
        itself never prints — it only formats messages and hands them
        to the callback.

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

    if max_percentile is not None:
        monoallelic_ms_max_percentile = max_percentile
        sample_allele_ms_max_percentile = max_percentile
        unrestricted_ms_max_percentile = max_percentile

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

    from .ms_evidence import load_public_ms_hits

    unique_peptides = pep_df["peptide"].unique().tolist()
    hits = load_public_ms_hits(
        unique_peptides,
        iedb_path=iedb_path,
        cedar_path=cedar_path,
        mhc_class="I",
        mhc_species="Homo sapiens",
    )

    from .scoring import score_presentation

    score_alleles = _score_alleles_for_panel(allele_list, hits)
    n_predictions = len(unique_peptides) * len(score_alleles)
    if on_progress is not None:
        on_progress(
            f"Scoring {len(unique_peptides)} peptides x {len(score_alleles)} alleles "
            f"(~{n_predictions} predictions) via {predictor}..."
        )

    scoring_start = time.perf_counter()
    scores = score_presentation(
        peptides=unique_peptides, alleles=score_alleles, predictor=predictor
    )
    scoring_elapsed = time.perf_counter() - scoring_start

    if on_progress is not None:
        on_progress(f"Scored {n_predictions} predictions in {scoring_elapsed:.1f}s")

    evidence_stats = _build_evidence_stats(hits, allele_list, _score_lookup(scores))
    candidates = _candidate_rows(
        scores=scores,
        pep_df=pep_df,
        allele_list=allele_list,
        evidence_stats=evidence_stats,
        cutoffs={
            "monoallelic_ms": monoallelic_ms_max_percentile,
            "sample_allele_ms": sample_allele_ms_max_percentile,
            "unrestricted_ms": unrestricted_ms_max_percentile,
            "predicted_only": predicted_only_max_percentile,
        },
        include_predicted_only=include_predicted_only,
    )
    best = _best_per_cell(candidates)

    if on_progress is not None:
        on_progress(
            f"Panel: {len(cta_list)} CTAs x {len(allele_list)} alleles, "
            f"{len(best)} filled cells using MS tier cutoffs "
            f"{monoallelic_ms_max_percentile}/{sample_allele_ms_max_percentile}/"
            f"{unrestricted_ms_max_percentile}"
        )

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


def _score_alleles_for_panel(allele_list: list[str], hits: pd.DataFrame) -> list[str]:
    """Return panel alleles plus sample-genotype alleles needed for best-of-sample checks."""
    from .mhc import split_mhc_restrictions

    out = list(dict.fromkeys(allele_list))
    seen = set(out)
    if hits.empty or "mhc_allele_set" not in hits.columns:
        return out

    for _, row in hits.iterrows():
        if str(row.get("mhc_allele_provenance", "")).strip() != "sample_allele_match":
            continue
        value = row.get("mhc_allele_set", "")
        for allele in split_mhc_restrictions(value):
            if allele.startswith("HLA-") and "*" in allele and allele not in seen:
                out.append(allele)
                seen.add(allele)
    return out


def _score_lookup(scores: pd.DataFrame) -> dict[tuple[str, str], float]:
    if scores.empty or "presentation_percentile" not in scores.columns:
        return {}
    return {
        (str(row.peptide), str(row.allele)): float(row.presentation_percentile)
        for row in scores.itertuples(index=False)
        if pd.notna(row.presentation_percentile)
    }


def _new_evidence_bucket() -> dict[str, object]:
    return {
        "ms_hit_count": 0,
        "ms_alleles": set(),
        "ms_pmids": set(),
        "ms_samples": set(),
    }


def _add_evidence(
    stats: dict[tuple[str, str | None, str], dict[str, object]],
    key: tuple[str, str | None, str],
    row: pd.Series,
) -> None:
    bucket = stats.setdefault(key, _new_evidence_bucket())
    bucket["ms_hit_count"] = int(bucket["ms_hit_count"]) + 1
    for source_col, output_key in (
        ("mhc_restriction", "ms_alleles"),
        ("pmid", "ms_pmids"),
    ):
        value = row.get(source_col)
        if isinstance(value, str):
            stripped = value.strip()
            if stripped:
                bucket[output_key].add(stripped)
        elif pd.notna(value):
            bucket[output_key].add(str(value))
    for source_col in ("cell_line_name", "cell_name", "source_tissue"):
        value = row.get(source_col)
        if isinstance(value, str):
            stripped = value.strip()
            if stripped:
                bucket["ms_samples"].add(stripped)


def _finalize_evidence_bucket(bucket: dict[str, object]) -> dict[str, object]:
    return {
        "ms_hit_count": int(bucket["ms_hit_count"]),
        "ms_alleles": ";".join(sorted(bucket["ms_alleles"])),
        "ms_pmids": ";".join(sorted(bucket["ms_pmids"])),
        "ms_samples": ";".join(sorted(bucket["ms_samples"])),
    }


def _is_truthy(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {"true", "1", "yes"}


def _is_best_sample_allele(
    peptide: str,
    allele: str,
    sample_alleles: set[str],
    score_lookup: dict[tuple[str, str], float],
) -> bool:
    allele_score = score_lookup.get((peptide, allele))
    if allele_score is None:
        return False
    sample_scores = [
        score_lookup[(peptide, sample_allele)]
        for sample_allele in sample_alleles
        if (peptide, sample_allele) in score_lookup
    ]
    if not sample_scores:
        return False
    return allele_score <= min(sample_scores) + 1e-12


def _build_evidence_stats(
    hits: pd.DataFrame,
    allele_list: list[str],
    score_lookup: dict[tuple[str, str], float],
) -> dict[tuple[str, str | None, str], dict[str, object]]:
    """Build evidence buckets keyed by peptide, allele, and evidence tier."""
    from .mhc import split_mhc_restrictions

    stats: dict[tuple[str, str | None, str], dict[str, object]] = {}
    if hits.empty or "peptide" not in hits.columns:
        return stats

    panel_alleles = set(allele_list)
    for _, row in hits.iterrows():
        peptide_value = row.get("peptide")
        if not isinstance(peptide_value, str) or not peptide_value.strip():
            continue
        peptide = peptide_value.strip()
        restrictions = set(split_mhc_restrictions(row.get("mhc_restriction", "")))
        exact_panel_alleles = {
            allele for allele in restrictions if allele in panel_alleles and "*" in allele
        }
        is_monoallelic = _is_truthy(row.get("is_monoallelic", False))
        matched = False

        if is_monoallelic and exact_panel_alleles:
            for allele in exact_panel_alleles:
                _add_evidence(stats, (peptide, allele, "monoallelic_ms"), row)
            matched = True
        elif exact_panel_alleles:
            for allele in exact_panel_alleles:
                _add_evidence(stats, (peptide, allele, "sample_allele_ms"), row)
            matched = True
        else:
            sample_alleles = set(split_mhc_restrictions(row.get("mhc_allele_set", "")))
            provenance = str(row.get("mhc_allele_provenance", "")).strip()
            if provenance == "sample_allele_match" and sample_alleles:
                for allele in sorted(panel_alleles & sample_alleles):
                    if _is_best_sample_allele(peptide, allele, sample_alleles, score_lookup):
                        _add_evidence(stats, (peptide, allele, "sample_allele_ms"), row)
                        matched = True

        has_specific_restriction = any(
            allele.startswith("HLA-") and "*" in allele for allele in restrictions
        )
        if not matched and not has_specific_restriction:
            _add_evidence(stats, (peptide, None, "unrestricted_ms"), row)

    return {
        key: _finalize_evidence_bucket(bucket)
        for key, bucket in stats.items()
        if int(bucket["ms_hit_count"]) > 0
    }


def _candidate_rows(
    scores: pd.DataFrame,
    pep_df: pd.DataFrame,
    allele_list: list[str],
    evidence_stats: dict[tuple[str, str | None, str], dict[str, object]],
    cutoffs: dict[str, float],
    include_predicted_only: bool,
) -> pd.DataFrame:
    columns = [*_LONG_OUTPUT_COLUMNS, "_tier_rank"]
    if scores.empty or "presentation_percentile" not in scores.columns:
        return pd.DataFrame(columns=columns)

    pep_meta = (
        pep_df[["peptide", "gene_name", "length"]]
        .drop_duplicates(subset=["peptide", "gene_name", "length"])
        .rename(columns={"gene_name": "cta"})
    )
    scored = scores[scores["allele"].isin(allele_list)].merge(pep_meta, on="peptide", how="inner")
    if scored.empty:
        return pd.DataFrame(columns=columns)

    rows = []
    for row in scored.itertuples(index=False):
        peptide = str(row.peptide)
        allele = str(row.allele)
        percentile = float(row.presentation_percentile)
        chosen_tier = None
        evidence = None
        for tier in ("monoallelic_ms", "sample_allele_ms", "unrestricted_ms"):
            key = (peptide, allele, tier)
            if tier == "unrestricted_ms" and key not in evidence_stats:
                key = (peptide, None, tier)
            tier_evidence = evidence_stats.get(key)
            if tier_evidence is not None and percentile < cutoffs[tier]:
                chosen_tier = tier
                evidence = tier_evidence
                break

        if chosen_tier is None:
            if not include_predicted_only or percentile >= cutoffs["predicted_only"]:
                continue
            chosen_tier = "predicted_only"
            evidence = {
                "ms_hit_count": 0,
                "ms_alleles": "",
                "ms_pmids": "",
                "ms_samples": "",
            }

        rows.append(
            {
                "cta": row.cta,
                "allele": allele,
                "peptide": peptide,
                "length": row.length,
                "evidence_tier": chosen_tier,
                "ms_hit_count": evidence["ms_hit_count"],
                "ms_alleles": evidence["ms_alleles"],
                "ms_pmids": evidence["ms_pmids"],
                "ms_samples": evidence["ms_samples"],
                "presentation_percentile": percentile,
                "presentation_score": getattr(row, "presentation_score", pd.NA),
                "affinity_nm": getattr(row, "affinity_nm", pd.NA),
                "_tier_rank": _EVIDENCE_TIER_RANKS[chosen_tier],
            }
        )

    return pd.DataFrame(rows, columns=columns)


def _best_per_cell(candidates: pd.DataFrame) -> pd.DataFrame:
    """For each (cta, allele), pick highest evidence tier, then best percentile."""
    if candidates.empty:
        return pd.DataFrame(columns=list(_LONG_OUTPUT_COLUMNS))

    sort_cols = ["_tier_rank", "presentation_percentile", "peptide", "allele"]
    best = (
        candidates.sort_values(sort_cols, kind="mergesort")
        .groupby(["cta", "allele"], as_index=False)
        .first()
    )
    return best[list(_LONG_OUTPUT_COLUMNS)]


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
