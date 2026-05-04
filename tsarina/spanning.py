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
        panel="global53_abc",
        lengths=(8, 9, 10, 11),
    )
"""

from __future__ import annotations

import sys
import time
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import TextIO

import pandas as pd

_DEFAULT_RANK_COLUMN = "ms_cancer_peptide_count"
_DEFAULT_PANEL = "global53_abc"
_DEFAULT_LENGTHS = (8, 9, 10, 11)
_DEFAULT_MONOALLELIC_MS_MAX_PERCENTILE = 2.0
_DEFAULT_SAMPLE_ALLELE_MS_MAX_PERCENTILE = 1.0
_DEFAULT_UNRESTRICTED_MS_MAX_PERCENTILE = 0.5
_DEFAULT_PREDICTED_ONLY_MAX_PERCENTILE = 0.1
_DEFAULT_PEPTIDES_PER_CELL = 3
_DEFAULT_SELECTION_ALLOWLIST = ("PRAME", "NY-ESO-1", "MAGEA4")
_DEFAULT_VITAL_TISSUE_MAX_NTPM = 2.0

_CTA_GROUPS: dict[str, tuple[str, ...]] = {
    "NY-ESO-1": ("CTAG1A", "CTAG1B"),
}

_VITAL_TISSUE_RNA_COLUMNS: tuple[str, ...] = (
    "rna_brain_max_ntpm",
    "rna_heart_max_ntpm",
    "rna_lung_max_ntpm",
    "rna_liver_max_ntpm",
    "rna_pancreas_max_ntpm",
)
_VITAL_TISSUE_MS_NAMES: frozenset[str] = frozenset(
    {
        "brain",
        "central nervous system (cns)",
        "cerebellum",
        "heart",
        "lung",
        "liver",
        "pancreas",
    }
)

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
    "ms_source_count",
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
    selection_allowlist: Iterable[str] | None = _DEFAULT_SELECTION_ALLOWLIST,
    exclude_vital_tissue_expression: bool = True,
    vital_tissue_max_ntpm: float = _DEFAULT_VITAL_TISSUE_MAX_NTPM,
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
    peptides_per_cell: int = _DEFAULT_PEPTIDES_PER_CELL,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    output_format: str = "wide",
    on_progress: Callable[[str], None] | None = None,
    progress_bar: bool = False,
    score_chunk_size: int | None = None,
    progress_file: TextIO | None = None,
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
        ranking. CTAs in ``selection_allowlist`` bypass this gate.
    restriction_levels
        Optional set of ``restriction`` values to keep (e.g.
        ``("TESTIS", "PLACENTAL")``).  ``None`` (default) keeps all.
    selection_allowlist
        CTA names allowed through automatic safety/confidence gates. Defaults
        to clinically anchored CTAs ``("PRAME", "NY-ESO-1", "MAGEA4")``.
        Aliases such as ``"MAGE-A4"`` and ``"CTAG1B"`` are normalized.
    exclude_vital_tissue_expression
        If True (default), automatic CTA selection excludes genes with RNA
        above ``vital_tissue_max_ntpm`` or unique public healthy-MS
        observations in brain/CNS/cerebellum, heart, lung, liver, or pancreas,
        unless the CTA is in ``selection_allowlist``.
    vital_tissue_max_ntpm
        Maximum allowed RNA nTPM in vital tissues for automatic CTA selection.
        Default 2.0.
    alleles
        Explicit allele list.  Overrides ``panel``.
    panel
        Named panel from :mod:`tsarina.alleles`.  Default
        ``"global53_abc"`` (53 globally broad HLA-A/B/C alleles
        constrained to MHCflurry affinity percentile-rank calibrated alleles
        and augmented with CTA-MS supported alleles).
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
    peptides_per_cell
        Maximum number of ranked peptides to keep for each CTA x HLA cell.
        Default 3. Candidates are ranked by MS source count first, then
        MS hit count, then prediction percentile.
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
        each pipeline stage.  Used by the CLI handler to emit
        stderr progress lines on multi-minute runs; library callers
        leave it ``None`` (default) for silent operation.  The library
        itself never prints — it only formats messages and hands them
        to the callback.
    progress_bar
        If True, score alleles in chunks and render a tqdm progress bar for
        the scoring stage. Default False so library calls remain quiet.
    score_chunk_size
        Number of HLA alleles per scoring chunk when ``progress_bar`` is
        enabled. Defaults to 8 for MHCflurry and one allele per chunk for
        other predictors.
    progress_file
        File handle for tqdm progress bars. Defaults to ``sys.stderr`` when
        ``progress_bar`` is True.

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
    if peptides_per_cell < 1:
        raise ValueError(f"peptides_per_cell must be >= 1, got {peptides_per_cell!r}")

    if max_percentile is not None:
        monoallelic_ms_max_percentile = max_percentile
        sample_allele_ms_max_percentile = max_percentile
        unrestricted_ms_max_percentile = max_percentile

    from .version import __version__

    _report_progress(on_progress, f"tsarina v{__version__}")

    allele_list = _resolve_alleles(alleles, panel)
    allele_frequencies = _weighted_allele_frequencies(allele_list)
    if alleles is None:
        allele_list = _order_alleles_by_frequency(allele_list, allele_frequencies)
    cta_list = _resolve_ctas(
        ctas=ctas,
        cta_count=cta_count,
        cta_rank_by=cta_rank_by,
        min_restriction_confidence=min_restriction_confidence,
        restriction_levels=restriction_levels,
        selection_allowlist=selection_allowlist,
        exclude_vital_tissue_expression=exclude_vital_tissue_expression,
        vital_tissue_max_ntpm=vital_tissue_max_ntpm,
    )
    cta_rank_values = _cta_rank_values(cta_list, cta_rank_by)

    _report_progress(
        on_progress,
        f"Panel inputs: {len(cta_list)} CTAs x {len(allele_list)} HLA alleles; "
        f"lengths {','.join(str(length) for length in lengths)}; "
        f"keeping up to {peptides_per_cell} peptides per CTA x HLA cell.",
    )

    if not cta_list or not allele_list:
        return _empty_output(
            cta_list, allele_list, output_format, allele_frequencies, cta_rank_values
        )

    _report_progress(
        on_progress,
        "Enumerating CTA-exclusive peptides..."
        if require_cta_exclusive
        else "Enumerating CTA peptides without the exclusivity gate...",
    )
    pep_df = _resolve_peptides(
        cta_list=cta_list,
        lengths=lengths,
        ensembl_release=ensembl_release,
        require_cta_exclusive=require_cta_exclusive,
        on_progress=on_progress,
        progress_bar=progress_bar,
        progress_file=progress_file,
    )
    if pep_df.empty:
        _report_progress(on_progress, "No CTA peptides survived peptide enumeration.")
        return _empty_output(
            cta_list, allele_list, output_format, allele_frequencies, cta_rank_values
        )

    from .ms_evidence import load_public_ms_hits

    unique_peptides = pep_df["peptide"].unique().tolist()
    _report_progress(
        on_progress,
        f"Enumerated {len(pep_df)} CTA peptide rows ({len(unique_peptides)} unique peptides).",
    )
    _report_progress(on_progress, "Loading public MS evidence for enumerated peptides...")
    hits = load_public_ms_hits(
        unique_peptides,
        iedb_path=iedb_path,
        cedar_path=cedar_path,
        mhc_class="I",
        mhc_species="Homo sapiens",
    )
    _report_progress(on_progress, f"Loaded {len(hits)} public MS evidence rows.")

    score_alleles = _score_alleles_for_panel(allele_list, hits)
    n_predictions = len(unique_peptides) * len(score_alleles)
    _report_progress(
        on_progress,
        f"Scoring {len(unique_peptides)} peptides x {len(score_alleles)} alleles "
        f"({n_predictions} predictions) via {predictor}...",
    )

    scoring_start = time.perf_counter()
    scores = _score_presentations(
        peptides=unique_peptides,
        alleles=score_alleles,
        predictor=predictor,
        progress_bar=progress_bar,
        score_chunk_size=score_chunk_size,
        progress_file=progress_file,
    )
    scoring_elapsed = time.perf_counter() - scoring_start

    _report_progress(on_progress, f"Scored {n_predictions} predictions in {scoring_elapsed:.1f}s.")

    evidence_stats = _build_evidence_stats(hits, allele_list, _score_lookup(scores))
    tier_counts = _evidence_tier_counts(evidence_stats)
    _report_progress(
        on_progress,
        "Built MS evidence tiers: "
        + ", ".join(f"{tier}={tier_counts.get(tier, 0)}" for tier in _EVIDENCE_TIER_RANKS),
    )
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
    selected = _top_per_cell(candidates, peptides_per_cell=peptides_per_cell)
    summary = panel_summary(
        selected=selected,
        cta_list=cta_list,
        allele_list=allele_list,
        allele_frequencies=allele_frequencies,
    )

    _report_progress(
        on_progress,
        f"Panel selected {summary['selected_peptide_count']} unique peptides across "
        f"{summary['filled_cell_count']}/{summary['possible_cell_count']} CTA x HLA cells "
        f"using MS tier cutoffs {monoallelic_ms_max_percentile}/"
        f"{sample_allele_ms_max_percentile}/{unrestricted_ms_max_percentile}.",
    )

    if output_format == "wide":
        out = _to_wide(selected, cta_list, allele_list)
    else:
        out = _to_long(selected, cta_list, allele_list)
    return _attach_panel_attrs(
        out, cta_list, allele_list, allele_frequencies, cta_rank_values, summary
    )


# ── Helpers ────────────────────────────────────────────────────────────


def _report_progress(on_progress: Callable[[str], None] | None, message: str) -> None:
    if on_progress is not None:
        on_progress(message)


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


def _weighted_allele_frequencies(allele_list: list[str]) -> dict[str, float]:
    """Return population-weighted allele frequencies for known panel alleles.

    The region table can contain multiple proxy populations per region. Use the
    maximum frequency per region and allele so an extra validation proxy does
    not overweight that region, then weight regions by population.
    """
    try:
        from .regions import REGION_POPULATIONS, region_allele_frequencies
    except ImportError:
        return dict.fromkeys(allele_list, 0.0)

    freqs = region_allele_frequencies()
    if freqs.empty or "frequency" not in freqs.columns:
        return dict.fromkeys(allele_list, 0.0)

    freqs = freqs[freqs["allele"].isin(allele_list)].copy()
    freqs = freqs[pd.notna(freqs["frequency"])]
    if freqs.empty:
        return dict.fromkeys(allele_list, 0.0)

    per_region = freqs.groupby(["region", "allele"], as_index=False)["frequency"].max()
    total_population = sum(float(population) for population in REGION_POPULATIONS.values())
    weighted: dict[str, float] = dict.fromkeys(allele_list, 0.0)
    for row in per_region.itertuples(index=False):
        population = REGION_POPULATIONS.get(row.region)
        if population is None:
            continue
        weighted[row.allele] += float(row.frequency) * float(population) / total_population
    return weighted


def _order_alleles_by_frequency(
    allele_list: list[str],
    allele_frequencies: dict[str, float],
) -> list[str]:
    return sorted(allele_list, key=lambda allele: (-allele_frequencies.get(allele, 0.0), allele))


def _resolve_ctas(
    ctas: Iterable[str] | None,
    cta_count: int,
    cta_rank_by: str,
    min_restriction_confidence: Iterable[str] | None,
    restriction_levels: Iterable[str] | None,
    selection_allowlist: Iterable[str] | None = _DEFAULT_SELECTION_ALLOWLIST,
    exclude_vital_tissue_expression: bool = True,
    vital_tissue_max_ntpm: float = _DEFAULT_VITAL_TISSUE_MAX_NTPM,
) -> list[str]:
    """Pick the CTA gene set per the resolution order: explicit override
    wins; else filter the bundled CTA CSV and take top-N by rank."""
    from .gene_sets import CTA_gene_names
    from .loader import cta_dataframe

    valid = CTA_gene_names()
    allowlist = _normalize_cta_labels(selection_allowlist, valid)

    if ctas is not None:
        return _normalize_cta_labels(ctas, valid)

    df = cta_dataframe()
    if "Symbol" not in df.columns:
        return []
    df = df.copy()
    df["_cta_display"] = df["Symbol"].map(_cta_display_name)
    df["_is_selection_allowlisted"] = df["_cta_display"].isin(allowlist)
    if "filtered" in df.columns:
        df = df[df["filtered"].astype(str).str.lower() == "true"]
    if "never_expressed" in df.columns:
        df = df[~(df["never_expressed"].astype(str).str.lower() == "true")]

    if min_restriction_confidence is not None and "restriction_confidence" in df.columns:
        wanted = {c.upper() for c in min_restriction_confidence}
        confidence_ok = df["restriction_confidence"].astype(str).str.upper().isin(wanted)
        df = df[confidence_ok | df["_is_selection_allowlisted"]]

    if restriction_levels is not None and "restriction" in df.columns:
        wanted_levels = set(restriction_levels)
        df = df[df["restriction"].isin(wanted_levels)]

    if exclude_vital_tissue_expression:
        vital_rna = df.apply(
            lambda row: _has_vital_tissue_rna_expression(row, vital_tissue_max_ntpm),
            axis=1,
        )
        vital_ms_labels = _unique_vital_healthy_ms_cta_labels(df[~df["_is_selection_allowlisted"]])
        vital_ms = df["_cta_display"].isin(vital_ms_labels)
        df = df[~(vital_rna | vital_ms) | df["_is_selection_allowlisted"]]

    if cta_rank_by in df.columns and df[cta_rank_by].notna().any():
        df = df.sort_values(cta_rank_by, ascending=False, na_position="last")
    else:
        df = df.sort_values("Symbol")

    df = df[df["Symbol"].isin(valid)]
    df = df.drop_duplicates(subset=["_cta_display"], keep="first")
    return df["_cta_display"].head(cta_count).tolist()


def _compact_cta_name(value: object) -> str:
    return "".join(ch for ch in str(value).upper() if ch.isalnum())


def _cta_display_name(symbol: object) -> str:
    token = str(symbol).strip()
    compact = _compact_cta_name(token)
    if compact in {"NYESO1", "NYESO", "CTAG1", "CTAG1A", "CTAG1B"}:
        return "NY-ESO-1"
    if compact.startswith("MAGEA") and compact[5:].isdigit():
        return compact
    return token


def _normalize_cta_labels(values: Iterable[str] | None, valid_symbols: set[str]) -> list[str]:
    if values is None:
        return []

    valid_by_compact = {_compact_cta_name(symbol): symbol for symbol in valid_symbols}
    labels: list[str] = []
    seen: set[str] = set()
    for value in values:
        compact = _compact_cta_name(value)
        if not compact:
            continue
        if compact in {"NYESO1", "NYESO", "CTAG1", "CTAG1A", "CTAG1B"}:
            label = "NY-ESO-1"
        elif compact in valid_by_compact:
            label = _cta_display_name(valid_by_compact[compact])
        else:
            continue
        if label in seen:
            continue
        if not any(symbol in valid_symbols for symbol in _cta_member_symbols(label)):
            continue
        labels.append(label)
        seen.add(label)
    return labels


def _cta_member_symbols(label: str) -> tuple[str, ...]:
    return _CTA_GROUPS.get(label, (label,))


def _expand_cta_symbols(labels: Iterable[str]) -> list[str]:
    expanded: list[str] = []
    seen: set[str] = set()
    for label in labels:
        for symbol in _cta_member_symbols(label):
            if symbol not in seen:
                expanded.append(symbol)
                seen.add(symbol)
    return expanded


def _split_semicolon_values(value: object) -> set[str]:
    if not isinstance(value, str):
        return set()
    return {part.strip() for part in value.split(";") if part.strip()}


def _has_vital_tissue_rna_expression(row: pd.Series, max_ntpm: float) -> bool:
    for column in _VITAL_TISSUE_RNA_COLUMNS:
        value = row.get(column)
        try:
            ntpm = float(value)
        except (TypeError, ValueError):
            continue
        if pd.notna(ntpm) and ntpm > max_ntpm:
            return True
    return False


def _has_vital_healthy_ms_tissue(row: pd.Series) -> bool:
    healthy_ms_tissues = {
        tissue.lower() for tissue in _split_semicolon_values(row.get("ms_healthy_somatic_tissues"))
    }
    return bool(healthy_ms_tissues & _VITAL_TISSUE_MS_NAMES)


def _unique_vital_healthy_ms_cta_labels(candidate_df: pd.DataFrame) -> set[str]:
    """CTA labels with unique-gene public healthy-MS evidence in vital tissue.

    The bundled CTA table carries gene-level healthy-MS tissue aggregates, but
    those aggregates can be driven by peptides shared across paralogous CTA
    families. Use that table only as a cheap suspect list, then verify against
    hitlist observation rows and require ``gene_names`` to map uniquely to the
    candidate CTA symbol before vetoing the whole CTA.
    """
    if candidate_df.empty or "Symbol" not in candidate_df.columns:
        return set()

    suspect_symbols: set[str] = set()
    for _, row in candidate_df.iterrows():
        if not _has_vital_healthy_ms_tissue(row):
            continue
        for symbol in _cta_member_symbols(str(row.get("_cta_display", row.get("Symbol")))):
            suspect_symbols.add(symbol)
    if not suspect_symbols:
        return set()

    from .indexing import load_ms_evidence

    hits = load_ms_evidence(
        gene_name=sorted(suspect_symbols),
        mhc_class="I",
        mhc_species="Homo sapiens",
        columns=[
            "peptide",
            "gene_names",
            "source_tissue",
            "src_healthy_tissue",
            "is_binding_assay",
        ],
        drop_binding_assays=True,
    )
    if hits.empty:
        return set()

    veto_labels: set[str] = set()
    for _, row in hits.iterrows():
        if not _is_truthy(row.get("src_healthy_tissue", False)):
            continue
        source_tissue = str(row.get("source_tissue", "")).strip().lower()
        if source_tissue not in _VITAL_TISSUE_MS_NAMES:
            continue
        row_symbols = _split_semicolon_values(row.get("gene_names"))
        if len(row_symbols) != 1:
            continue
        symbol = next(iter(row_symbols))
        if symbol in suspect_symbols:
            veto_labels.add(_cta_display_name(symbol))
    return veto_labels


def _cta_rank_values(cta_list: list[str], cta_rank_by: str) -> dict[str, float | str]:
    from .loader import cta_dataframe

    if not cta_list:
        return {}

    df = cta_dataframe()
    if "Symbol" not in df.columns or cta_rank_by not in df.columns:
        return {}

    rank_values: dict[str, float | str] = {}
    for cta in cta_list:
        values = df[df["Symbol"].isin(_cta_member_symbols(cta))][cta_rank_by].dropna()
        if values.empty:
            continue
        numeric = pd.to_numeric(values, errors="coerce")
        if numeric.notna().any():
            rank_values[cta] = float(numeric.max())
        else:
            rank_values[cta] = str(values.iloc[0])
    return rank_values


def _resolve_peptides(
    cta_list: list[str],
    lengths: tuple[int, ...],
    ensembl_release: int,
    require_cta_exclusive: bool,
    on_progress: Callable[[str], None] | None = None,
    progress_bar: bool = False,
    progress_file: TextIO | None = None,
) -> pd.DataFrame:
    cta_symbols = _expand_cta_symbols(cta_list)
    if require_cta_exclusive:
        from .peptides import cta_exclusive_peptides

        all_peps = cta_exclusive_peptides(
            ensembl_release=ensembl_release,
            lengths=tuple(lengths),
            gene_names=cta_symbols,
            on_progress=on_progress,
            progress_bar=progress_bar,
            progress_file=progress_file,
        )
    else:
        from .peptides import cta_peptides

        all_peps = cta_peptides(
            ensembl_release=ensembl_release,
            lengths=tuple(lengths),
            gene_names=cta_symbols,
            on_progress=on_progress,
            progress_bar=progress_bar,
            progress_file=progress_file,
        )
    if all_peps.empty:
        return all_peps
    out = all_peps[all_peps["gene_name"].isin(cta_symbols)].copy()
    out["gene_name"] = out["gene_name"].map(_cta_display_name)
    return out[out["gene_name"].isin(cta_list)].copy()


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


def _score_presentations(
    peptides: list[str],
    alleles: list[str],
    predictor: str,
    progress_bar: bool,
    score_chunk_size: int | None,
    progress_file: TextIO | None,
) -> pd.DataFrame:
    from .scoring import score_presentation

    if not progress_bar:
        return score_presentation(peptides=peptides, alleles=alleles, predictor=predictor)

    chunk_size = score_chunk_size or (
        8 if predictor.lower().replace("-", "_") == "mhcflurry" else 1
    )
    if chunk_size < 1:
        raise ValueError(f"score_chunk_size must be >= 1, got {score_chunk_size!r}")

    chunks = [alleles[i : i + chunk_size] for i in range(0, len(alleles), chunk_size)]
    frames: list[pd.DataFrame] = []
    from tqdm import tqdm

    for allele_chunk in tqdm(
        chunks,
        desc="Scoring pMHCs",
        unit="allele-chunk",
        file=progress_file or sys.stderr,
        leave=False,
    ):
        frames.append(
            score_presentation(peptides=peptides, alleles=allele_chunk, predictor=predictor)
        )
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


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
    ms_source_count = len(bucket["ms_pmids"])
    if ms_source_count == 0:
        ms_source_count = len(bucket["ms_samples"])
    if ms_source_count == 0 and int(bucket["ms_hit_count"]) > 0:
        ms_source_count = 1
    return {
        "ms_hit_count": int(bucket["ms_hit_count"]),
        "ms_source_count": ms_source_count,
        "ms_alleles": ";".join(sorted(bucket["ms_alleles"])),
        "ms_pmids": ";".join(sorted(bucket["ms_pmids"])),
        "ms_samples": ";".join(sorted(bucket["ms_samples"])),
    }


def _evidence_tier_counts(
    evidence_stats: dict[tuple[str, str | None, str], dict[str, object]],
) -> dict[str, int]:
    counts = dict.fromkeys(_EVIDENCE_TIER_RANKS, 0)
    for _, _, tier in evidence_stats:
        counts[tier] = counts.get(tier, 0) + 1
    return counts


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
                "ms_source_count": 0,
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
                "ms_source_count": evidence["ms_source_count"],
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


def _top_per_cell(candidates: pd.DataFrame, peptides_per_cell: int) -> pd.DataFrame:
    """For each (cta, allele), keep the top-N peptides by evidence support."""
    if candidates.empty:
        return pd.DataFrame(columns=list(_LONG_OUTPUT_COLUMNS))

    ranked = candidates.copy()
    cell_has_ms = ranked.groupby(["cta", "allele"])["ms_hit_count"].transform("max") > 0
    ranked = ranked[(ranked["ms_hit_count"].astype(int) > 0) | ~cell_has_ms].copy()
    ranked["_source_rank"] = -ranked["ms_source_count"].astype(int)
    ranked["_hit_rank"] = -ranked["ms_hit_count"].astype(int)
    sort_cols = [
        "cta",
        "allele",
        "_source_rank",
        "_hit_rank",
        "presentation_percentile",
        "_tier_rank",
        "peptide",
    ]
    ranked = ranked.sort_values(sort_cols, kind="mergesort")
    ranked["peptide_rank_in_cell"] = ranked.groupby(["cta", "allele"]).cumcount() + 1
    selected = ranked[ranked["peptide_rank_in_cell"] <= peptides_per_cell].copy()
    columns = [*_LONG_OUTPUT_COLUMNS, "peptide_rank_in_cell"]
    return selected[columns].reset_index(drop=True)


def _to_wide(
    selected: pd.DataFrame,
    cta_list: list[str],
    allele_list: list[str],
) -> pd.DataFrame:
    if selected.empty:
        wide = pd.DataFrame(index=cta_list, columns=allele_list, dtype=object)
    else:
        cell_values = (
            selected.sort_values(["cta", "allele", "peptide_rank_in_cell"])
            .groupby(["cta", "allele"])["peptide"]
            .apply("; ".join)
            .reset_index()
        )
        wide = cell_values.pivot(index="cta", columns="allele", values="peptide")
        wide = wide.reindex(index=cta_list, columns=allele_list)
    wide.index.name = "cta"
    wide = wide.reset_index()
    return wide


def _to_long(
    selected: pd.DataFrame,
    cta_list: list[str],
    allele_list: list[str],
) -> pd.DataFrame:
    columns = [*_LONG_OUTPUT_COLUMNS, "peptide_rank_in_cell"]
    if selected.empty:
        return pd.DataFrame(columns=columns)
    cta_order = {gene: i for i, gene in enumerate(cta_list)}
    allele_order = {a: i for i, a in enumerate(allele_list)}
    out = selected.copy()
    out["_cta_order"] = out["cta"].map(cta_order)
    out["_allele_order"] = out["allele"].map(allele_order)
    out = out.sort_values(["_cta_order", "_allele_order", "peptide_rank_in_cell"]).drop(
        columns=["_cta_order", "_allele_order"]
    )
    return out[columns].reset_index(drop=True)


def _carrier_probability(allele_frequency: float) -> float:
    allele_frequency = max(0.0, min(float(allele_frequency), 1.0))
    return 1.0 - (1.0 - allele_frequency) ** 2


def _combined_carrier_probability(allele_frequencies: Iterable[float]) -> float:
    miss_probability = 1.0
    for frequency in allele_frequencies:
        miss_probability *= 1.0 - _carrier_probability(float(frequency))
    return max(0.0, min(1.0, 1.0 - miss_probability))


def panel_summary(
    selected: pd.DataFrame,
    cta_list: list[str],
    allele_list: list[str],
    allele_frequencies: dict[str, float] | None = None,
) -> dict[str, object]:
    """Summarize a selected panel table.

    Population coverage estimates use weighted regional allele frequencies and
    a simple independent-carrier approximation. The values are intended for
    ranking and sanity-checking panel breadth, not for clinical-grade
    population-genetics inference.
    """
    allele_frequencies = allele_frequencies or dict.fromkeys(allele_list, 0.0)
    possible_cell_count = len(cta_list) * len(allele_list)

    if selected.empty:
        filled_pairs = pd.DataFrame(columns=["cta", "allele"])
        selected_peptide_count = 0
    else:
        filled_pairs = selected[["cta", "allele"]].drop_duplicates()
        selected_peptide_count = int(selected["peptide"].nunique())

    filled_cell_count = len(filled_pairs)
    avg_peptides = 0.0 if filled_cell_count == 0 else float(len(selected) / filled_cell_count)

    cta_rows = []
    for cta in cta_list:
        covered_alleles = (
            filled_pairs.loc[filled_pairs["cta"] == cta, "allele"].drop_duplicates().tolist()
            if not filled_pairs.empty
            else []
        )
        cta_selected = selected[selected["cta"] == cta] if not selected.empty else selected
        coverage = _combined_carrier_probability(
            allele_frequencies.get(allele, 0.0) for allele in covered_alleles
        )
        cta_rows.append(
            {
                "cta": cta,
                "covered_hla_count": len(covered_alleles),
                "selected_peptide_count": int(cta_selected["peptide"].nunique())
                if not cta_selected.empty
                else 0,
                "estimated_population_coverage": coverage,
            }
        )

    hla_rows = []
    total_ctas = len(cta_list)
    for allele in allele_list:
        covered_ctas = (
            filled_pairs.loc[filled_pairs["allele"] == allele, "cta"].drop_duplicates().tolist()
            if not filled_pairs.empty
            else []
        )
        allele_selected = selected[selected["allele"] == allele] if not selected.empty else selected
        hla_rows.append(
            {
                "allele": allele,
                "weighted_allele_frequency": allele_frequencies.get(allele, 0.0),
                "covered_cta_count": len(covered_ctas),
                "covered_cta_fraction": 0.0
                if total_ctas == 0
                else float(len(covered_ctas) / total_ctas),
                "selected_peptide_count": int(allele_selected["peptide"].nunique())
                if not allele_selected.empty
                else 0,
            }
        )

    tier_counts: dict[str, int] = {}
    if not selected.empty and "evidence_tier" in selected.columns:
        tier_counts = {
            str(tier): int(count)
            for tier, count in selected["evidence_tier"].value_counts(sort=False).items()
        }

    return {
        "hla_allele_count": len(allele_list),
        "cta_count": len(cta_list),
        "selected_peptide_count": selected_peptide_count,
        "selected_row_count": len(selected),
        "filled_cell_count": filled_cell_count,
        "possible_cell_count": possible_cell_count,
        "filled_cell_fraction": 0.0
        if possible_cell_count == 0
        else float(filled_cell_count / possible_cell_count),
        "average_peptides_per_filled_cell": avg_peptides,
        "evidence_tier_counts": tier_counts,
        "cta_coverage": cta_rows,
        "hla_coverage": hla_rows,
        "coverage_note": (
            "Estimated from weighted regional HLA allele frequencies using an "
            "independent-carrier approximation."
        ),
    }


def _attach_panel_attrs(
    out: pd.DataFrame,
    cta_list: list[str],
    allele_list: list[str],
    allele_frequencies: dict[str, float],
    cta_rank_values: dict[str, float | str],
    summary: dict[str, object],
) -> pd.DataFrame:
    out.attrs["cta_order"] = list(cta_list)
    out.attrs["allele_order"] = list(allele_list)
    out.attrs["allele_frequencies"] = dict(allele_frequencies)
    out.attrs["cta_rank_values"] = dict(cta_rank_values)
    out.attrs["panel_summary"] = summary
    return out


def _empty_output(
    cta_list: list[str],
    allele_list: list[str],
    output_format: str,
    allele_frequencies: dict[str, float] | None = None,
    cta_rank_values: dict[str, float | str] | None = None,
) -> pd.DataFrame:
    if output_format == "wide":
        wide = pd.DataFrame(index=cta_list, columns=allele_list, dtype=object)
        wide.index.name = "cta"
        out = wide.reset_index()
    else:
        out = pd.DataFrame(columns=[*_LONG_OUTPUT_COLUMNS, "peptide_rank_in_cell"])
    allele_frequencies = allele_frequencies or dict.fromkeys(allele_list, 0.0)
    summary = panel_summary(
        selected=pd.DataFrame(columns=[*_LONG_OUTPUT_COLUMNS, "peptide_rank_in_cell"]),
        cta_list=cta_list,
        allele_list=allele_list,
        allele_frequencies=allele_frequencies,
    )
    return _attach_panel_attrs(
        out,
        cta_list,
        allele_list,
        allele_frequencies,
        cta_rank_values or {},
        summary,
    )
