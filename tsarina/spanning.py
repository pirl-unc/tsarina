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
from collections import defaultdict
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import TextIO

import pandas as pd

_DEFAULT_RANK_COLUMN = "ms_cta_exclusive_cancer_peptide_count"
_DEFAULT_PANEL = "global53_abc"
_DEFAULT_LENGTHS = (8, 9, 10, 11)
_DEFAULT_MONOALLELIC_MS_MAX_PERCENTILE = 2.0
_DEFAULT_SAMPLE_ALLELE_MS_MAX_PERCENTILE = 1.0
_DEFAULT_UNRESTRICTED_MS_MAX_PERCENTILE = 0.5
_DEFAULT_PREDICTED_ONLY_MAX_PERCENTILE = 0.1
_DEFAULT_PEPTIDES_PER_CELL = 3
_CTAG1_GROUP_LABEL = "CTAG1A/CTAG1B"
_CTAG1_ALIASES = {
    "NYESO1",
    "NYESO",
    "CTAG1",
    "CTAG1A",
    "CTAG1B",
    "CTAG1AB",
    "CTAG1ACTAG1B",
    "CTAG1BCTAG1A",
}
_DEFAULT_SELECTION_ALLOWLIST = ("PRAME", _CTAG1_GROUP_LABEL, "MAGEA4")
_DEFAULT_VITAL_TISSUE_MAX_NTPM = 2.0
_DEFAULT_EXCLUDE_NON_MAGEA4_MAGE_FAMILY = True
_DEFAULT_GROUP_IDENTICAL_CTA_PMHCS = True
_DEFAULT_GROUP_IDENTICAL_CTA_PEPTIDE_SETS = True
_DEFAULT_ANNOTATE_NETMHCPAN_AFFINITY = False

_CTA_GROUPS: dict[str, tuple[str, ...]] = {
    _CTAG1_GROUP_LABEL: ("CTAG1A", "CTAG1B"),
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
    "cta_members",
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
    "affinity_percentile",
    "netmhcpan_affinity_nm",
    "netmhcpan_affinity_percentile",
)

_PMHC_SIGNATURE_COLUMNS: tuple[str, ...] = (
    "allele",
    "peptide",
    "length",
    "evidence_tier",
    "ms_hit_count",
    "ms_source_count",
    "ms_alleles",
    "ms_pmids",
    "ms_samples",
    "peptide_rank_in_cell",
)

_INTERNAL_SELECTED_COLUMNS: tuple[str, ...] = (
    "cta_peptide_set_label",
    "cta_peptide_set_members",
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
    exclude_non_magea4_mage_family: bool = _DEFAULT_EXCLUDE_NON_MAGEA4_MAGE_FAMILY,
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
    group_identical_cta_pmhcs: bool = _DEFAULT_GROUP_IDENTICAL_CTA_PMHCS,
    group_identical_cta_peptide_sets: bool = _DEFAULT_GROUP_IDENTICAL_CTA_PEPTIDE_SETS,
    annotate_netmhcpan_affinity: bool = _DEFAULT_ANNOTATE_NETMHCPAN_AFFINITY,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    output_format: str = "wide",
    include_empty_ctas: bool | None = None,
    on_progress: Callable[[str], None] | None = None,
    progress_bar: bool = False,
    score_chunk_size: int | None = None,
    progress_file: TextIO | None = None,
) -> pd.DataFrame:
    """Build a CTA x HLA pMHC panel matrix.

    Parameters
    ----------
    cta_count
        Number of downstream non-empty CTAs to include when ``ctas`` is not
        supplied. Automatic selection scans lower-ranked candidates as needed
        to backfill CTAs that produce no selected pMHCs. Passing
        ``include_empty_ctas=True`` restores the audit view of the top
        ``cta_count`` ranked candidates, including failures.
    cta_rank_by
        Column in the bundled CTA CSV used to rank candidates.  Default
        ``"ms_cta_exclusive_cancer_peptide_count"`` (CTA-exclusive
        peptides observed in cancer MS first). Falls back to alphabetical
        ``Symbol`` order when the column is missing or all-NaN.
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
        to clinically anchored CTAs ``("PRAME", "CTAG1A/CTAG1B",
        "MAGEA4")``. Aliases such as ``"NY-ESO-1"``, ``"MAGE-A4"``,
        and ``"CTAG1B"`` are normalized.
    exclude_vital_tissue_expression
        If True (default), automatic CTA selection excludes genes with RNA
        above ``vital_tissue_max_ntpm`` or unique public healthy-MS
        observations in brain/CNS/cerebellum, heart, lung, liver, or pancreas,
        unless the CTA is in ``selection_allowlist``.
    vital_tissue_max_ntpm
        Maximum allowed RNA nTPM in vital tissues for automatic CTA selection.
        Default 2.0.
    exclude_non_magea4_mage_family
        If True (default), automatic CTA selection excludes MAGE-family targets
        other than ``MAGEA4`` unless the target is in ``selection_allowlist``.
        Explicit ``ctas`` always win and bypass this gate.
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
    group_identical_cta_pmhcs
        If True (default), non-empty CTAs with identical final selected pMHC
        rows are collapsed into one output target. The long output preserves
        semicolon-separated source labels in ``cta_members``.
    group_identical_cta_peptide_sets
        If True (default), non-empty CTAs whose enumerated peptide sets are
        identical are collapsed before final pMHC-signature grouping. This
        combines paralog targets such as XAGE/PAGE family members when they
        genuinely present the same CTA-exclusive peptide set.
    annotate_netmhcpan_affinity
        If True, run a NetMHCpan binding-affinity pass on selected pMHC rows
        and add ``netmhcpan_affinity_nm`` plus
        ``netmhcpan_affinity_percentile``. This is opt-in because NetMHCpan is
        an external dependency and slower than the default MHCflurry path.
    iedb_path, cedar_path
        Optional explicit raw public-MS files. Defaults use hitlist's cached
        observations path.
    output_format
        ``"wide"`` (default) — pivot with CTA rows, allele columns, and
        the chosen peptide as the cell value.
        ``"long"`` — one row per filled cell with peptide, length,
        percentile, score, and affinity columns.
    include_empty_ctas
        Whether to keep CTA rows with no selected pMHCs after peptide,
        exclusivity, public-MS, and prediction gates. ``None`` (default)
        keeps explicit ``ctas`` requests verbatim but hides empty rows for
        automatic top-N selection.
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
    if cta_count < 1:
        raise ValueError(f"cta_count must be >= 1, got {cta_count!r}")
    if peptides_per_cell < 1:
        raise ValueError(f"peptides_per_cell must be >= 1, got {peptides_per_cell!r}")
    keep_empty_ctas = ctas is not None if include_empty_ctas is None else include_empty_ctas

    if max_percentile is not None:
        monoallelic_ms_max_percentile = max_percentile
        sample_allele_ms_max_percentile = max_percentile
        unrestricted_ms_max_percentile = max_percentile

    from .version import __version__

    _report_progress(on_progress, f"tsarina v{__version__}")

    allele_list = _resolve_alleles(alleles, panel)
    frequency_audit = _frequency_audit_for_alleles(allele_list)
    allele_frequencies = _frequencies_from_audit(frequency_audit, allele_list)
    allele_frequency_sources = _frequency_sources_from_audit(frequency_audit, allele_list)
    if alleles is None:
        allele_list = _order_alleles_by_frequency(allele_list, allele_frequencies)
    automatic_backfill = ctas is None and not keep_empty_ctas
    cta_candidates = _resolve_ctas(
        ctas=ctas,
        cta_count=None if automatic_backfill else cta_count,
        cta_rank_by=cta_rank_by,
        min_restriction_confidence=min_restriction_confidence,
        restriction_levels=restriction_levels,
        selection_allowlist=selection_allowlist,
        exclude_vital_tissue_expression=exclude_vital_tissue_expression,
        vital_tissue_max_ntpm=vital_tissue_max_ntpm,
        exclude_non_magea4_mage_family=exclude_non_magea4_mage_family,
    )
    cta_rank_values = _cta_rank_values(cta_candidates, cta_rank_by)

    cta_input_text = (
        f"up to {cta_count} non-empty CTAs from {len(cta_candidates)} ranked candidates"
        if automatic_backfill
        else f"{len(cta_candidates)} CTAs"
    )
    _report_progress(
        on_progress,
        f"Panel inputs: {cta_input_text} x {len(allele_list)} HLA alleles; "
        f"lengths {','.join(str(length) for length in lengths)}; "
        f"keeping up to {peptides_per_cell} peptides per CTA x HLA cell.",
    )

    if not cta_candidates or not allele_list:
        output_cta_list, empty_ctas = _output_cta_list(
            cta_candidates,
            _empty_selected_frame(),
            include_empty_ctas=keep_empty_ctas,
        )
        return _empty_output(
            output_cta_list,
            allele_list,
            output_format,
            allele_frequencies,
            cta_rank_values,
            allele_frequency_sources,
            frequency_audit,
            input_cta_list=cta_candidates,
            empty_ctas=empty_ctas,
            include_empty_ctas=keep_empty_ctas,
        )

    cutoffs = {
        "monoallelic_ms": monoallelic_ms_max_percentile,
        "sample_allele_ms": sample_allele_ms_max_percentile,
        "unrestricted_ms": unrestricted_ms_max_percentile,
        "predicted_only": predicted_only_max_percentile,
    }
    selected_frames: list[pd.DataFrame] = []
    processed_cta_list: list[str] = []
    batch_size = cta_count if automatic_backfill else len(cta_candidates)
    batch_size = max(1, batch_size)

    for start in range(0, len(cta_candidates), batch_size):
        batch_ctas = cta_candidates[start : start + batch_size]
        processed_cta_list.extend(batch_ctas)
        if automatic_backfill:
            _report_progress(
                on_progress,
                f"Evaluating CTA candidates {start + 1}-{start + len(batch_ctas)} "
                f"of {len(cta_candidates)} for downstream backfill...",
            )

        batch_selected = _select_cta_batch(
            cta_list=batch_ctas,
            allele_list=allele_list,
            lengths=lengths,
            ensembl_release=ensembl_release,
            require_cta_exclusive=require_cta_exclusive,
            iedb_path=iedb_path,
            cedar_path=cedar_path,
            predictor=predictor,
            cutoffs=cutoffs,
            include_predicted_only=include_predicted_only,
            peptides_per_cell=peptides_per_cell,
            group_identical_cta_peptide_sets=group_identical_cta_peptide_sets,
            on_progress=on_progress,
            progress_bar=progress_bar,
            score_chunk_size=score_chunk_size,
            progress_file=progress_file,
        )
        if not batch_selected.empty:
            selected_frames.append(batch_selected)

        selected_so_far = (
            pd.concat(selected_frames, ignore_index=True)
            if selected_frames
            else _empty_selected_frame()
        )
        nonempty_ctas = _nonempty_ctas_in_order(processed_cta_list, selected_so_far)
        if automatic_backfill:
            nonempty_cta_groups, _ = _cta_group_order(
                nonempty_ctas,
                selected_so_far,
                group_identical_cta_pmhcs=group_identical_cta_pmhcs,
                group_identical_cta_peptide_sets=group_identical_cta_peptide_sets,
            )
            _report_progress(
                on_progress,
                f"Backfill status: {len(nonempty_cta_groups)}/{cta_count} non-empty CTA targets.",
            )
            if len(nonempty_cta_groups) >= cta_count:
                break

    selected = (
        pd.concat(selected_frames, ignore_index=True)
        if selected_frames
        else _empty_selected_frame()
    )
    if automatic_backfill:
        nonempty_ctas = _nonempty_ctas_in_order(processed_cta_list, selected)
        grouped_cta_list, cta_groups = _cta_group_order(
            nonempty_ctas,
            selected,
            group_identical_cta_pmhcs=group_identical_cta_pmhcs,
            group_identical_cta_peptide_sets=group_identical_cta_peptide_sets,
        )
        output_cta_list = grouped_cta_list[:cta_count]
        output_cta_set = _cta_group_source_set(cta_groups[:cta_count])
        nonempty_cta_set = set(nonempty_ctas)
        empty_ctas = [cta for cta in processed_cta_list if cta not in nonempty_cta_set]
        selected = selected[selected["cta"].isin(output_cta_set)].reset_index(drop=True)
        cta_groups = cta_groups[:cta_count]
    else:
        ungrouped_output_cta_list, empty_ctas = _output_cta_list(
            processed_cta_list,
            selected,
            include_empty_ctas=keep_empty_ctas,
        )
        output_cta_list, cta_groups = _cta_group_order(
            ungrouped_output_cta_list,
            selected,
            group_identical_cta_pmhcs=group_identical_cta_pmhcs,
            group_identical_cta_peptide_sets=group_identical_cta_peptide_sets,
        )
        selected = selected[selected["cta"].isin(_cta_group_source_set(cta_groups))].reset_index(
            drop=True
        )
    selected = _apply_cta_groups(selected, cta_groups)
    cta_rank_values = _cta_rank_values_with_groups(cta_rank_values, cta_groups)
    if annotate_netmhcpan_affinity:
        selected = _annotate_netmhcpan_affinity(
            selected=selected,
            predictor=predictor,
            on_progress=on_progress,
            progress_bar=progress_bar,
            score_chunk_size=score_chunk_size,
            progress_file=progress_file,
        )
    summary = panel_summary(
        selected=selected,
        cta_list=output_cta_list,
        allele_list=allele_list,
        allele_frequencies=allele_frequencies,
        allele_frequency_sources=allele_frequency_sources,
    )
    summary["candidate_cta_count"] = len(cta_candidates)
    summary["input_cta_count"] = len(processed_cta_list)
    summary["empty_cta_count"] = len(empty_ctas)
    summary["empty_ctas"] = empty_ctas
    summary["include_empty_ctas"] = keep_empty_ctas
    summary["group_identical_cta_pmhcs"] = group_identical_cta_pmhcs
    summary["group_identical_cta_peptide_sets"] = group_identical_cta_peptide_sets
    summary["cta_groups"] = _public_cta_groups(cta_groups)
    summary["grouped_cta_member_count"] = sum(
        max(0, len(group["members"]) - 1) for group in cta_groups
    )
    summary["annotate_netmhcpan_affinity"] = annotate_netmhcpan_affinity

    _report_progress(
        on_progress,
        f"Panel selected {summary['selected_peptide_count']} unique peptides across "
        f"{summary['filled_cell_count']}/{summary['possible_cell_count']} CTA x HLA cells "
        f"using MS tier cutoffs {monoallelic_ms_max_percentile}/"
        f"{sample_allele_ms_max_percentile}/{unrestricted_ms_max_percentile}.",
    )
    if empty_ctas and not keep_empty_ctas:
        _report_progress(
            on_progress,
            f"Omitted {len(empty_ctas)} scanned CTA candidates with no selected pMHCs.",
        )
    grouped_member_count = int(summary.get("grouped_cta_member_count", 0))
    if grouped_member_count:
        _report_progress(
            on_progress,
            f"Grouped {grouped_member_count} redundant CTA member(s) with identical selected pMHCs.",
        )

    if output_format == "wide":
        out = _to_wide(selected, output_cta_list, allele_list)
    else:
        out = _to_long(selected, output_cta_list, allele_list)
    return _attach_panel_attrs(
        out,
        output_cta_list,
        allele_list,
        allele_frequencies,
        cta_rank_values,
        summary,
        allele_frequency_sources,
        frequency_audit,
        input_cta_list=processed_cta_list,
        empty_ctas=empty_ctas,
        cta_groups=cta_groups,
    )


# ── Helpers ────────────────────────────────────────────────────────────


def _report_progress(on_progress: Callable[[str], None] | None, message: str) -> None:
    if on_progress is not None:
        on_progress(message)


def _selected_columns() -> list[str]:
    return [*_LONG_OUTPUT_COLUMNS, *_INTERNAL_SELECTED_COLUMNS, "peptide_rank_in_cell"]


def _empty_selected_frame() -> pd.DataFrame:
    return pd.DataFrame(columns=_selected_columns())


def _nonempty_ctas_in_order(cta_list: list[str], selected: pd.DataFrame) -> list[str]:
    if selected.empty or "cta" not in selected.columns:
        return []
    nonempty = set(selected["cta"].dropna().astype(str))
    return [cta for cta in cta_list if cta in nonempty]


def _signature_value(value: object) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def _selected_pmhc_signature(selected: pd.DataFrame, cta: str) -> tuple[tuple[str, ...], ...]:
    if selected.empty or "cta" not in selected.columns:
        return ()
    cta_selected = selected[selected["cta"] == cta]
    if cta_selected.empty:
        return ()
    columns = [column for column in _PMHC_SIGNATURE_COLUMNS if column in cta_selected.columns]
    records = []
    for row in cta_selected[columns].itertuples(index=False, name=None):
        records.append(tuple(_signature_value(value) for value in row))
    return tuple(sorted(records))


def _split_member_text(value: object) -> list[str]:
    if not isinstance(value, str):
        return []
    return [part.strip() for part in value.split(";") if part.strip()]


def _selected_first_text(selected: pd.DataFrame, cta: str, column: str) -> str:
    if selected.empty or column not in selected.columns:
        return ""
    values = selected.loc[selected["cta"] == cta, column].dropna().astype(str)
    if values.empty:
        return ""
    return str(values.iloc[0])


def _selected_member_labels(selected: pd.DataFrame, cta: str) -> list[str]:
    members = _split_member_text(_selected_first_text(selected, cta, "cta_peptide_set_members"))
    return members or [cta]


def _cta_group_label(members: tuple[str, ...]) -> str:
    if len(members) == 1:
        return members[0]
    return "/".join(members)


def _cta_group_order(
    cta_list: list[str],
    selected: pd.DataFrame,
    group_identical_cta_pmhcs: bool,
    group_identical_cta_peptide_sets: bool,
) -> tuple[list[str], list[dict[str, object]]]:
    if not cta_list:
        return [], []

    groups: list[list[str]] = []
    signature_to_group_index: dict[tuple[str, object], int] = {}
    for cta in cta_list:
        peptide_set_label = _selected_first_text(selected, cta, "cta_peptide_set_label")
        peptide_set_members = _split_member_text(
            _selected_first_text(selected, cta, "cta_peptide_set_members")
        )
        if group_identical_cta_peptide_sets and peptide_set_label and len(peptide_set_members) > 1:
            signature = ("peptide_set", peptide_set_label)
        else:
            pmhc_signature = _selected_pmhc_signature(selected, cta)
            signature = ("pmhc", pmhc_signature) if group_identical_cta_pmhcs else ("cta", cta)
        if signature[1]:
            group_index = signature_to_group_index.get(signature)
            if group_index is not None:
                groups[group_index].append(cta)
                continue
            signature_to_group_index[signature] = len(groups)
        groups.append([cta])

    grouped = []
    labels = []
    for source_ctas in groups:
        source_tuple = tuple(source_ctas)
        peptide_set_label = _selected_first_text(selected, source_ctas[0], "cta_peptide_set_label")
        peptide_set_members = _split_member_text(
            _selected_first_text(selected, source_ctas[0], "cta_peptide_set_members")
        )
        if peptide_set_label and len(peptide_set_members) > 1:
            label = peptide_set_label
        else:
            label = _cta_group_label(source_tuple)
        members: list[str] = []
        seen_members: set[str] = set()
        for source_cta in source_ctas:
            for member in _selected_member_labels(selected, source_cta):
                if member not in seen_members:
                    members.append(member)
                    seen_members.add(member)
        labels.append(label)
        grouped.append({"cta": label, "members": members, "source_ctas": list(source_tuple)})
    return labels, grouped


def _cta_group_source_set(cta_groups: list[dict[str, object]]) -> set[str]:
    source_ctas: set[str] = set()
    for group in cta_groups:
        values = group.get("source_ctas", group.get("members", []))
        for source_cta in values:
            source_ctas.add(str(source_cta))
    return source_ctas


def _apply_cta_groups(selected: pd.DataFrame, cta_groups: list[dict[str, object]]) -> pd.DataFrame:
    out = selected.copy()
    for column in [*_LONG_OUTPUT_COLUMNS, *_INTERNAL_SELECTED_COLUMNS]:
        if column not in out.columns:
            out[column] = pd.NA
    if "peptide_rank_in_cell" not in out.columns:
        out["peptide_rank_in_cell"] = pd.NA
    if out.empty or "cta" not in out.columns:
        return out[[*_LONG_OUTPUT_COLUMNS, "peptide_rank_in_cell"]]

    source_to_label: dict[str, str] = {}
    source_to_members: dict[str, str] = {}
    for group in cta_groups:
        label = str(group["cta"])
        members = [str(member) for member in group.get("members", [])]
        member_text = ";".join(members) if members else label
        for source_cta in group.get("source_ctas", members):
            source_to_label[str(source_cta)] = label
            source_to_members[str(source_cta)] = member_text

    out["cta_members"] = out["cta"].map(source_to_members).fillna(out["cta"])
    out["cta"] = out["cta"].map(source_to_label).fillna(out["cta"])
    dedupe_columns = [
        column
        for column in [*_LONG_OUTPUT_COLUMNS, "peptide_rank_in_cell"]
        if column in out.columns
    ]
    out = out.drop_duplicates(subset=dedupe_columns).reset_index(drop=True)
    return out[[*_LONG_OUTPUT_COLUMNS, "peptide_rank_in_cell"]]


def _cta_rank_values_with_groups(
    cta_rank_values: dict[str, float | str],
    cta_groups: list[dict[str, object]],
) -> dict[str, float | str]:
    out = dict(cta_rank_values)
    for group in cta_groups:
        label = str(group["cta"])
        values = [
            cta_rank_values[source_cta]
            for source_cta in group.get("source_ctas", group.get("members", []))
            if source_cta in cta_rank_values
        ]
        numeric_values = []
        for value in values:
            try:
                numeric_values.append(float(value))
            except (TypeError, ValueError):
                continue
        if numeric_values:
            out[label] = max(numeric_values)
        elif values:
            out[label] = values[0]
    return out


def _predictor_key(predictor: str) -> str:
    return predictor.lower().replace("-", "_")


def _ensure_score_columns(scores: pd.DataFrame) -> pd.DataFrame:
    out = scores.copy()
    for column in (
        "peptide",
        "allele",
        "presentation_score",
        "presentation_percentile",
        "affinity_nm",
        "affinity_percentile",
    ):
        if column not in out.columns:
            out[column] = pd.NA
    return out


def _annotate_netmhcpan_affinity(
    selected: pd.DataFrame,
    predictor: str,
    on_progress: Callable[[str], None] | None,
    progress_bar: bool,
    score_chunk_size: int | None,
    progress_file: TextIO | None,
) -> pd.DataFrame:
    out = selected.copy()
    for column in ("netmhcpan_affinity_nm", "netmhcpan_affinity_percentile"):
        if column not in out.columns:
            out[column] = pd.NA
    if out.empty:
        return out

    if _predictor_key(predictor) in {"netmhcpan", "netmhcpan_ba"}:
        out["netmhcpan_affinity_nm"] = out.get("affinity_nm", pd.NA)
        out["netmhcpan_affinity_percentile"] = out.get("affinity_percentile", pd.NA)
        return out

    peptide_list = out["peptide"].dropna().astype(str).drop_duplicates().tolist()
    allele_list = out["allele"].dropna().astype(str).drop_duplicates().tolist()
    n_predictions = len(peptide_list) * len(allele_list)
    _report_progress(
        on_progress,
        f"Annotating {len(out)} selected pMHC row(s) with NetMHCpan BA "
        f"({n_predictions} prediction grid entries)...",
    )
    scoring_start = time.perf_counter()
    scores = _ensure_score_columns(
        _score_presentations(
            peptides=peptide_list,
            alleles=allele_list,
            predictor="netmhcpan",
            progress_bar=progress_bar,
            score_chunk_size=score_chunk_size,
            progress_file=progress_file,
        )
    )
    scoring_elapsed = time.perf_counter() - scoring_start
    _report_progress(on_progress, f"NetMHCpan BA annotation finished in {scoring_elapsed:.1f}s.")

    annotation = (
        scores[["peptide", "allele", "affinity_nm", "affinity_percentile"]]
        .drop_duplicates(subset=["peptide", "allele"])
        .rename(
            columns={
                "affinity_nm": "netmhcpan_affinity_nm",
                "affinity_percentile": "netmhcpan_affinity_percentile",
            }
        )
    )
    out["_row_order"] = range(len(out))
    out = out.drop(columns=["netmhcpan_affinity_nm", "netmhcpan_affinity_percentile"]).merge(
        annotation,
        on=["peptide", "allele"],
        how="left",
    )
    out = out.sort_values("_row_order").drop(columns=["_row_order"]).reset_index(drop=True)
    return out


def _select_cta_batch(
    cta_list: list[str],
    allele_list: list[str],
    lengths: tuple[int, ...],
    ensembl_release: int,
    require_cta_exclusive: bool,
    iedb_path: str | Path | None,
    cedar_path: str | Path | None,
    predictor: str,
    cutoffs: dict[str, float],
    include_predicted_only: bool,
    peptides_per_cell: int,
    group_identical_cta_peptide_sets: bool,
    on_progress: Callable[[str], None] | None = None,
    progress_bar: bool = False,
    score_chunk_size: int | None = None,
    progress_file: TextIO | None = None,
) -> pd.DataFrame:
    if not cta_list:
        return _empty_selected_frame()

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
        return _empty_selected_frame()
    pep_df = _annotate_identical_cta_peptide_sets(
        pep_df,
        cta_list=cta_list,
        group_identical_cta_peptide_sets=group_identical_cta_peptide_sets,
        on_progress=on_progress,
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
    scores = _ensure_score_columns(
        _score_presentations(
            peptides=unique_peptides,
            alleles=score_alleles,
            predictor=predictor,
            progress_bar=progress_bar,
            score_chunk_size=score_chunk_size,
            progress_file=progress_file,
        )
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
        cutoffs=cutoffs,
        include_predicted_only=include_predicted_only,
    )
    return _top_per_cell(candidates, peptides_per_cell=peptides_per_cell)


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
    not overweight that region, then weight regions by population. If an allele
    has no numeric regional proxy, fall back to published global frequencies
    when available so known-frequency panel alleles do not report artificial
    zero coverage.
    """
    return _frequencies_from_audit(_frequency_audit_for_alleles(allele_list), allele_list)


def _frequency_audit_for_alleles(allele_list: list[str]) -> pd.DataFrame:
    try:
        from .regions import allele_frequency_audit
    except ImportError:
        return pd.DataFrame(
            {
                "allele": allele_list,
                "coverage_frequency": [0.0] * len(allele_list),
                "coverage_frequency_source": ["missing"] * len(allele_list),
            }
        )
    return allele_frequency_audit(allele_list)


def _frequencies_from_audit(audit: pd.DataFrame, allele_list: list[str]) -> dict[str, float]:
    frequencies: dict[str, float] = dict.fromkeys(allele_list, 0.0)
    if audit.empty or not {"allele", "coverage_frequency"} <= set(audit.columns):
        return frequencies
    for row in audit.itertuples(index=False):
        frequency = getattr(row, "coverage_frequency", 0.0)
        if pd.notna(frequency):
            frequencies[str(row.allele)] = float(frequency)
    return frequencies


def _frequency_sources_from_audit(audit: pd.DataFrame, allele_list: list[str]) -> dict[str, str]:
    sources: dict[str, str] = dict.fromkeys(allele_list, "missing")
    if audit.empty or not {"allele", "coverage_frequency_source"} <= set(audit.columns):
        return sources
    for row in audit.itertuples(index=False):
        sources[str(row.allele)] = str(getattr(row, "coverage_frequency_source", "missing"))
    return sources


def _order_alleles_by_frequency(
    allele_list: list[str],
    allele_frequencies: dict[str, float],
) -> list[str]:
    return sorted(allele_list, key=lambda allele: (-allele_frequencies.get(allele, 0.0), allele))


def _output_cta_list(
    cta_list: list[str],
    selected: pd.DataFrame,
    include_empty_ctas: bool,
) -> tuple[list[str], list[str]]:
    if include_empty_ctas:
        return list(cta_list), []
    if selected.empty or "cta" not in selected.columns:
        return [], list(cta_list)
    nonempty = set(selected["cta"].dropna().unique())
    output = [cta for cta in cta_list if cta in nonempty]
    empty = [cta for cta in cta_list if cta not in nonempty]
    return output, empty


def _resolve_ctas(
    ctas: Iterable[str] | None,
    cta_count: int | None,
    cta_rank_by: str,
    min_restriction_confidence: Iterable[str] | None,
    restriction_levels: Iterable[str] | None,
    selection_allowlist: Iterable[str] | None = _DEFAULT_SELECTION_ALLOWLIST,
    exclude_vital_tissue_expression: bool = True,
    vital_tissue_max_ntpm: float = _DEFAULT_VITAL_TISSUE_MAX_NTPM,
    exclude_non_magea4_mage_family: bool = _DEFAULT_EXCLUDE_NON_MAGEA4_MAGE_FAMILY,
) -> list[str]:
    """Pick the CTA gene set per the resolution order: explicit override
    wins; else filter the bundled CTA CSV and take top-N by rank. If
    ``cta_count`` is ``None``, return all ranked automatic candidates."""
    from .gene_sets import CTA_gene_names
    from .loader import cta_dataframe, passes_filters_mask

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
    df = df[passes_filters_mask(df)]
    if "never_expressed" in df.columns:
        df = df[~(df["never_expressed"].astype(str).str.lower() == "true")]

    if min_restriction_confidence is not None and "restriction_confidence" in df.columns:
        wanted = {c.upper() for c in min_restriction_confidence}
        confidence_ok = df["restriction_confidence"].astype(str).str.upper().isin(wanted)
        df = df[confidence_ok | df["_is_selection_allowlisted"]]

    if restriction_levels is not None and "restriction" in df.columns:
        wanted_levels = set(restriction_levels)
        df = df[df["restriction"].isin(wanted_levels)]

    if exclude_non_magea4_mage_family:
        non_magea4_mage = df["_cta_display"].map(_is_non_magea4_mage_family)
        df = df[~non_magea4_mage | df["_is_selection_allowlisted"]]

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
    if allowlist:
        allowlisted = df[df["_cta_display"].isin(allowlist)]
        ranked = df[~df["_cta_display"].isin(allowlist)]
        df = pd.concat([allowlisted, ranked], ignore_index=True)
    if cta_count is None:
        return df["_cta_display"].tolist()
    return df["_cta_display"].head(cta_count).tolist()


def _compact_cta_name(value: object) -> str:
    return "".join(ch for ch in str(value).upper() if ch.isalnum())


def _cta_display_name(symbol: object) -> str:
    token = str(symbol).strip()
    compact = _compact_cta_name(token)
    if compact in _CTAG1_ALIASES:
        return _CTAG1_GROUP_LABEL
    if compact.startswith("MAGEA") and compact[5:].isdigit():
        return compact
    return token


def _is_non_magea4_mage_family(label: object) -> bool:
    compact = _compact_cta_name(label)
    return compact.startswith("MAGE") and compact != "MAGEA4"


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
        if compact in _CTAG1_ALIASES:
            label = _CTAG1_GROUP_LABEL
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
    if len(cta_symbols) != len(cta_list):
        _report_progress(
            on_progress,
            f"Peptide enumeration expands {len(cta_list)} CTA target labels to "
            f"{len(cta_symbols)} Ensembl genes after alias/group expansion.",
        )
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


def _peptide_set_signature(pep_df: pd.DataFrame, cta: str) -> tuple[str, ...]:
    if pep_df.empty or "gene_name" not in pep_df.columns or "peptide" not in pep_df.columns:
        return ()
    peptides = pep_df.loc[pep_df["gene_name"] == cta, "peptide"].dropna().astype(str)
    return tuple(sorted(peptides.unique()))


def _annotate_identical_cta_peptide_sets(
    pep_df: pd.DataFrame,
    cta_list: list[str],
    group_identical_cta_peptide_sets: bool,
    on_progress: Callable[[str], None] | None = None,
) -> pd.DataFrame:
    out = pep_df.copy()
    if "gene_name" not in out.columns:
        out["cta_peptide_set_label"] = pd.NA
        out["cta_peptide_set_members"] = pd.NA
        return out

    label_by_cta = {cta: cta for cta in cta_list}
    members_by_cta = {cta: cta for cta in cta_list}
    if group_identical_cta_peptide_sets:
        signature_to_ctas: dict[tuple[str, ...], list[str]] = {}
        for cta in cta_list:
            signature = _peptide_set_signature(out, cta)
            if signature:
                signature_to_ctas.setdefault(signature, []).append(cta)

        grouped_labels = []
        for members in signature_to_ctas.values():
            if len(members) < 2:
                continue
            member_tuple = tuple(members)
            group_label = _cta_group_label(member_tuple)
            member_text = ";".join(members)
            grouped_labels.append(group_label)
            for member in members:
                label_by_cta[member] = group_label
                members_by_cta[member] = member_text

        if grouped_labels:
            _report_progress(
                on_progress,
                "Grouped CTA targets with identical peptide sets: "
                + ", ".join(sorted(grouped_labels)),
            )

    out["cta_peptide_set_label"] = out["gene_name"].map(label_by_cta).fillna(out["gene_name"])
    out["cta_peptide_set_members"] = out["gene_name"].map(members_by_cta).fillna(out["gene_name"])
    return out


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
        sample_alleles = set(split_mhc_restrictions(row.get("mhc_allele_set", "")))
        provenance = str(row.get("mhc_allele_provenance", "")).strip()
        matched = False

        if is_monoallelic and exact_panel_alleles:
            for allele in exact_panel_alleles:
                _add_evidence(stats, (peptide, allele, "monoallelic_ms"), row)
            matched = True
        elif provenance == "sample_allele_match" and sample_alleles:
            for allele in sorted(panel_alleles & sample_alleles):
                if _is_best_sample_allele(peptide, allele, sample_alleles, score_lookup):
                    _add_evidence(stats, (peptide, allele, "sample_allele_ms"), row)
                    matched = True
        elif exact_panel_alleles:
            for allele in exact_panel_alleles:
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
    columns = [*_LONG_OUTPUT_COLUMNS, *_INTERNAL_SELECTED_COLUMNS, "_tier_rank"]
    if scores.empty or "presentation_percentile" not in scores.columns:
        return pd.DataFrame(columns=columns)

    pep_df = pep_df.copy()
    for column in _INTERNAL_SELECTED_COLUMNS:
        if column not in pep_df.columns:
            pep_df[column] = pep_df.get("gene_name", pd.NA)
    pep_meta = (
        pep_df[["peptide", "gene_name", "length", *_INTERNAL_SELECTED_COLUMNS]]
        .drop_duplicates(subset=["peptide", "gene_name", "length", *_INTERNAL_SELECTED_COLUMNS])
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
                "cta_members": row.cta,
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
                "affinity_percentile": getattr(row, "affinity_percentile", pd.NA),
                "netmhcpan_affinity_nm": pd.NA,
                "netmhcpan_affinity_percentile": pd.NA,
                "cta_peptide_set_label": getattr(row, "cta_peptide_set_label", row.cta),
                "cta_peptide_set_members": getattr(row, "cta_peptide_set_members", row.cta),
                "_tier_rank": _EVIDENCE_TIER_RANKS[chosen_tier],
            }
        )

    return pd.DataFrame(rows, columns=columns)


def _top_per_cell(candidates: pd.DataFrame, peptides_per_cell: int) -> pd.DataFrame:
    """For each (cta, allele), keep the top-N peptides by evidence support."""
    if candidates.empty:
        return pd.DataFrame(columns=_selected_columns())

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
    columns = _selected_columns()
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


def _allele_locus(allele: object) -> str:
    token = str(allele).strip().upper()
    if token.startswith("HLA-"):
        token = token.removeprefix("HLA-")
    if "*" in token:
        locus = token.split("*", 1)[0]
    elif token and token[0] in {"A", "B", "C"}:
        locus = token[0]
    else:
        locus = token
    return locus or str(allele)


def _combined_carrier_probability(allele_frequencies: dict[str, float]) -> float:
    locus_frequencies: dict[str, float] = defaultdict(float)
    for allele, frequency in allele_frequencies.items():
        locus_frequencies[_allele_locus(allele)] += max(0.0, min(float(frequency), 1.0))

    miss_probability = 1.0
    for frequency in locus_frequencies.values():
        miss_probability *= 1.0 - _carrier_probability(float(frequency))
    return max(0.0, min(1.0, 1.0 - miss_probability))


def _selected_tier_row_counts(selected: pd.DataFrame) -> dict[str, int]:
    counts = dict.fromkeys(_EVIDENCE_TIER_RANKS, 0)
    if selected.empty or "evidence_tier" not in selected.columns:
        return counts
    counts.update(
        {
            str(tier): int(count)
            for tier, count in selected["evidence_tier"].value_counts(sort=False).items()
        }
    )
    return counts


def panel_summary(
    selected: pd.DataFrame,
    cta_list: list[str],
    allele_list: list[str],
    allele_frequencies: dict[str, float] | None = None,
    allele_frequency_sources: dict[str, str] | None = None,
) -> dict[str, object]:
    """Summarize a selected panel table.

    Population coverage estimates use weighted regional allele frequencies,
    sum covered allele frequencies within each HLA locus, convert each locus
    sum to a carrier probability, then combine loci. The values are intended
    for ranking and sanity-checking panel breadth, not for clinical-grade
    population-genetics inference.
    """
    allele_frequencies = allele_frequencies or dict.fromkeys(allele_list, 0.0)
    allele_frequency_sources = allele_frequency_sources or dict.fromkeys(allele_list, "missing")
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
        if not cta_selected.empty and "cta_members" in cta_selected.columns:
            cta_members = str(cta_selected["cta_members"].dropna().astype(str).iloc[0])
        else:
            cta_members = cta
        tier_row_counts = _selected_tier_row_counts(cta_selected)
        coverage = _combined_carrier_probability(
            {allele: allele_frequencies.get(allele, 0.0) for allele in covered_alleles}
        )
        cta_rows.append(
            {
                "cta": cta,
                "cta_members": cta_members,
                "covered_hla_count": len(covered_alleles),
                "selected_peptide_count": int(cta_selected["peptide"].nunique())
                if not cta_selected.empty
                else 0,
                "monoallelic_ms_pmhc_count": tier_row_counts.get("monoallelic_ms", 0),
                "sample_allele_ms_pmhc_count": tier_row_counts.get("sample_allele_ms", 0),
                "unrestricted_ms_pmhc_count": tier_row_counts.get("unrestricted_ms", 0),
                "predicted_only_pmhc_count": tier_row_counts.get("predicted_only", 0),
                "estimated_population_coverage": coverage,
            }
        )
    cta_rows = sorted(
        cta_rows,
        key=lambda row: (
            -int(row["selected_peptide_count"]),
            -int(row["covered_hla_count"]),
            -float(row["estimated_population_coverage"]),
            str(row["cta"]),
        ),
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
                "frequency_source": allele_frequency_sources.get(allele, "missing"),
                "covered_cta_count": len(covered_ctas),
                "covered_cta_fraction": 0.0
                if total_ctas == 0
                else float(len(covered_ctas) / total_ctas),
                "selected_peptide_count": int(allele_selected["peptide"].nunique())
                if not allele_selected.empty
                else 0,
            }
        )

    tier_counts = {
        tier: count for tier, count in _selected_tier_row_counts(selected).items() if count
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
            "Estimated from 0-1 HLA allele frequencies. Coverage uses population-weighted "
            "regional proxy frequencies when available and published global averages only "
            "for alleles with no numeric regional proxy. Covered allele frequencies are "
            "summed within each HLA locus, converted to carrier probability per locus, "
            "then combined across loci."
        ),
    }


def _attach_panel_attrs(
    out: pd.DataFrame,
    cta_list: list[str],
    allele_list: list[str],
    allele_frequencies: dict[str, float],
    cta_rank_values: dict[str, float | str],
    summary: dict[str, object],
    allele_frequency_sources: dict[str, str] | None = None,
    frequency_audit: pd.DataFrame | None = None,
    input_cta_list: list[str] | None = None,
    empty_ctas: list[str] | None = None,
    cta_groups: list[dict[str, object]] | None = None,
) -> pd.DataFrame:
    out.attrs["cta_order"] = list(cta_list)
    out.attrs["input_cta_order"] = list(input_cta_list or cta_list)
    out.attrs["empty_ctas"] = list(empty_ctas or [])
    out.attrs["cta_groups"] = _public_cta_groups(cta_groups or [])
    out.attrs["allele_order"] = list(allele_list)
    out.attrs["allele_frequencies"] = dict(allele_frequencies)
    out.attrs["allele_frequency_sources"] = dict(allele_frequency_sources or {})
    out.attrs["allele_frequency_audit"] = (
        [] if frequency_audit is None else frequency_audit.to_dict("records")
    )
    out.attrs["cta_rank_values"] = dict(cta_rank_values)
    out.attrs["panel_summary"] = summary
    return out


def _public_cta_groups(cta_groups: list[dict[str, object]]) -> list[dict[str, object]]:
    return [
        {"cta": str(group["cta"]), "members": [str(member) for member in group.get("members", [])]}
        for group in cta_groups
    ]


def _empty_output(
    cta_list: list[str],
    allele_list: list[str],
    output_format: str,
    allele_frequencies: dict[str, float] | None = None,
    cta_rank_values: dict[str, float | str] | None = None,
    allele_frequency_sources: dict[str, str] | None = None,
    frequency_audit: pd.DataFrame | None = None,
    input_cta_list: list[str] | None = None,
    empty_ctas: list[str] | None = None,
    include_empty_ctas: bool = True,
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
        allele_frequency_sources=allele_frequency_sources,
    )
    input_cta_list = input_cta_list or cta_list
    empty_ctas = empty_ctas or []
    summary["input_cta_count"] = len(input_cta_list)
    summary["empty_cta_count"] = len(empty_ctas)
    summary["empty_ctas"] = empty_ctas
    summary["include_empty_ctas"] = include_empty_ctas
    cta_groups = [{"cta": cta, "members": [cta]} for cta in cta_list]
    summary["group_identical_cta_pmhcs"] = False
    summary["group_identical_cta_peptide_sets"] = False
    summary["cta_groups"] = _public_cta_groups(cta_groups)
    summary["grouped_cta_member_count"] = 0
    summary["annotate_netmhcpan_affinity"] = False
    return _attach_panel_attrs(
        out,
        cta_list,
        allele_list,
        allele_frequencies,
        cta_rank_values or {},
        summary,
        allele_frequency_sources,
        frequency_audit,
        input_cta_list=input_cta_list,
        empty_ctas=empty_ctas,
        cta_groups=cta_groups,
    )
