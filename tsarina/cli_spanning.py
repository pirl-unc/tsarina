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

"""``tsarina panel`` CLI — CTA x HLA pMHC panel matrix.

Wraps :func:`tsarina.spanning.spanning_pmhc_set` with argparse.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_SUPPORTED_PREDICTORS = ("mhcflurry", "netmhcpan", "netmhcpan_el")
_SUPPORTED_FORMATS = ("table", "wide", "long")
_SUPPORTED_PANELS = (
    "iedb27_ab",
    "iedb36_abc",
    "global44_abc",
    "global48_abc",
    "global51_abc_ssa",
    "global51_abc",
    "global53_abc",
)
_DEFAULT_PANEL = "global53_abc"
_DEFAULT_LENGTHS = (8, 9, 10, 11)
_DEFAULT_SELECTION_ALLOWLIST = "PRAME,NY-ESO-1,MAGEA4"
_DEFAULT_VITAL_TISSUE_MAX_NTPM = 2.0


def _split_csv(value: str) -> list[str]:
    return [s.strip() for s in value.split(",") if s.strip()]


def _parse_lengths(value: str) -> tuple[int, ...]:
    try:
        return tuple(int(x) for x in _split_csv(value))
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"--lengths must be integers: {value!r}") from e


def _configure_parser(p: argparse.ArgumentParser) -> argparse.ArgumentParser:
    p.add_argument(
        "--cta-count",
        type=int,
        default=25,
        help=(
            "Maximum downstream non-empty CTAs to include when --ctas is not supplied; "
            "lower-ranked candidates are scanned as needed, while --show-empty-ctas "
            "restores the top-candidate audit view (default 25)."
        ),
    )
    p.add_argument(
        "--cta-rank-by",
        default="ms_cancer_peptide_count",
        help=(
            "Column in the bundled CTA CSV used to rank candidates "
            "(default 'ms_cancer_peptide_count'; falls back to alphabetical "
            "Symbol order if the column is missing or all-NaN)."
        ),
    )
    p.add_argument(
        "--ctas",
        type=_split_csv,
        default=None,
        help=(
            "Explicit gene symbols, comma-separated (e.g. 'MAGEA4,PRAME,NY-ESO-1'). "
            "Overrides --cta-count / --cta-rank-by."
        ),
    )
    p.add_argument(
        "--min-restriction-confidence",
        type=_split_csv,
        default=["HIGH", "MODERATE"],
        help=(
            "Allowed restriction_confidence bins, comma-separated (default "
            "'HIGH,MODERATE'; pass 'ANY' to disable)."
        ),
    )
    p.add_argument(
        "--restriction-levels",
        type=_split_csv,
        default=None,
        help=(
            "Restriction levels to keep, comma-separated (e.g. 'TESTIS,PLACENTAL'). "
            "Default keeps all levels."
        ),
    )
    p.add_argument(
        "--selection-allowlist",
        type=_split_csv,
        default=_split_csv(_DEFAULT_SELECTION_ALLOWLIST),
        help=(
            "CTA names allowed through automatic safety/confidence gates "
            f"(default '{_DEFAULT_SELECTION_ALLOWLIST}'). Aliases such as NY-ESO-1 "
            "and MAGE-A4 are normalized."
        ),
    )
    p.add_argument(
        "--no-vital-tissue-filter",
        dest="exclude_vital_tissue_expression",
        action="store_false",
        help=(
            "Disable the automatic vital-tissue expression filter. By default, panel "
            "excludes CTAs with RNA above --vital-tissue-max-ntpm or unique public "
            "healthy-MS observations in brain/CNS/cerebellum, heart, lung, liver, or "
            "pancreas, unless allowlisted."
        ),
    )
    p.add_argument(
        "--vital-tissue-max-ntpm",
        type=float,
        default=_DEFAULT_VITAL_TISSUE_MAX_NTPM,
        help=(
            "Maximum allowed RNA nTPM in vital tissues for automatic CTA selection "
            f"(default {_DEFAULT_VITAL_TISSUE_MAX_NTPM})."
        ),
    )
    p.add_argument(
        "--allow-non-magea4-mage-family",
        dest="exclude_non_magea4_mage_family",
        action="store_false",
        help=(
            "Allow automatic panel selection to include MAGE-family CTAs other than "
            "MAGE-A4. By default, non-MAGE-A4 MAGE-family CTAs are excluded unless "
            "explicitly requested with --ctas or included in --selection-allowlist."
        ),
    )
    p.add_argument(
        "--alleles",
        type=_split_csv,
        default=None,
        help="Explicit allele list (e.g. 'HLA-A*02:01,HLA-A*24:02'). Overrides --panel.",
    )
    p.add_argument(
        "--panel",
        choices=_SUPPORTED_PANELS,
        default=_DEFAULT_PANEL,
        help=f"Named panel from tsarina.alleles (default '{_DEFAULT_PANEL}').",
    )
    p.add_argument(
        "--lengths",
        type=_parse_lengths,
        default=_DEFAULT_LENGTHS,
        help="Peptide lengths to enumerate (default 8,9,10,11).",
    )
    p.add_argument(
        "--ensembl-release",
        type=int,
        default=112,
        help="Ensembl release (default 112).",
    )
    p.add_argument(
        "--no-require-cta-exclusive",
        dest="require_cta_exclusive",
        action="store_false",
        help="Allow CTA peptides that also appear in non-CTA proteins.",
    )
    p.add_argument(
        "--predictor",
        choices=_SUPPORTED_PREDICTORS,
        default="mhcflurry",
        help="mhctools predictor (default mhcflurry).",
    )
    p.add_argument(
        "--monoallelic-ms-max-percentile",
        type=float,
        default=2.0,
        help="Max percentile for mono-allelic MS-supported pMHCs (default 2.0).",
    )
    p.add_argument(
        "--sample-allele-ms-max-percentile",
        type=float,
        default=1.0,
        help=(
            "Max percentile for multi-allelic sample-genotype MS support where the "
            "panel allele is best among sample alleles (default 1.0)."
        ),
    )
    p.add_argument(
        "--unrestricted-ms-max-percentile",
        type=float,
        default=0.5,
        help="Max percentile for peptide-level MS support without usable allele (default 0.5).",
    )
    p.add_argument(
        "--include-predicted-only",
        action="store_true",
        help="Allow prediction-only candidates as last-priority fillers.",
    )
    p.add_argument(
        "--predicted-only-max-percentile",
        type=float,
        default=0.1,
        help="Max percentile for prediction-only candidates when enabled (default 0.1).",
    )
    p.add_argument(
        "--max-percentile",
        type=float,
        default=None,
        help=(
            "Deprecated shortcut: set all MS-supported tier cutoffs to this value. "
            "Does not enable predicted-only candidates."
        ),
    )
    p.add_argument(
        "--peptides-per-cell",
        type=int,
        default=3,
        help=(
            "Maximum peptides to keep per CTA x HLA cell, ranked by MS source count, "
            "MS hit count, then prediction percentile (default 3)."
        ),
    )
    p.add_argument(
        "--iedb-path",
        default=None,
        help="Optional explicit IEDB ligand CSV; default uses cached hitlist observations.",
    )
    p.add_argument(
        "--cedar-path",
        default=None,
        help="Optional explicit CEDAR CSV; default uses cached hitlist observations.",
    )
    p.add_argument(
        "--format",
        choices=_SUPPORTED_FORMATS,
        default="table",
        help=(
            "Output format (default 'table'). 'table' = readable terminal report; "
            "'wide' = CTA rows x allele cols with selected peptides as cell values; "
            "'long' = one row per selected peptide with evidence tier / MS provenance / "
            "percentile / score / affinity columns."
        ),
    )
    p.add_argument(
        "--no-summary",
        dest="summary",
        action="store_false",
        help="Do not append the panel coverage summary.",
    )
    p.add_argument(
        "--show-empty-ctas",
        action="store_true",
        help=(
            "Show automatically selected CTAs that have no selected pMHCs after "
            "peptide, exclusivity, public-MS, and prediction gates. Explicit --ctas "
            "requests are always preserved."
        ),
    )
    p.add_argument(
        "--no-progress",
        dest="progress",
        action="store_false",
        help="Suppress progress messages and scoring progress bars.",
    )
    p.add_argument(
        "--progress-bars",
        dest="progress_bars",
        action="store_true",
        default=None,
        help="Force tqdm progress bars for chunked scoring.",
    )
    p.add_argument(
        "--no-progress-bars",
        dest="progress_bars",
        action="store_false",
        help="Disable tqdm progress bars while keeping progress messages.",
    )
    p.add_argument(
        "--score-chunk-size",
        type=int,
        default=None,
        help=(
            "Number of HLA alleles per scoring chunk when progress bars are enabled "
            "(default 8 for MHCflurry, 1 for other predictors)."
        ),
    )
    p.add_argument(
        "-o",
        "--output",
        default=None,
        help="Write output to this path (CSV for wide/long, text for table; default stdout).",
    )
    return p


def build_parser(sub: argparse._SubParsersAction) -> argparse.ArgumentParser:
    p = sub.add_parser(
        "panel",
        help="Build a CTA x HLA pMHC matrix for a population HLA panel.",
        description=(
            "Produce a CTA x HLA pivot table where each cell is the best MS-supported "
            "peptide-HLA candidate for that CTA and allele. Defaults to up to 25 "
            "non-empty CTAs crossed with the Global-53 HLA-A/B/C panel, 8-11mers, "
            "and tier-specific presentation-percentile cutoffs: mono-allelic MS <2.0, "
            "multi-allelic "
            "sample-genotype MS <1.0, unrestricted MS <0.5. Prediction-only "
            "candidates are excluded unless --include-predicted-only is supplied. "
            "The default output is a readable table plus coverage summary."
        ),
    )
    _configure_parser(p)

    return p


def handle(args: argparse.Namespace) -> None:
    from .spanning import spanning_pmhc_set

    if getattr(args, "deprecated_spanning", False):
        print("Warning: `tsarina spanning` is deprecated; use `tsarina panel`.", file=sys.stderr)

    min_restriction_confidence: tuple[str, ...] | None
    if any(v.upper() == "ANY" for v in args.min_restriction_confidence):
        min_restriction_confidence = None
    else:
        min_restriction_confidence = tuple(v.upper() for v in args.min_restriction_confidence)

    def _on_progress(msg: str) -> None:
        print(msg, file=sys.stderr)

    progress_bar = False
    if args.progress:
        progress_bar = sys.stderr.isatty() if args.progress_bars is None else args.progress_bars

    output_format = "long" if args.format == "table" else args.format
    df = spanning_pmhc_set(
        cta_count=args.cta_count,
        cta_rank_by=args.cta_rank_by,
        ctas=args.ctas,
        min_restriction_confidence=min_restriction_confidence,
        restriction_levels=tuple(args.restriction_levels) if args.restriction_levels else None,
        selection_allowlist=tuple(args.selection_allowlist),
        exclude_vital_tissue_expression=args.exclude_vital_tissue_expression,
        vital_tissue_max_ntpm=args.vital_tissue_max_ntpm,
        exclude_non_magea4_mage_family=args.exclude_non_magea4_mage_family,
        alleles=tuple(args.alleles) if args.alleles else None,
        panel=args.panel,
        lengths=args.lengths,
        ensembl_release=args.ensembl_release,
        require_cta_exclusive=args.require_cta_exclusive,
        predictor=args.predictor,
        monoallelic_ms_max_percentile=args.monoallelic_ms_max_percentile,
        sample_allele_ms_max_percentile=args.sample_allele_ms_max_percentile,
        unrestricted_ms_max_percentile=args.unrestricted_ms_max_percentile,
        include_predicted_only=args.include_predicted_only,
        predicted_only_max_percentile=args.predicted_only_max_percentile,
        max_percentile=args.max_percentile,
        peptides_per_cell=args.peptides_per_cell,
        iedb_path=args.iedb_path,
        cedar_path=args.cedar_path,
        output_format=output_format,
        include_empty_ctas=True if args.ctas is not None else args.show_empty_ctas,
        on_progress=_on_progress if args.progress else None,
        progress_bar=progress_bar,
        score_chunk_size=args.score_chunk_size,
        progress_file=sys.stderr,
    )

    if args.format == "table":
        text = format_panel_table(df)
        if args.summary:
            text = f"{text}\n\n{format_panel_summary(df)}"
        if args.output:
            Path(args.output).write_text(f"{text}\n")
            print(f"Wrote panel report to {args.output}", file=sys.stderr)
        else:
            print(text)
        return

    if args.output:
        df.to_csv(args.output, index=False)
        print(f"Wrote {len(df)} rows to {args.output}", file=sys.stderr)
    else:
        df.to_csv(sys.stdout, index=False)

    if args.summary:
        print(format_panel_summary(df), file=sys.stderr)


def _format_percent(value: object) -> str:
    try:
        return f"{100.0 * float(value):.1f}%"
    except (TypeError, ValueError):
        return "n/a"


def _format_rank_value(value: object) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "" if value is None else str(value)
    if number.is_integer():
        return str(int(number))
    return f"{number:.2f}"


def _format_optional_float(value: object, digits: int) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    return f"{number:.{digits}f}"


def _pmid_count(value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        return "0"
    return str(len([part for part in value.split(";") if part.strip()]))


def _plain_table(headers: list[str], rows: list[list[str]]) -> str:
    widths = [
        max(len(str(header)), *(len(str(row[i])) for row in rows)) if rows else len(str(header))
        for i, header in enumerate(headers)
    ]
    header_line = "  ".join(str(header).ljust(widths[i]) for i, header in enumerate(headers))
    rule_line = "  ".join("-" * width for width in widths)
    body = ["  ".join(str(value).ljust(widths[i]) for i, value in enumerate(row)) for row in rows]
    return "\n".join([header_line, rule_line, *body])


def format_panel_table(df) -> str:
    if df.empty:
        return "No panel peptides selected."

    cta_order = {cta: i for i, cta in enumerate(df.attrs.get("cta_order", []))}
    allele_order = {allele: i for i, allele in enumerate(df.attrs.get("allele_order", []))}
    cta_rank_values = df.attrs.get("cta_rank_values", {})
    allele_frequencies = df.attrs.get("allele_frequencies", {})

    out = df.copy()
    out["_cta_order"] = out["cta"].map(cta_order).fillna(len(cta_order))
    out["_allele_order"] = out["allele"].map(allele_order).fillna(len(allele_order))
    out = out.sort_values(["_cta_order", "_allele_order", "peptide_rank_in_cell"])

    rows: list[list[str]] = []
    previous_cta = None
    previous_cell: tuple[str, str] | None = None
    for row in out.itertuples(index=False):
        cta = str(row.cta)
        allele = str(row.allele)
        cell = (cta, allele)
        show_cta = cta != previous_cta
        show_allele = cell != previous_cell
        rows.append(
            [
                cta if show_cta else "",
                _format_rank_value(cta_rank_values.get(cta)) if show_cta else "",
                allele if show_allele else "",
                _format_percent(allele_frequencies.get(allele)) if show_allele else "",
                str(row.peptide),
                str(row.evidence_tier),
                str(row.ms_source_count),
                str(row.ms_hit_count),
                _pmid_count(row.ms_pmids),
                _format_optional_float(row.presentation_percentile, 3),
                _format_optional_float(row.presentation_score, 3),
                _format_optional_float(row.affinity_nm, 1),
            ]
        )
        previous_cta = cta
        previous_cell = cell

    return _plain_table(
        [
            "CTA",
            "CTA rank",
            "HLA",
            "HLA freq",
            "Peptide",
            "Tier",
            "Sources",
            "MS hits",
            "PMIDs",
            "%Rank",
            "Score",
            "nM",
        ],
        rows,
    )


def format_panel_summary(df) -> str:
    summary = df.attrs.get("panel_summary")
    if not summary:
        return "Summary unavailable."

    lines = [
        "Summary",
        f"  HLA alleles: {summary['hla_allele_count']}",
        f"  CTAs: {summary['cta_count']}",
        f"  Selected unique peptides: {summary['selected_peptide_count']}",
        (
            "  Filled CTA x HLA cells: "
            f"{summary['filled_cell_count']}/{summary['possible_cell_count']} "
            f"({_format_percent(summary['filled_cell_fraction'])})"
        ),
        (
            "  Mean peptides per filled cell: "
            f"{float(summary['average_peptides_per_filled_cell']):.2f}"
        ),
    ]

    tier_counts = summary.get("evidence_tier_counts", {})
    if tier_counts:
        lines.append(
            "  Evidence tiers: "
            + ", ".join(f"{tier}={count}" for tier, count in sorted(tier_counts.items()))
        )
    empty_count = int(summary.get("empty_cta_count", 0))
    if empty_count and not summary.get("include_empty_ctas", True):
        lines.append(f"  Omitted CTAs with no selected pMHCs: {empty_count}")
    lines.append(f"  Coverage note: {summary['coverage_note']}")

    cta_rows = [
        [
            str(row["cta"]),
            str(row["covered_hla_count"]),
            str(row["selected_peptide_count"]),
            str(row.get("monoallelic_ms_pmhc_count", 0)),
            str(row.get("sample_allele_ms_pmhc_count", 0)),
            str(row.get("unrestricted_ms_pmhc_count", 0)),
            _format_percent(row["estimated_population_coverage"]),
        ]
        for row in summary["cta_coverage"]
    ]
    lines.extend(
        [
            "",
            "Expected Population Coverage Per CTA",
            _plain_table(
                [
                    "CTA",
                    "HLA hits",
                    "Peptides",
                    "Mono MS",
                    "Sample MS",
                    "Unres MS",
                    "Est. coverage",
                ],
                cta_rows,
            ),
        ]
    )

    hla_rows = [
        [
            str(row["allele"]),
            _format_percent(row["weighted_allele_frequency"]),
            str(row.get("frequency_source", "missing")),
            str(row["covered_cta_count"]),
            _format_percent(row["covered_cta_fraction"]),
            str(row["selected_peptide_count"]),
        ]
        for row in summary["hla_coverage"]
    ]
    lines.extend(
        [
            "",
            "CTA Coverage Per HLA",
            _plain_table(
                ["HLA", "HLA freq", "Freq source", "CTAs", "CTA frac", "Peptides"],
                hla_rows,
            ),
        ]
    )
    return "\n".join(lines)
