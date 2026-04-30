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

_SUPPORTED_PREDICTORS = ("mhcflurry", "netmhcpan", "netmhcpan_el")
_SUPPORTED_FORMATS = ("wide", "long")
_SUPPORTED_PANELS = (
    "iedb27_ab",
    "iedb36_abc",
    "global44_abc",
    "global48_abc",
    "global51_abc_ssa",
)
_DEFAULT_PANEL = "global51_abc_ssa"
_DEFAULT_LENGTHS = (8, 9, 10, 11)


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
        help="Top N CTAs to include when --ctas is not supplied (default 25).",
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
        default="wide",
        help=(
            "Output format (default 'wide'). 'wide' = CTA rows x allele cols with "
            "peptide as cell value. 'long' = one row per filled cell with peptide / "
            "evidence tier / MS provenance / percentile / score / affinity columns."
        ),
    )
    p.add_argument(
        "-o",
        "--output",
        default=None,
        help="Write CSV to this path (default: stdout).",
    )
    return p


def build_parser(sub: argparse._SubParsersAction) -> argparse.ArgumentParser:
    p = sub.add_parser(
        "panel",
        help="Build a CTA x HLA pMHC matrix for a population HLA panel.",
        description=(
            "Produce a CTA x HLA pivot table where each cell is the best MS-supported "
            "peptide-HLA candidate for that CTA and allele. Defaults to top-25 CTAs "
            "crossed with the Global-51 HLA-A/B/C panel, 8-11mers, and tier-specific "
            "presentation-percentile cutoffs: mono-allelic MS <=2.0, multi-allelic "
            "sample-genotype MS <=1.0, unrestricted MS <=0.5. Prediction-only "
            "candidates are excluded unless --include-predicted-only is supplied."
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

    df = spanning_pmhc_set(
        cta_count=args.cta_count,
        cta_rank_by=args.cta_rank_by,
        ctas=args.ctas,
        min_restriction_confidence=min_restriction_confidence,
        restriction_levels=tuple(args.restriction_levels) if args.restriction_levels else None,
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
        iedb_path=args.iedb_path,
        cedar_path=args.cedar_path,
        output_format=args.format,
        on_progress=_on_progress,
    )

    if args.output:
        df.to_csv(args.output, index=False)
        print(f"Wrote {len(df)} rows to {args.output}", file=sys.stderr)
    else:
        df.to_csv(sys.stdout, index=False)
