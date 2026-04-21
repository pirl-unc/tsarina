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

"""``tsarina spanning`` CLI — population-spanning CTA x HLA pMHC table.

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


def _split_csv(value: str) -> list[str]:
    return [s.strip() for s in value.split(",") if s.strip()]


def _parse_lengths(value: str) -> tuple[int, ...]:
    try:
        return tuple(int(x) for x in _split_csv(value))
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"--lengths must be integers: {value!r}") from e


def build_parser(sub: argparse._SubParsersAction) -> argparse.ArgumentParser:
    p = sub.add_parser(
        "spanning",
        help="Build a population-spanning CTA x HLA pMHC table.",
        description=(
            "Produce a CTA x HLA pivot table where each cell is the best-presented "
            "peptide for that CTA on that HLA allele. Defaults to top-25 CTAs (ranked "
            "by ms_cancer_peptide_count) crossed with the IEDB-27 panel at 9-mer "
            "length and a 2.0%% presentation-percentile cutoff — typically yielding "
            "200-400 filled cells, suitable as a baseline spanning set for off-the-"
            "shelf vaccine / TCR programs."
        ),
    )
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
        default="iedb27_ab",
        help="Named panel from tsarina.alleles (default 'iedb27_ab').",
    )
    p.add_argument(
        "--lengths",
        type=_parse_lengths,
        default=(9,),
        help="Peptide lengths to enumerate (default 9).",
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
        "--max-percentile",
        type=float,
        default=2.0,
        help="Maximum presentation percentile for a cell to be filled (default 2.0).",
    )
    p.add_argument(
        "--format",
        choices=_SUPPORTED_FORMATS,
        default="wide",
        help=(
            "Output format (default 'wide'). 'wide' = CTA rows x allele cols with "
            "peptide as cell value. 'long' = one row per filled cell with peptide / "
            "length / percentile / score / affinity columns."
        ),
    )
    p.add_argument(
        "-o",
        "--output",
        default=None,
        help="Write CSV to this path (default: stdout).",
    )
    return p


def handle(args: argparse.Namespace) -> None:
    from .spanning import spanning_pmhc_set

    min_restriction_confidence: tuple[str, ...] | None
    if any(v.upper() == "ANY" for v in args.min_restriction_confidence):
        min_restriction_confidence = None
    else:
        min_restriction_confidence = tuple(v.upper() for v in args.min_restriction_confidence)

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
        max_percentile=args.max_percentile,
        output_format=args.format,
    )

    if args.output:
        df.to_csv(args.output, index=False)
        print(f"Wrote {len(df)} rows to {args.output}", file=sys.stderr)
    else:
        df.to_csv(sys.stdout, index=False)
