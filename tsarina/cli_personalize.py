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

"""`tsarina personalize` CLI — patient-level target prioritization.

Wraps :func:`tsarina.personalize.personalized_targets` with argparse.  IEDB/CEDAR paths
default to whatever the user registered with ``tsarina data register …``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .cli_common import add_iedb_cedar_args, add_predictor_arg
from .cli_common import flatten_multi as _flatten_multi
from .cli_common import parse_lengths as _parse_lengths
from .cli_common import split_csv as _split_csv


def _parse_cta(raw: list[str]) -> dict[str, float]:
    """Parse ``--cta`` tokens (already collected by ``nargs="+"``) into a
    GENE -> TPM dict. Accepts ``GENE=TPM`` entries mixed freely across
    comma-separated, space-separated, or quoted-string tokens.

    The ``=TPM`` half is optional -- a bare gene name maps to ``NaN``,
    which :func:`tsarina.personalize.personalized_targets` treats as "no
    measured expression given, include it regardless of --min-cta-tpm"
    rather than as zero expression."""
    out: dict[str, float] = {}
    for pair in _flatten_multi(raw):
        gene, sep, tpm_s = pair.partition("=")
        if not sep:
            out[gene.strip()] = float("nan")
            continue
        try:
            out[gene.strip()] = float(tpm_s)
        except ValueError as e:
            raise argparse.ArgumentTypeError(f"--cta TPM '{tpm_s}' is not a number") from e
    return out


#: --output extension -> output format. Anything not listed falls back to
#: csv (the safe interchange default for a file); an explicit --format
#: always wins over the extension.
_FORMAT_BY_SUFFIX: dict[str, str] = {
    ".csv": "csv",
    ".tsv": "tsv",
    ".tab": "tsv",
    ".txt": "table",
}


def _resolve_format(explicit: str | None, output: str | None) -> str:
    """Pick the output format: explicit --format, else the --output
    extension, else table for a terminal / csv for an unrecognized file."""
    if explicit:
        return explicit
    if not output:
        return "table"
    return _FORMAT_BY_SUFFIX.get(Path(output).suffix.lower(), "csv")


def _parse_hla(raw: list[str]) -> list[str]:
    """Parse ``--hla`` tokens into canonical allele strings.

    Each token is normalized through mhcgnomes (:func:`tsarina.mhc.parse_mhc`),
    so ``HLA-A*02:01``, ``HLA-A02:01``, and ``A0201`` all resolve to the same
    canonical ``HLA-A*02:01`` -- the ``*`` is a shell glob character, so
    accepting a form without it means ``--hla`` never strictly requires
    quoting."""
    from .mhc import parse_mhc

    out: list[str] = []
    for token in _flatten_multi(raw):
        parsed = parse_mhc(token, expect="allele")
        if parsed is None:
            raise argparse.ArgumentTypeError(f"--hla entry '{token}' is not a recognized allele")
        out.append(parsed.to_string())
    return out


def build_parser(sub: argparse._SubParsersAction) -> argparse.ArgumentParser:
    p = sub.add_parser(
        "personalize",
        help="Build a personalized pMHC target list for a patient.",
        description=(
            "Build a prioritized peptide-MHC target list from a patient's HLA "
            "type, CTA expression, detected mutations, and viral status. "
            "IEDB/CEDAR paths auto-resolve from the hitlist data registry."
        ),
    )
    p.add_argument(
        "--hla",
        required=True,
        nargs="+",
        help=(
            "HLA alleles: comma- and/or space-separated, quoted or not "
            "(e.g. HLA-A*02:01,HLA-B*07:02 or HLA-A*02:01 HLA-B*07:02). "
            "The '*' is optional -- HLA-A02:01 and A0201 both work, which "
            "lets you skip quoting the whole thing (a bare '*' is a shell "
            "glob character)."
        ),
    )
    p.add_argument(
        "--cta",
        nargs="+",
        default=[],
        help=(
            "GENE or GENE=TPM entries: comma- and/or space-separated, "
            "quoted or not (e.g. MAGEA4=142.5,PRAME=87.3 or "
            "MAGEA4=142.5 PRAME=87.3). '=TPM' is optional -- a bare gene "
            "name is included regardless of --min-cta-tpm."
        ),
    )
    p.add_argument(
        "--mutations",
        type=_split_csv,
        default=[],
        help="Comma-separated hotspot labels (e.g. 'KRAS G12D,TP53 R175H').",
    )
    p.add_argument(
        "--viruses",
        type=_split_csv,
        default=[],
        help="Comma-separated virus keys (e.g. hpv16,ebv).",
    )
    p.add_argument(
        "--lengths",
        type=_parse_lengths,
        default=(8, 9, 10, 11),
        help="Peptide lengths (default 8,9,10,11).",
    )
    p.add_argument(
        "--ensembl-release",
        type=int,
        default=112,
        help="Ensembl release (default 112).",
    )
    p.add_argument(
        "--mhc-class",
        choices=("I", "II"),
        default="I",
        help="MHC class filter for IEDB scanning (default I).",
    )
    p.add_argument(
        "--min-cta-tpm",
        type=float,
        default=2.0,
        help="Minimum CTA expression in TPM to include. Defaults to 2.0.",
    )
    p.add_argument(
        "--min-restriction-confidence",
        type=_split_csv,
        default=["HIGH", "MODERATE"],
        help=(
            "Allowed CTA restriction_confidence bins, comma-separated. "
            "Defaults to 'HIGH,MODERATE'; pass 'ANY' to disable."
        ),
    )
    p.add_argument(
        "--mtec-matrix-path",
        dest="mtec_matrix_path",
        default=None,
        help="Path to mTEC gene TPM matrix (TSV); gates CTAs by thymic expression.",
    )
    p.add_argument(
        "--mtec-max-tpm",
        type=float,
        default=1.0,
        help="Maximum mean mTEC TPM when --mtec-matrix-path is given. Defaults to 1.0.",
    )
    p.add_argument(
        "--no-require-human-exclusive-viral",
        dest="require_human_exclusive_viral",
        action="store_false",
        help="Include viral k-mers even if they match a human protein.",
    )
    p.add_argument(
        "--no-enforce-tumor-specificity",
        dest="enforce_tumor_specificity",
        action="store_false",
        help="Keep peptides observed on healthy tissue (default: drop them).",
    )
    p.add_argument(
        "--keep-weak-tier",
        dest="drop_weak_tier",
        action="store_false",
        help="Retain tier-4 (WEAK/unscored) rows in the output.",
    )
    p.add_argument(
        "--no-score",
        action="store_true",
        help="Disable MHC presentation scoring (skip topiary).",
    )
    add_predictor_arg(p, context="presentation scoring")
    add_iedb_cedar_args(p)
    p.add_argument(
        "--skip-ms-evidence",
        action="store_true",
        help="Do not look up IEDB/CEDAR evidence (useful when no data registered).",
    )
    p.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress stage-progress messages on stderr (shown by default).",
    )
    p.add_argument(
        "--no-proteoform-rollup",
        dest="proteoform_rollup",
        action="store_false",
        help=(
            "Report one row per gene symbol instead of collapsing CTAs "
            "that translate to a byte-identical protein into one group "
            "(NY-ESO-1's CTAG1A+CTAG1B, XAGE1A+XAGE1B, SSX2+SSX2B, ...). "
            "Rollup is on by default."
        ),
    )
    p.add_argument(
        "--format",
        choices=("csv", "tsv", "table"),
        default=None,
        help=(
            "Output shape. Inferred when omitted: from --output's "
            "extension (.csv/.tsv/.txt), else 'csv' for any other "
            "--output path, else 'table' when printing to the terminal. "
            "Pass explicitly to override."
        ),
    )
    p.add_argument(
        "-o",
        "--output",
        default=None,
        help=(
            "Write to this path (default: stdout). The extension picks "
            "the format unless --format says otherwise: .csv, .tsv, .txt."
        ),
    )
    return p


def handle(args: argparse.Namespace) -> None:
    from .datasources import DatasetNotRegisteredError
    from .personalize import personalized_targets

    min_restriction_confidence: tuple[str, ...] | None
    if any(v.upper() == "ANY" for v in args.min_restriction_confidence):
        min_restriction_confidence = None
    else:
        min_restriction_confidence = tuple(v.upper() for v in args.min_restriction_confidence)

    try:
        hla_alleles = _parse_hla(args.hla)
        cta_expression = _parse_cta(args.cta)
    except argparse.ArgumentTypeError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        df = personalized_targets(
            hla_alleles=hla_alleles,
            cta_expression=cta_expression or None,
            mutations=args.mutations or None,
            viruses=args.viruses or None,
            lengths=args.lengths,
            ensembl_release=args.ensembl_release,
            iedb_path=args.iedb_path,
            cedar_path=args.cedar_path,
            mhc_class=args.mhc_class,
            min_cta_tpm=args.min_cta_tpm,
            min_restriction_confidence=min_restriction_confidence,
            mtec_matrix_path=args.mtec_matrix_path,
            mtec_max_tpm=args.mtec_max_tpm,
            require_human_exclusive_viral=args.require_human_exclusive_viral,
            enforce_tumor_specificity=args.enforce_tumor_specificity,
            score_presentation=not args.no_score,
            skip_ms_evidence=args.skip_ms_evidence,
            predictor=args.predictor,
            drop_weak_tier=args.drop_weak_tier,
            show_progress=not args.quiet,
            proteoform_rollup=args.proteoform_rollup,
        )
    except DatasetNotRegisteredError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    fmt = _resolve_format(args.format, args.output)

    if fmt == "table":
        from .personalize import format_table

        rendered = format_table(df) + "\n"
    elif fmt == "tsv":
        rendered = df.to_csv(index=False, sep="\t")
    else:
        rendered = df.to_csv(index=False)

    if args.output:
        with open(args.output, "w") as f:
            f.write(rendered)
        print(f"Wrote {len(df)} rows to {args.output}", file=sys.stderr)
    else:
        sys.stdout.write(rendered)
