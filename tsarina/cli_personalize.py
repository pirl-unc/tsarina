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

Wraps :func:`tsarina.personalize.personalize` with argparse.  IEDB/CEDAR paths
default to whatever the user registered with ``tsarina data register …``.
"""

from __future__ import annotations

import argparse
import sys

_SUPPORTED_PREDICTORS = ("mhcflurry", "netmhcpan", "netmhcpan_el")


def _split_csv(value: str) -> list[str]:
    return [s.strip() for s in value.split(",") if s.strip()]


def _parse_cta(value: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for pair in _split_csv(value):
        if "=" not in pair:
            raise argparse.ArgumentTypeError(
                f"--cta entry '{pair}' must be GENE=TPM (e.g. PRAME=87.3)"
            )
        gene, _, tpm_s = pair.partition("=")
        try:
            out[gene.strip()] = float(tpm_s)
        except ValueError as e:
            raise argparse.ArgumentTypeError(f"--cta TPM '{tpm_s}' is not a number") from e
    return out


def _parse_lengths(value: str) -> tuple[int, ...]:
    try:
        return tuple(int(x) for x in _split_csv(value))
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"--lengths must be integers: {value!r}") from e


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
        type=_split_csv,
        help="Comma-separated HLA alleles (e.g. HLA-A*02:01,HLA-B*07:02).",
    )
    p.add_argument(
        "--cta",
        type=_parse_cta,
        default={},
        help="Comma-separated GENE=TPM pairs (e.g. MAGEA4=142.5,PRAME=87.3).",
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
        help="Minimum CTA expression in TPM to include (default 2.0).",
    )
    p.add_argument(
        "--min-restriction-confidence",
        type=_split_csv,
        default=["HIGH", "MODERATE"],
        help=(
            "Allowed CTA restriction_confidence bins, comma-separated "
            "(default 'HIGH,MODERATE'; pass 'ANY' to disable)."
        ),
    )
    p.add_argument(
        "--mtec-matrix",
        dest="mtec_matrix_path",
        default=None,
        help="Path to mTEC gene TPM matrix (TSV); gates CTAs by thymic expression.",
    )
    p.add_argument(
        "--mtec-max-tpm",
        type=float,
        default=1.0,
        help="Maximum mean mTEC TPM when --mtec-matrix is given (default 1.0).",
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
    p.add_argument(
        "--predictor",
        choices=_SUPPORTED_PREDICTORS,
        default="mhcflurry",
        help="mhctools predictor used for presentation scoring (default mhcflurry).",
    )
    p.add_argument(
        "--iedb",
        dest="iedb_path",
        default=None,
        help="Override IEDB path (default: hitlist registry).",
    )
    p.add_argument(
        "--cedar",
        dest="cedar_path",
        default=None,
        help="Override CEDAR path (default: hitlist registry).",
    )
    p.add_argument(
        "--skip-ms-evidence",
        action="store_true",
        help="Do not look up IEDB/CEDAR evidence (useful when no data registered).",
    )
    p.add_argument(
        "-o",
        "--output",
        default=None,
        help="Write CSV to this path (default: stdout).",
    )
    return p


def handle(args: argparse.Namespace) -> None:
    from .datasources import DatasetNotRegisteredError
    from .personalize import personalize

    min_restriction_confidence: tuple[str, ...] | None
    if any(v.upper() == "ANY" for v in args.min_restriction_confidence):
        min_restriction_confidence = None
    else:
        min_restriction_confidence = tuple(v.upper() for v in args.min_restriction_confidence)

    try:
        df = personalize(
            hla_alleles=args.hla,
            cta_expression=args.cta or None,
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
        )
    except DatasetNotRegisteredError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if args.output:
        df.to_csv(args.output, index=False)
        print(f"Wrote {len(df)} rows to {args.output}", file=sys.stderr)
    else:
        df.to_csv(sys.stdout, index=False)
