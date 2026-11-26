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


def _parse_fusion_specs(raw: list[str]) -> list[dict[str, object]]:
    """Parse inline ``--fusion`` specs of form ``LABEL:NTERM|CTERM`` or ``NTERM|CTERM``."""
    records: list[dict[str, object]] = []
    for token in _flatten_multi(raw):
        if "|" not in token:
            raise argparse.ArgumentTypeError(
                f"--fusion entry '{token}' must contain '|' separating N-terminal and C-terminal "
                "sequences, e.g. EWSR1-FLI1:SYGQQSSY|NPSYDSVR"
            )
        if ":" in token:
            label, seq_part = token.split(":", 1)
            label = label.strip()
        else:
            label, seq_part = "", token
        n_seq, c_seq = seq_part.split("|", 1)
        n_seq, c_seq = n_seq.strip().upper(), c_seq.strip().upper()
        if not n_seq or not c_seq:
            raise argparse.ArgumentTypeError(
                f"--fusion entry '{token}' must have non-empty sequences on both sides of '|'"
            )
        rec: dict[str, object] = {
            "n_terminal_seq": n_seq,
            "c_terminal_seq": c_seq,
        }
        if label:
            rec["label"] = label
            if "-" in label:
                g5, g3 = label.split("-", 1)
                if g5.strip() and g3.strip():
                    rec["gene_5prime"] = g5.strip()
                    rec["gene_3prime"] = g3.strip()
        records.append(rec)
    return records


def _load_fusions_file(path_str: str) -> list[dict[str, object]]:
    """Load fusion breakpoint records from a CSV, TSV, or JSON file."""
    import json

    import pandas as pd

    p = Path(path_str)
    if not p.exists():
        raise argparse.ArgumentTypeError(f"--fusions file not found: {p}")
    suffix = p.suffix.lower()
    if suffix == ".json":
        data = json.loads(p.read_text())
        if isinstance(data, dict):
            return [data]
        if isinstance(data, list):
            return [dict(item) for item in data]
        raise argparse.ArgumentTypeError(
            f"--fusions JSON file '{p}' must contain an object or list of objects"
        )
    sep = "\t" if suffix in (".tsv", ".tab") else ","
    df = pd.read_csv(p, sep=sep)
    return [dict(row) for row in df.to_dict(orient="records")]


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
        help="Build a personalized pMHC target list for a patient",
        description=(
            "Build a prioritized peptide-MHC target list from a patient's HLA "
            "type, CTA expression, hotspot mutations, cancer fusions, and viral "
            "status. IEDB/CEDAR paths auto-resolve from the hitlist data registry."
        ),
    )

    inputs = p.add_argument_group("Target inputs")
    inputs.add_argument(
        "--hla",
        "--alleles",
        dest="hla",
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
    inputs.add_argument(
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
    inputs.add_argument(
        "--mutations",
        type=_split_csv,
        default=[],
        help="Comma-separated hotspot labels (e.g. 'KRAS G12D,TP53 R175H').",
    )
    inputs.add_argument(
        "--viruses",
        type=_split_csv,
        default=[],
        help="Comma-separated virus keys (e.g. hpv16,ebv).",
    )
    inputs.add_argument(
        "--fusion",
        nargs="+",
        default=[],
        help=(
            "Inline fusion breakpoint specs of form LABEL:NTERM|CTERM or NTERM|CTERM, "
            "e.g. EWSR1-FLI1:SYGQQSSY|NPSYDSVR."
        ),
    )
    inputs.add_argument(
        "--fusions",
        default=None,
        help=(
            "Path to a CSV, TSV, or JSON file of fusion breakpoint records with "
            "n_terminal_seq/c_terminal_seq or sequence/junction_pos columns."
        ),
    )
    inputs.add_argument(
        "--lengths",
        type=_parse_lengths,
        default=(8, 9, 10, 11),
        help="Peptide lengths (default 8,9,10,11).",
    )
    inputs.add_argument(
        "--ensembl-release",
        type=int,
        default=None,
        help="Ensembl release. Defaults to 112 for human, 115 for mouse, and 98 for dog.",
    )

    safety = p.add_argument_group("Species, sex, and safety gates")
    safety.add_argument(
        "--species",
        default="human",
        help=(
            "Host species common or scientific name resolved via PyEnsembl "
            "(default: human / Homo sapiens; e.g. mouse, dog)."
        ),
    )
    safety.add_argument(
        "--mhc-species",
        default=None,
        help=(
            "Restricting MHC species when distinct from --species, such as "
            "'Homo sapiens' for HLA-transgenic mouse hosts."
        ),
    )
    safety.add_argument(
        "--sex",
        choices=("female", "male"),
        default=None,
        help=(
            "Patient biological sex for OncoRef sex-restricted antigen eligibility: "
            "filters CTAs to sex-eligible targets and admits sex-safe reproductive "
            "antigens such as male-restricted or placental genes in female patients."
        ),
    )
    safety.add_argument(
        "--sex-eligibility",
        choices=("all_adult", "sex_organ_only"),
        default="all_adult",
        help=(
            "Sex-restriction safety policy when --sex is set: 'all_adult' (default, "
            "zero adult normal-tissue expression in that sex) or 'sex_organ_only'."
        ),
    )
    safety.add_argument(
        "--min-cta-tpm",
        type=float,
        default=2.0,
        help="Minimum CTA expression in TPM to include. Defaults to 2.0.",
    )
    safety.add_argument(
        "--min-restriction-confidence",
        type=_split_csv,
        default=["ANY"],
        help=(
            "Allowed CTA restriction_confidence bins, comma-separated. "
            "Default 'ANY' uses the current OncoRef policy without a legacy HPA confidence gate."
        ),
    )
    safety.add_argument(
        "--mtec-matrix-path",
        dest="mtec_matrix_path",
        default=None,
        help="Path to mTEC gene TPM matrix (TSV); gates CTAs by thymic expression.",
    )
    safety.add_argument(
        "--mtec-max-tpm",
        type=float,
        default=1.0,
        help="Maximum mean mTEC TPM when --mtec-matrix-path is given. Defaults to 1.0.",
    )
    safety.add_argument(
        "--no-require-human-exclusive-viral",
        dest="require_human_exclusive_viral",
        action="store_false",
        help="Include viral k-mers even if they match a human protein.",
    )
    safety.add_argument(
        "--no-require-host-exclusive-fusions",
        dest="require_host_exclusive_fusions",
        action="store_false",
        help="Include fusion junction k-mers even if they match a reference host protein.",
    )
    safety.add_argument(
        "--no-enforce-tumor-specificity",
        dest="enforce_tumor_specificity",
        action="store_false",
        help="Keep peptides observed on healthy tissue; default drops them.",
    )
    safety.add_argument(
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

    evidence = p.add_argument_group("MS evidence and presentation scoring")
    evidence.add_argument(
        "--mhc-class",
        choices=("I", "II"),
        default="I",
        help="MHC class filter for IEDB scanning (default I).",
    )
    evidence.add_argument(
        "--ms-evidence-scope",
        choices=("all", "requested_mhc", "any_mhc"),
        default="all",
        help=(
            "MS attribution filter: 'all' retains all candidates with tier "
            "labels; 'requested_mhc' requires exact/donor-set/predicted "
            "support for requested alleles; 'any_mhc' requires any exact-peptide MS."
        ),
    )
    evidence.add_argument(
        "--no-score",
        action="store_true",
        help="Disable MHC presentation scoring (skip topiary).",
    )
    add_predictor_arg(evidence, context="presentation scoring")
    add_iedb_cedar_args(evidence)
    evidence.add_argument(
        "--skip-ms-evidence",
        action="store_true",
        help="Do not look up IEDB/CEDAR evidence (useful when no data registered).",
    )

    output = p.add_argument_group("Output and formatting")
    output.add_argument(
        "--keep-weak-tier",
        dest="drop_weak_tier",
        action="store_false",
        help="Retain tier-4 (WEAK/unscored) rows in the output.",
    )
    output.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress stage-progress messages on stderr; shown by default.",
    )
    output.add_argument(
        "--top",
        type=int,
        default=None,
        help="Limit output to the top N ranked targets; useful with --format table.",
    )
    output.add_argument(
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
    output.add_argument(
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
        fusion_records: list[dict[str, object]] = []
        if getattr(args, "fusion", None):
            fusion_records.extend(_parse_fusion_specs(args.fusion))
        if getattr(args, "fusions", None):
            fusion_records.extend(_load_fusions_file(args.fusions))
    except argparse.ArgumentTypeError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    extra_kwargs: dict[str, object] = {}
    if fusion_records:
        extra_kwargs["fusions"] = fusion_records
        extra_kwargs["require_host_exclusive_fusions"] = getattr(
            args, "require_host_exclusive_fusions", True
        )

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
            evidence_package=getattr(args, "evidence_package", None),
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
            host_species=getattr(args, "species", "human"),
            mhc_species=getattr(args, "mhc_species", None),
            patient_sex=getattr(args, "sex", None),
            sex_eligibility=getattr(args, "sex_eligibility", "all_adult"),
            ms_evidence_scope=getattr(args, "ms_evidence_scope", "all"),
            **extra_kwargs,
        )
    except (DatasetNotRegisteredError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if args.top is not None:
        if args.top <= 0:
            print("Error: --top must be a positive integer.", file=sys.stderr)
            sys.exit(1)
        df = df.head(args.top)

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
