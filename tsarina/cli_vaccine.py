"""Auditable CTA vaccine construction CLI, aligned with Vaxrank naming."""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

from .alleles import panel_names
from .vaccine_construct import VaccineConfig


def build_parser(sub):
    p = sub.add_parser("vaccine", help="Design a mortality-prioritized CTA vaccine antigen")
    p.add_argument("-o", "--output-dir", required=True)
    p.add_argument("-k", "--top-k", type=int, default=10)
    p.add_argument(
        "--selection-mode",
        choices=["ranked", "supported", "budget"],
        default="ranked",
        help="Ranked counts candidates; supported requires k contributors; budget ignores k and fills a length cap by marginal gain",
    )
    p.add_argument(
        "--exclude-gene-pattern",
        action="append",
        default=[],
        help="Repeatable gene-symbol glob, e.g. 'MAGE*'; does not alter CTA membership",
    )
    p.add_argument(
        "--allow-gene",
        action="append",
        default=[],
        help="Exact gene-symbol exception to exclusion patterns, e.g. MAGEA4",
    )
    p.add_argument(
        "--ms-support-mode",
        choices=["presentation", "sample-affinity"],
        default="presentation",
        help="Presentation tier cutoffs, or affinity to any typed sample allele",
    )
    p.add_argument("--ms-affinity-nm", type=float, default=1000)
    p.add_argument(
        "--normal-ms-policy",
        choices=["audit", "exclude"],
        default="audit",
        help="Audit normal MS overlaps, or subtract 8-mers from healthy primary MS outside the CTA tissue scope",
    )
    p.add_argument(
        "--normal-ms-atlas-dir",
        type=Path,
        help="Primary Atlas peptides/donors/sample_hits TSV(.gz) snapshot",
    )
    p.add_argument(
        "--normal-ms-min-donors",
        type=int,
        default=1,
        help="Distinct verified nonmalignant donors for sequence exclusion",
    )
    p.add_argument(
        "--allow-untyped-ms",
        action="store_true",
        help="In sample-affinity mode, permit untyped MS with predicted panel binding",
    )
    p.add_argument("--cta-definition", choices=["strict", "loose", "both"], default="strict")
    p.add_argument("--panel", choices=panel_names(), default="global54_abc")
    p.add_argument(
        "--alleles", "--hla", help="Comma-separated exact HLA-A/B/C alleles; overrides panel"
    )
    p.add_argument(
        "--cancer-cohorts",
        type=Path,
        help="JSON mapping mortality categories to disjoint cohort lists",
    )
    p.add_argument("--ensembl-release", type=int, default=112)
    p.add_argument("--auto-fetch", action="store_true", help="Allow OncoRef expression downloads")
    p.add_argument("--min-padding", type=int, default=0)
    p.add_argument("--max-padding", type=int, default=10)
    p.add_argument(
        "--padding-step", type=int, default=5, help="Search step; maximum is always included"
    )
    p.add_argument("--lengths", default="8,9,10,11", help="Class-I ligand lengths")
    p.add_argument(
        "--predictor", choices=["mhcflurry", "netmhcpan", "netmhcpan_el"], default="mhcflurry"
    )
    p.add_argument("--junction-affinity-nm", type=float, default=1000)
    p.add_argument(
        "--linkers",
        default="AAY",
        help="Comma-separated candidates; direct joins always considered",
    )
    p.add_argument("--beam-width", type=int, default=4)
    p.add_argument("--optimization-rounds", type=int, default=3)
    p.add_argument("--require-clean-junctions", action="store_true")
    p.add_argument("--vaccine-type", choices=["dna", "rna", "mrna"], default="rna")
    p.add_argument("--include-utrs", action="store_true")
    p.add_argument(
        "--utr-5p", default="HBB", help="Vaxrank name, custom nucleotide sequence, or none"
    )
    p.add_argument("--utr-3p", default="HBB_FI")
    p.add_argument("--poly-a-length", type=int, default=0)
    p.add_argument("--max-length-aa", "--max-construct-length-aa", type=int)
    p.add_argument(
        "--max-length-nt",
        "--max-construct-length-nt",
        type=int,
        help="Total nt including UTRs, stop and polyA",
    )
    p.add_argument("--codon-species", default="h_sapiens")


def handle(args):
    from .vaccine import design_vaccine

    try:
        config = VaccineConfig(
            top_k=args.top_k,
            selection_mode=args.selection_mode,
            exclude_gene_patterns=tuple(args.exclude_gene_pattern),
            allow_genes=tuple(args.allow_gene),
            ms_support_mode=args.ms_support_mode.replace("-", "_"),
            ms_affinity_nm=args.ms_affinity_nm,
            allow_untyped_ms=args.allow_untyped_ms,
            normal_ms_policy=args.normal_ms_policy,
            normal_ms_atlas_dir=str(args.normal_ms_atlas_dir) if args.normal_ms_atlas_dir else None,
            normal_ms_min_donors=args.normal_ms_min_donors,
            definition="strict" if args.cta_definition == "both" else args.cta_definition,
            panel=args.panel,
            alleles=tuple(a.strip() for a in args.alleles.split(",")) if args.alleles else None,
            ensembl_release=args.ensembl_release,
            min_padding=args.min_padding,
            max_padding=args.max_padding,
            padding_step=args.padding_step,
            lengths=tuple(map(int, args.lengths.split(","))),
            predictor=args.predictor,
            junction_affinity_nm=args.junction_affinity_nm,
            linkers=tuple(dict.fromkeys(["", *args.linkers.split(",")])),
            beam_width=args.beam_width,
            optimization_rounds=args.optimization_rounds,
            require_clean_junctions=args.require_clean_junctions,
            vaccine_type="rna" if args.vaccine_type == "mrna" else args.vaccine_type,
            include_utrs=args.include_utrs,
            utr_5p=args.utr_5p,
            utr_3p=args.utr_3p,
            poly_a_length=args.poly_a_length,
            max_length_aa=args.max_length_aa,
            max_length_nt=args.max_length_nt,
            codon_species=args.codon_species,
        )
        config.validate()
        cohorts = json.loads(args.cancer_cohorts.read_text()) if args.cancer_cohorts else None
        definitions = ["strict", "loose"] if args.cta_definition == "both" else [config.definition]
        failed = False
        reports = {}
        for definition in definitions:
            out = (
                Path(args.output_dir) / definition
                if len(definitions) > 1
                else Path(args.output_dir)
            )
            try:
                result = design_vaccine(
                    replace(config, definition=definition),
                    cohorts=cohorts,
                    auto_fetch=args.auto_fetch,
                    output_dir=out,
                    on_progress=lambda message: print(message, file=sys.stderr),
                )
            except (ValueError, KeyError, ImportError, FileNotFoundError) as error:
                failed = True
                print(f"{definition} vaccine design failed: {error}", file=sys.stderr)
                continue
            design = result["design"]
            reports[definition] = out
            print(
                f"{definition}: {design['length_aa']} aa, {design['length_nt']} nt; audit: {out / 'report.md'}"
            )
        if len(reports) > 1:
            from .vaccine_website import render_saved_reports

            render_saved_reports(reports, Path(args.output_dir) / "website")
        if failed:
            sys.exit(1)
    except (ValueError, KeyError, ImportError, FileNotFoundError) as error:
        print(f"Vaccine design failed: {error}", file=sys.stderr)
        sys.exit(1)


def build_report_parser(sub):
    p = sub.add_parser(
        "vaccine-report", help="Render a friendly website from verified saved vaccine reports"
    )
    p.add_argument(
        "--report",
        action="append",
        required=True,
        metavar="NAME=PATH",
        help="Repeat for comparisons, e.g. strict=results/strict loose=results/loose",
    )
    p.add_argument("-o", "--output-dir", required=True)
    p.add_argument(
        "--analysis-date", help="Explicit ISO analysis date for reproducible presentation"
    )


def handle_report(args):
    from .vaccine_website import render_saved_reports

    try:
        reports = {}
        for value in args.report:
            name, path = value.split("=", 1)
            if name in reports:
                raise ValueError(f"Duplicate report name: {name}")
            reports[name] = Path(path)
        index = render_saved_reports(reports, args.output_dir, analysis_date=args.analysis_date)
        print(f"Website: {index}")
    except (ValueError, KeyError, FileNotFoundError) as error:
        print(f"Vaccine report failed: {error}", file=sys.stderr)
        sys.exit(1)
