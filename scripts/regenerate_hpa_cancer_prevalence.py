#!/usr/bin/env python
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Regenerate bundled HPA cancer-expression prevalence tables for CTA candidates.

The HPA RNA cancer sample table is large (~1.36 GB gzipped as of HPA v24/v25).
This script streams it and keeps only rows for the genes listed in the supplied
CTA evidence CSV, so future larger CTA candidate sets only require rerunning the
script with the updated CSV.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
import urllib.request
import zipfile
from collections import defaultdict
from collections.abc import Iterator
from contextlib import ExitStack, contextmanager
from pathlib import Path

import pandas as pd

DEFAULT_RNA_URL = "https://www.proteinatlas.org/download/tsv/rna_cancer_sample.tsv.gz"
DEFAULT_IHC_URL = "https://www.proteinatlas.org/download/tsv/cancer_data.tsv.zip"
DEFAULT_THRESHOLDS = (0.1, 1.0, 2.0, 5.0)
DEFAULT_CTA_CSV = Path("tsarina/data/cancer-testis-antigens.csv")
DEFAULT_OUTPUT_DIR = Path("tsarina/data")

_COHORT_RE = re.compile(r"^(?P<cancer>.+?) \((?P<cohort>TCGA|validation)\)$")


@contextmanager
def _open_binary(path_or_url: str | Path) -> Iterator[object]:
    value = str(path_or_url)
    if value.startswith(("http://", "https://")):
        response = urllib.request.urlopen(value, timeout=120)
        try:
            yield response
        finally:
            response.close()
    else:
        with open(value, "rb") as f:
            yield f


@contextmanager
def _optional_gzip_writer(path: Path | None) -> Iterator[object | None]:
    if path is None:
        yield None
        return
    with gzip.open(path, "wb") as f:
        yield f


def _threshold_label(value: float) -> str:
    text = f"{value:g}".replace(".", "_")
    return text.replace("-", "minus_")


def _split_cancer_label(value: str) -> tuple[str, str]:
    match = _COHORT_RE.match(value)
    if not match:
        return value, ""
    return match.group("cancer"), match.group("cohort")


def _load_cta_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"Symbol", "Ensembl_Gene_ID"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required column(s): {sorted(missing)}")

    out = df.copy()
    out["gene_id"] = out["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
    out["symbol"] = out["Symbol"].astype(str)
    return out.drop_duplicates(subset=["gene_id"])


def _metadata_columns(cta: pd.DataFrame) -> list[str]:
    wanted = [
        "gene_id",
        "symbol",
        "source_databases",
        "passes_filters",
        "filtered",
        "never_expressed",
        "restriction",
        "restriction_confidence",
    ]
    return [column for column in wanted if column in cta.columns]


def regenerate_rna_prevalence(
    *,
    cta: pd.DataFrame,
    source: str | Path,
    output_path: Path,
    thresholds: tuple[float, ...] = DEFAULT_THRESHOLDS,
    sample_subset_output: Path | None = None,
    progress_every: int = 20_000_000,
) -> pd.DataFrame:
    """Stream HPA RNA cancer samples and write per-gene/cancer prevalence."""
    gene_ids = set(cta["gene_id"].astype(str))
    gene_id_bytes = {gene_id.encode("utf-8") for gene_id in gene_ids}
    meta_cols = [column for column in _metadata_columns(cta) if column != "gene_id"]
    meta = cta.set_index("gene_id")[meta_cols].to_dict(orient="index")

    stats = defaultdict(lambda: [0, [0] * len(thresholds), 0.0, 0.0])
    scanned = 0
    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with ExitStack() as stack:
        subset_fh = stack.enter_context(_optional_gzip_writer(sample_subset_output))
        source_fh = stack.enter_context(_open_binary(source))
        gz = stack.enter_context(gzip.GzipFile(fileobj=source_fh))
        header = gz.readline()
        if subset_fh is not None:
            subset_fh.write(header)
        for line in gz:
            scanned += 1
            gene_id = line.split(b"\t", 1)[0]
            if gene_id not in gene_id_bytes:
                if progress_every and scanned % progress_every == 0:
                    print(f"scanned={scanned:,} kept={kept:,}", file=sys.stderr)
                continue

            kept += 1
            if subset_fh is not None:
                subset_fh.write(line)

            parts = line.rstrip(b"\n").split(b"\t")
            if len(parts) != 4:
                continue
            gene_text = parts[0].decode("utf-8")
            cancer_label = parts[2].decode("utf-8")
            try:
                ptpm = float(parts[3])
            except ValueError:
                continue
            cancer_type, cohort = _split_cancer_label(cancer_label)
            rec = stats[(gene_text, cancer_label, cancer_type, cohort)]
            rec[0] += 1
            for i, threshold in enumerate(thresholds):
                rec[1][i] += int(ptpm >= threshold)
            rec[2] += ptpm
            rec[3] = max(rec[3], ptpm)

            if progress_every and scanned % progress_every == 0:
                print(f"scanned={scanned:,} kept={kept:,}", file=sys.stderr)

    records: list[dict[str, object]] = []
    for (gene_id, cancer_label, cancer_type, cohort), (n, counts, total, max_ptpm) in sorted(
        stats.items()
    ):
        row: dict[str, object] = {
            "gene_id": gene_id,
            "symbol": str(meta.get(gene_id, {}).get("symbol", "")),
            "cancer": cancer_label,
            "cancer_type": cancer_type,
            "cohort": cohort,
            "samples": n,
            "mean_ptpm": total / n if n else 0.0,
            "max_ptpm": max_ptpm,
        }
        for threshold, count in zip(thresholds, counts):
            label = _threshold_label(threshold)
            row[f"expressed_samples_ptpm_ge_{label}"] = count
            row[f"prevalence_ptpm_ge_{label}"] = count / n if n else 0.0
        for column, value in meta.get(gene_id, {}).items():
            if column not in row:
                row[column] = value
        records.append(row)

    out = pd.DataFrame(records)
    out.to_csv(output_path, index=False)
    matched = set(out["gene_id"]) if not out.empty else set()
    missing = sorted(gene_ids - matched)
    if missing:
        print(
            f"Warning: {len(missing)} CTA gene(s) had no HPA RNA cancer rows.",
            file=sys.stderr,
        )
    print(
        f"RNA: scanned={scanned:,} kept={kept:,} genes={len(matched):,} rows={len(out):,}",
        file=sys.stderr,
    )
    return out


def regenerate_ihc_prevalence(
    *,
    cta: pd.DataFrame,
    source: str | Path,
    output_path: Path,
) -> pd.DataFrame:
    """Subset HPA cancer IHC count table to CTA candidates."""
    gene_ids = set(cta["gene_id"].astype(str))
    meta_cols = _metadata_columns(cta)
    meta = cta[meta_cols].drop_duplicates(subset=["gene_id"])

    with _open_binary(source) as source_fh, zipfile.ZipFile(source_fh) as zf:
        names = zf.namelist()
        if not names:
            raise ValueError(f"{source} is empty")
        ihc = pd.read_csv(zf.open(names[0]), sep="\t")

    subset = ihc[ihc["Gene"].astype(str).isin(gene_ids)].copy()
    subset = subset.merge(meta, left_on="Gene", right_on="gene_id", how="left")
    subset = subset.drop(columns=["gene_id"])
    subset = subset.rename(
        columns={
            "Gene": "gene_id",
            "Gene name": "hpa_gene_name",
            "Cancer": "cancer",
            "High": "high",
            "Medium": "medium",
            "Low": "low",
            "Not detected": "not_detected",
        }
    )
    count_cols = ["high", "medium", "low", "not_detected"]
    subset["total"] = subset[count_cols].sum(axis=1)
    subset["detected"] = subset[["high", "medium", "low"]].sum(axis=1)
    subset["medium_high"] = subset[["high", "medium"]].sum(axis=1)
    denom = subset["total"].where(subset["total"] != 0, pd.NA)
    subset["prevalence_detected"] = subset["detected"] / denom
    subset["prevalence_medium_high"] = subset["medium_high"] / denom

    columns = [
        "gene_id",
        "symbol",
        "hpa_gene_name",
        "cancer",
        "high",
        "medium",
        "low",
        "not_detected",
        "total",
        "detected",
        "medium_high",
        "prevalence_detected",
        "prevalence_medium_high",
        *[
            column
            for column in meta_cols
            if column not in {"gene_id", "symbol"} and column in subset.columns
        ],
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    subset[columns].to_csv(output_path, index=False)
    print(
        "IHC: "
        f"source_rows={len(ihc):,} kept_rows={len(subset):,} genes={subset['gene_id'].nunique():,}",
        file=sys.stderr,
    )
    return subset[columns]


def _parse_thresholds(value: str) -> tuple[float, ...]:
    thresholds = tuple(float(part.strip()) for part in value.split(",") if part.strip())
    if not thresholds:
        raise argparse.ArgumentTypeError("at least one threshold is required")
    return thresholds


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cta-csv", type=Path, default=DEFAULT_CTA_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--rna-source", default=DEFAULT_RNA_URL)
    parser.add_argument("--ihc-source", default=DEFAULT_IHC_URL)
    parser.add_argument("--thresholds", type=_parse_thresholds, default=DEFAULT_THRESHOLDS)
    parser.add_argument(
        "--rna-output",
        type=Path,
        default=None,
        help="Default: <output-dir>/hpa-cancer-rna-prevalence.csv",
    )
    parser.add_argument(
        "--ihc-output",
        type=Path,
        default=None,
        help="Default: <output-dir>/hpa-cancer-ihc-prevalence.csv",
    )
    parser.add_argument(
        "--rna-sample-subset-output",
        type=Path,
        default=None,
        help="Optional gzip path for the CTA-only sample-level RNA rows.",
    )
    parser.add_argument("--skip-rna", action="store_true")
    parser.add_argument("--skip-ihc", action="store_true")
    parser.add_argument("--progress-every", type=int, default=20_000_000)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    cta = _load_cta_metadata(args.cta_csv)
    rna_output = args.rna_output or args.output_dir / "hpa-cancer-rna-prevalence.csv"
    ihc_output = args.ihc_output or args.output_dir / "hpa-cancer-ihc-prevalence.csv"

    if not args.skip_rna:
        regenerate_rna_prevalence(
            cta=cta,
            source=args.rna_source,
            output_path=rna_output,
            thresholds=args.thresholds,
            sample_subset_output=args.rna_sample_subset_output,
            progress_every=args.progress_every,
        )
    if not args.skip_ihc:
        regenerate_ihc_prevalence(cta=cta, source=args.ihc_source, output_path=ihc_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
