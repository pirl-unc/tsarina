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

"""`tsarina hits` CLI — all public MS hits for a protein, with filtering.

Enumerates k-mers from a gene's canonical protein (via
``hitlist.proteome.ProteomeIndex.from_ensembl``), scans IEDB/CEDAR for those
peptides, and applies allele / species / serotype filters.  Output can be
aggregated per-peptide, per-(peptide, allele), or returned as raw scan rows.
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd

_SUPPORTED_FORMATS = ("peptides", "pmhc", "raw")
_SUPPORTED_PREDICTORS = ("mhcflurry", "netmhcpan", "netmhcpan_el")
_SUPPORTED_RESOLUTIONS = ("four_digit", "two_digit", "serological", "class_only")


def _split_csv(value: str) -> list[str]:
    return [s.strip() for s in value.split(",") if s.strip()]


def _parse_lengths(value: str) -> tuple[int, ...]:
    try:
        return tuple(int(x) for x in _split_csv(value))
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"--lengths must be integers: {value!r}") from e


def build_parser(sub: argparse._SubParsersAction) -> argparse.ArgumentParser:
    p = sub.add_parser(
        "hits",
        help="List all public MS hits for a protein, with filtering.",
        description=(
            "Enumerate k-mers from a gene's canonical protein, scan IEDB/CEDAR "
            "for matches, and filter by allele / species / serotype / "
            "resolution.  IEDB/CEDAR paths auto-resolve from the hitlist "
            "data registry."
        ),
    )
    target = p.add_mutually_exclusive_group(required=True)
    target.add_argument("--gene", help="HGNC gene symbol (e.g. PRAME).")
    target.add_argument(
        "--uniprot",
        help="UniProt accession. Best-effort — maps via gene symbol cross-reference.",
    )

    p.add_argument(
        "--allele",
        type=_split_csv,
        default=[],
        help="Comma-separated MHC restrictions to keep (e.g. 'HLA-A*24:02,HLA-A*02:01').",
    )
    p.add_argument(
        "--species",
        default="Homo sapiens",
        help="MHC species filter (default 'Homo sapiens'). Use 'any' to disable.",
    )
    p.add_argument(
        "--serotype",
        type=_split_csv,
        default=[],
        help="Comma-separated serotype labels to keep (e.g. A2,A24).",
    )
    p.add_argument(
        "--min-resolution",
        choices=_SUPPORTED_RESOLUTIONS,
        default=None,
        help="Minimum allele resolution to keep.",
    )
    p.add_argument(
        "--mhc-class",
        choices=("I", "II"),
        default=None,
        help="MHC class filter.",
    )
    p.add_argument(
        "--lengths",
        type=_parse_lengths,
        default=(8, 9, 10, 11),
        help="Peptide lengths to enumerate (default 8,9,10,11).",
    )
    p.add_argument(
        "--ensembl-release",
        type=int,
        default=112,
        help="Ensembl release (default 112).",
    )
    p.add_argument(
        "--include-binding-assays",
        action="store_true",
        help="Keep IEDB binding-assay rows (default: MS only).",
    )
    p.add_argument(
        "--format",
        choices=_SUPPORTED_FORMATS,
        default="pmhc",
        help="Output aggregation (default pmhc).",
    )
    p.add_argument(
        "--predict",
        action="store_true",
        help="Also score observed (peptide, allele) pairs via topiary.",
    )
    p.add_argument(
        "--predictor",
        choices=_SUPPORTED_PREDICTORS,
        default="mhcflurry",
        help="mhctools predictor for --predict (default mhcflurry).",
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
        help="Skip IEDB/CEDAR; output the proteome-enumerated peptides only.",
    )
    p.add_argument(
        "-o",
        "--output",
        default=None,
        help="Write CSV to this path (default: stdout).",
    )
    return p


def _enumerate_gene_peptides(
    gene: str,
    ensembl_release: int,
    lengths: tuple[int, ...],
) -> pd.DataFrame:
    """Build a DataFrame of (peptide, protein_id, gene_name, gene_id, position, n_flank, c_flank)."""
    from hitlist.proteome import ProteomeIndex

    idx = ProteomeIndex.from_ensembl(release=ensembl_release, lengths=lengths, verbose=False)

    target_ids = [pid for pid, meta in idx.protein_meta.items() if meta.get("gene_name") == gene]
    if not target_ids:
        raise ValueError(
            f"No proteins found for gene '{gene}' in Ensembl release {ensembl_release}."
        )

    rows: list[dict] = []
    for pid in target_ids:
        seq = idx.proteins[pid]
        meta = idx.protein_meta[pid]
        for L in lengths:
            for i in range(len(seq) - L + 1):
                peptide = seq[i : i + L]
                rows.append(
                    {
                        "peptide": peptide,
                        "protein_id": pid,
                        "gene_name": meta.get("gene_name", ""),
                        "gene_id": meta.get("gene_id", ""),
                        "position": i,
                        "n_flank": seq[max(0, i - 5) : i],
                        "c_flank": seq[i + L : i + L + 5],
                    }
                )
    return pd.DataFrame(rows)


def _resolve_uniprot_to_gene(uniprot: str) -> str:
    """Best-effort UniProt → gene symbol via UniProt REST API."""
    import urllib.error
    import urllib.request

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot}.json"
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            import json

            data = json.loads(resp.read())
    except (urllib.error.URLError, TimeoutError) as e:
        raise ValueError(f"Could not resolve UniProt accession '{uniprot}': {e}") from e

    for entry in data.get("genes", []):
        name = entry.get("geneName", {}).get("value")
        if name:
            return name
    raise ValueError(f"No gene name found for UniProt '{uniprot}'.")


def _filter_by_allele(hits: pd.DataFrame, alleles: list[str]) -> pd.DataFrame:
    if not alleles:
        return hits
    wanted = {a.strip() for a in alleles}
    return hits[hits["mhc_restriction"].isin(wanted)].copy()


def _filter_by_serotype(hits: pd.DataFrame, serotypes: list[str]) -> pd.DataFrame:
    if not serotypes:
        return hits
    from hitlist.curation import allele_to_serotype

    wanted = {s.strip() for s in serotypes}
    mask = hits["mhc_restriction"].map(lambda r: allele_to_serotype(r) in wanted)
    return hits[mask].copy()


def handle(args: argparse.Namespace) -> None:
    # ── 1. Gene resolution ─────────────────────────────────────────────
    if args.gene:
        gene = args.gene
    else:
        try:
            gene = _resolve_uniprot_to_gene(args.uniprot)
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        print(f"UniProt {args.uniprot} → gene {gene}", file=sys.stderr)

    # ── 2. Enumerate peptides ──────────────────────────────────────────
    try:
        pep_df = _enumerate_gene_peptides(gene, args.ensembl_release, args.lengths)
    except (ValueError, ImportError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if args.skip_ms_evidence:
        out = pep_df
        _write(args.output, out)
        return

    # ── 3. Load MS evidence ────────────────────────────────────────────
    from hitlist.aggregate import aggregate_per_peptide, aggregate_per_pmhc

    mhc_species = None if args.species.lower() == "any" else args.species
    peptide_set = set(pep_df["peptide"].unique())

    if args.iedb_path is not None or args.cedar_path is not None:
        # Explicit raw-CSV override path.
        from hitlist.scanner import scan

        from .datasources import DatasetNotRegisteredError, resolve_dataset_paths

        try:
            iedb_path, cedar_path = resolve_dataset_paths(args.iedb_path, args.cedar_path)
        except DatasetNotRegisteredError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        hits = scan(
            peptides=peptide_set,
            iedb_path=str(iedb_path),
            cedar_path=str(cedar_path) if cedar_path else None,
            mhc_class=args.mhc_class,
            mhc_species=mhc_species,
            classify_source=True,
            min_allele_resolution=args.min_resolution,
        )
        if not args.include_binding_assays and "is_binding_assay" in hits.columns:
            hits = hits[~hits["is_binding_assay"]].copy()
    else:
        # Fast path: cached observations index.
        from .indexing import load_ms_evidence

        hits = load_ms_evidence(
            peptides=peptide_set,
            mhc_class=args.mhc_class,
            mhc_species=mhc_species,
            drop_binding_assays=not args.include_binding_assays,
        )
        if args.min_resolution is not None and not hits.empty:
            from hitlist.curation import allele_resolution_rank, classify_allele_resolution

            min_rank = allele_resolution_rank(args.min_resolution)
            mask = hits["mhc_restriction"].map(
                lambda r: allele_resolution_rank(classify_allele_resolution(r)) <= min_rank
            )
            hits = hits[mask].copy()

    hits = _filter_by_allele(hits, args.allele)
    hits = _filter_by_serotype(hits, args.serotype)

    # ── 4. Aggregate per requested format ──────────────────────────────
    if args.format == "peptides":
        out = aggregate_per_peptide(hits)
        gene_cols = pep_df[["peptide", "protein_id", "gene_name", "gene_id"]].drop_duplicates(
            subset="peptide"
        )
        out = gene_cols.merge(out, on="peptide", how="inner")
    elif args.format == "pmhc":
        out = aggregate_per_pmhc(hits)
        gene_cols = pep_df[["peptide", "protein_id", "gene_name", "gene_id"]].drop_duplicates(
            subset="peptide"
        )
        out = gene_cols.merge(out, on="peptide", how="inner")
    else:
        # raw: join gene context per peptide
        gene_cols = pep_df[
            ["peptide", "protein_id", "gene_name", "gene_id", "position", "n_flank", "c_flank"]
        ].drop_duplicates(subset=["peptide", "protein_id", "position"])
        out = hits.merge(gene_cols, on="peptide", how="left")

    # ── 5. Optional topiary scoring ────────────────────────────────────
    if args.predict and not out.empty:
        from .scoring import score_presentation

        if args.format == "pmhc":
            peptide_allele = out[["peptide", "mhc_restriction"]].drop_duplicates()
            scores = score_presentation(
                peptides=peptide_allele["peptide"].tolist(),
                alleles=peptide_allele["mhc_restriction"].unique().tolist(),
                predictor=args.predictor,
            ).rename(columns={"allele": "mhc_restriction"})
            out = out.merge(scores, on=["peptide", "mhc_restriction"], how="left")
        else:
            # peptides / raw: score each peptide against observed alleles collectively.
            alleles_col = "ms_alleles" if "ms_alleles" in out.columns else "mhc_restriction"
            if alleles_col in out.columns:
                alleles = sorted(
                    {
                        a
                        for cell in out[alleles_col].dropna()
                        for a in (cell.split(";") if isinstance(cell, str) else [cell])
                    }
                )
                if alleles:
                    scores = score_presentation(
                        peptides=out["peptide"].tolist(),
                        alleles=alleles,
                        predictor=args.predictor,
                    )
                    best = (
                        scores.sort_values("presentation_percentile")
                        .groupby("peptide", as_index=False)
                        .first()
                        .rename(
                            columns={
                                "allele": "best_allele",
                                "presentation_percentile": "best_presentation_percentile",
                                "presentation_score": "best_presentation_score",
                                "affinity_nm": "best_affinity_nm",
                            }
                        )
                    )
                    out = out.merge(best, on="peptide", how="left")

    _write(args.output, out)


def _write(output: str | None, df: pd.DataFrame) -> None:
    if output:
        df.to_csv(output, index=False)
        print(f"Wrote {len(df)} rows to {output}", file=sys.stderr)
    else:
        df.to_csv(sys.stdout, index=False)
