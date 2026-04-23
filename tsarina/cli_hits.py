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

Default path queries the hitlist observations index with
``load_observations(gene_name=...)`` — multi-mapping peptide → gene
attribution comes straight from the mappings sidecar, no Ensembl build
required.

The slower ProteomeIndex-based enumeration is kept for two niche cases:
``--skip-ms-evidence`` (output the theoretical peptide menu for a gene
from a fresh Ensembl release) and ``--iedb``/``--cedar`` (raw CSV
override that bypasses the index).
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd

_SUPPORTED_FORMATS = ("peptides", "pmhc", "refs", "raw")
_SUPPORTED_PREDICTORS = ("mhcflurry", "netmhcpan", "netmhcpan_el")
_SUPPORTED_RESOLUTIONS = ("four_digit", "two_digit", "serological", "class_only")


def _split_csv(value: str) -> list[str]:
    return [s.strip() for s in value.split(",") if s.strip()]


def _parse_lengths(value: str) -> tuple[int, ...]:
    try:
        return tuple(int(x) for x in _split_csv(value))
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"--lengths must be integers: {value!r}") from e


# Class-default peptide length windows used when --lengths is omitted.
# Intentionally wider than textbook 8-11 / 12-25 to catch non-canonical
# ligands hitlist actually curates (phospho-extended class I, long class
# II tails).  An explicit --lengths overrides these.
_CLASS_I_DEFAULT_LENGTHS = tuple(range(8, 16))  # 8-15
_CLASS_II_DEFAULT_LENGTHS = tuple(range(12, 46))  # 12-45


def _resolve_lengths(args: argparse.Namespace) -> tuple[int, ...] | None:
    """Resolve --lengths with --mhc-class fallbacks.

    - Explicit --lengths always wins.
    - Otherwise derive from --mhc-class: class I -> 8-15, class II -> 12-45.
    - If neither is set, return None so the cached observations path skips
      length bounds (preserves "load everything" semantics for unqualified
      queries rather than silently dropping off-window rows).
    """
    if args.lengths is not None:
        return args.lengths
    if args.mhc_class == "I":
        return _CLASS_I_DEFAULT_LENGTHS
    if args.mhc_class == "II":
        return _CLASS_II_DEFAULT_LENGTHS
    return None


def build_parser(sub: argparse._SubParsersAction) -> argparse.ArgumentParser:
    p = sub.add_parser(
        "hits",
        help="List all public MS hits for a protein, with filtering.",
        description=(
            "List all public MS hits for a gene or protein from the hitlist "
            "observations index, filtered by allele / species / serotype / "
            "resolution / MHC class.  Output includes semicolon-joined "
            "gene_names / gene_ids / protein_ids so paralogous attribution "
            "is preserved (e.g. MAGE family peptides).  Build the index once "
            "with `tsarina data build` — it auto-builds on first use otherwise."
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
        default=None,
        help=(
            "Peptide lengths. On the enumeration paths (--skip-ms-evidence / "
            "--iedb / --cedar) controls the k-mer walk; on the cached "
            "observations path bounds the loader's length filter "
            "(hitlist>=1.15.1). If omitted, defaults are derived from "
            "--mhc-class: 8-15 for class I, 12-45 for class II. When "
            "neither --lengths nor --mhc-class is set, the cached path "
            "skips length bounds entirely so nothing is silently dropped."
        ),
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
        help=(
            "Include IEDB/CEDAR binding-assay rows in addition to MS "
            "observations (default: MS only).  On the cached-index path this "
            "routes through hitlist's load_all_evidence() union; rows are "
            "tagged with an evidence_kind column ('ms' vs 'binding')."
        ),
    )
    p.add_argument(
        "--mono-allelic-only",
        action="store_true",
        help=(
            "Keep only hits with direct mono-allelic evidence (observed on a "
            "cell line expressing a single HLA allele).  Excludes multi-allelic "
            "studies where the allele was assigned by a predictor (NetMHCpan / "
            "MHCflurry) rather than observed directly."
        ),
    )
    p.add_argument(
        "--format",
        choices=_SUPPORTED_FORMATS,
        default="pmhc",
        help=(
            "Output aggregation (default pmhc). 'refs' adds per-pMHC columns "
            "for ms_pmhc_ref_count, ms_pmhc_pmids, ms_pmhc_tissues, "
            "ms_pmhc_diseases, ms_pmhc_cell_lines, and the "
            "ms_pmhc_in_cancer / ms_pmhc_in_healthy_tissue source flags. "
            "Backed by hitlist.aggregate.aggregate_per_pmhc_with_refs."
        ),
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
    """Filter observations to the given serotypes.

    Uses mhcgnomes to expand each requested serotype into its full allele
    list, then checks membership.  Handles three cases:

    1. Molecular restriction (``HLA-A*02:01``) — matches if the parsed allele
       is in the requested serotype's ``alleles`` list.
    2. Serological restriction (``HLA-A2``) — matches if the parsed serotype
       name equals the requested name.
    3. Unparseable / class-only — skipped.

    Workaround for hitlist#44 (``allele_to_serotype`` mis-reports A-locus
    Bw4-carrying alleles): we don't depend on the canonical-serotype label,
    we go straight to mhcgnomes membership.
    """
    if not serotypes:
        return hits
    from mhcgnomes import parse as mhc_parse
    from mhcgnomes.allele import Allele
    from mhcgnomes.serotype import Serotype

    wanted_allele_keys: set[tuple[str, tuple[str, ...]]] = set()
    wanted_serotype_names: set[str] = set()
    for raw in serotypes:
        s = raw.strip()
        if not s:
            continue
        wanted_serotype_names.add(s.removeprefix("HLA-"))
        try:
            parsed = mhc_parse(s if s.startswith("HLA-") else f"HLA-{s}")
        except Exception:
            continue
        if isinstance(parsed, Serotype):
            for a in parsed.alleles:
                wanted_allele_keys.add((a.gene.name, a.allele_fields))

    def _matches(mhc: str) -> bool:
        if not isinstance(mhc, str) or not mhc:
            return False
        try:
            parsed = mhc_parse(mhc)
        except Exception:
            return False
        if isinstance(parsed, Allele):
            return (parsed.gene.name, parsed.allele_fields) in wanted_allele_keys
        if isinstance(parsed, Serotype):
            return parsed.name in wanted_serotype_names
        return False

    return hits[hits["mhc_restriction"].map(_matches)].copy()


def _apply_min_resolution(hits: pd.DataFrame, min_resolution: str | None) -> pd.DataFrame:
    if min_resolution is None or hits.empty:
        return hits
    from hitlist.curation import allele_resolution_rank, classify_allele_resolution

    min_rank = allele_resolution_rank(min_resolution)
    mask = hits["mhc_restriction"].map(
        lambda r: allele_resolution_rank(classify_allele_resolution(r)) <= min_rank
    )
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

    mhc_species = None if args.species.lower() == "any" else args.species

    # ── 2. Proteome-enumeration fallbacks ──────────────────────────────
    # --skip-ms-evidence outputs the theoretical peptide menu.
    # --iedb / --cedar overrides still need a peptide list for raw scan.
    needs_peptide_enumeration = args.skip_ms_evidence or (
        args.iedb_path is not None or args.cedar_path is not None
    )
    # Single resolution pass used by both the cached path (for bounds
    # pushdown) and the enumeration path (for the k-mer walk).  May be
    # None when neither --lengths nor --mhc-class was given; enumeration
    # can't run with None, so fall back to the class-I window there.
    resolved_lengths = _resolve_lengths(args)
    enumeration_lengths = resolved_lengths or _CLASS_I_DEFAULT_LENGTHS
    pep_df: pd.DataFrame | None = None
    if needs_peptide_enumeration:
        try:
            pep_df = _enumerate_gene_peptides(gene, args.ensembl_release, enumeration_lengths)
        except (ValueError, ImportError) as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

    if args.skip_ms_evidence:
        _write(args.output, pep_df)
        return

    # ── 3. Load MS evidence ────────────────────────────────────────────
    from hitlist.aggregate import (
        aggregate_per_peptide,
        aggregate_per_pmhc,
        aggregate_per_pmhc_with_refs,
    )

    if args.iedb_path is not None or args.cedar_path is not None:
        # Explicit raw-CSV override path — scan only, no observations index.
        from hitlist.scanner import scan

        from .datasources import DatasetNotRegisteredError, resolve_dataset_paths

        try:
            iedb_path, cedar_path = resolve_dataset_paths(args.iedb_path, args.cedar_path)
        except DatasetNotRegisteredError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        assert pep_df is not None  # enumerated above
        hits = scan(
            peptides=set(pep_df["peptide"].unique()),
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
        # Fast path: direct gene_name pushdown through hitlist's mappings sidecar.
        from .indexing import ensure_index_built

        ensure_index_built()

        # hitlist 1.15.1+ length_min/length_max pushdown.  Applied when
        # lengths were explicitly given OR a class was specified (in
        # which case _resolve_lengths picked a class-appropriate window).
        # When neither was set resolved_lengths is None and no bounds
        # are pushed, preserving the "load everything" behavior so
        # unqualified queries don't silently drop off-window rows.
        # Non-contiguous --lengths are restored to exact-set semantics
        # by the post-filter below.
        length_min = min(resolved_lengths) if resolved_lengths else None
        length_max = max(resolved_lengths) if resolved_lengths else None

        if args.include_binding_assays:
            # hitlist 1.10.0+ exposes load_all_evidence() — UNION of MS
            # observations + binding-assay rows, tagged with evidence_kind.
            from hitlist.observations import load_all_evidence

            hits = load_all_evidence(
                gene_name=gene,
                mhc_class=args.mhc_class,
                species=mhc_species,
                length_min=length_min,
                length_max=length_max,
            )
        else:
            from hitlist.observations import load_observations

            hits = load_observations(
                gene_name=gene,
                mhc_class=args.mhc_class,
                species=mhc_species,
                length_min=length_min,
                length_max=length_max,
            )
        if resolved_lengths and not hits.empty:
            hits = hits[hits["peptide"].str.len().isin(resolved_lengths)].copy()
        hits = _apply_min_resolution(hits, args.min_resolution)

    hits = _filter_by_allele(hits, args.allele)
    hits = _filter_by_serotype(hits, args.serotype)

    if args.mono_allelic_only and not hits.empty and "is_monoallelic" in hits.columns:
        hits = hits[hits["is_monoallelic"]].copy()

    if hits.empty:
        _write(args.output, pd.DataFrame({"peptide": pd.Series(dtype=str)}))
        return

    # ── 4. Aggregate per requested format ──────────────────────────────
    #
    # The observations index already carries semicolon-joined
    # gene_names / gene_ids / protein_ids.  For the raw-CSV override path,
    # join identifiers from the enumeration frame so output shape stays
    # consistent across paths.
    gene_ident_cols = [c for c in ("gene_names", "gene_ids", "protein_ids") if c in hits.columns]
    if args.format == "peptides":
        out = aggregate_per_peptide(hits)
        if gene_ident_cols:
            ids = hits[["peptide", *gene_ident_cols]].drop_duplicates(subset="peptide")
            out = ids.merge(out, on="peptide", how="inner")
        elif pep_df is not None:
            legacy = pep_df[["peptide", "protein_id", "gene_name", "gene_id"]].drop_duplicates(
                subset="peptide"
            )
            out = legacy.merge(out, on="peptide", how="inner")
    elif args.format == "pmhc":
        out = aggregate_per_pmhc(hits)
        if gene_ident_cols:
            ids = hits[["peptide", *gene_ident_cols]].drop_duplicates(subset="peptide")
            out = ids.merge(out, on="peptide", how="inner")
        elif pep_df is not None:
            legacy = pep_df[["peptide", "protein_id", "gene_name", "gene_id"]].drop_duplicates(
                subset="peptide"
            )
            out = legacy.merge(out, on="peptide", how="inner")
    elif args.format == "refs":
        out = aggregate_per_pmhc_with_refs(hits)
        if gene_ident_cols:
            ids = hits[["peptide", *gene_ident_cols]].drop_duplicates(subset="peptide")
            out = ids.merge(out, on="peptide", how="inner")
        elif pep_df is not None:
            legacy = pep_df[["peptide", "protein_id", "gene_name", "gene_id"]].drop_duplicates(
                subset="peptide"
            )
            out = legacy.merge(out, on="peptide", how="inner")
    else:
        # raw: observations already carry gene columns; only the raw-CSV
        # path needs the ProteomeIndex join for flanking context.
        if pep_df is not None and not gene_ident_cols:
            legacy = pep_df[
                [
                    "peptide",
                    "protein_id",
                    "gene_name",
                    "gene_id",
                    "position",
                    "n_flank",
                    "c_flank",
                ]
            ].drop_duplicates(subset=["peptide", "protein_id", "position"])
            out = hits.merge(legacy, on="peptide", how="left")
        else:
            out = hits

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
