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
import functools
import sys
from typing import TYPE_CHECKING

from .cli_common import (
    add_iedb_cedar_args,
    add_predictor_arg,
)
from .cli_common import (
    parse_lengths as _parse_lengths,
)
from .cli_common import (
    split_csv as _split_csv,
)
from .cli_common import (
    warn_ignored_flags as _warn_ignored_flags,
)

if TYPE_CHECKING:
    import pandas as pd

_SUPPORTED_FORMATS = ("peptides", "pmhc", "refs", "raw", "table")
# Allele-resolution levels, most specific first.  Mirrors the leading prefix of
# hitlist's ``ALLELE_RESOLUTION_ORDER`` (``unresolved`` is omitted: asking for it
# keeps everything).  Kept as a literal so ``--help`` does not pay for importing
# hitlist.curation; a test guards it against upstream drift.
_SUPPORTED_RESOLUTIONS = ("four_digit", "donor_set", "two_digit", "serological", "class_only")


# Class-default peptide length windows used when --lengths is omitted.
# Intentionally wider than textbook 8-11 / 12-25 to catch non-canonical
# ligands hitlist actually curates (phospho-extended class I, long class
# II tails).  An explicit --lengths overrides these.
_CLASS_I_DEFAULT_LENGTHS = tuple(range(8, 16))  # 8-15
_CLASS_II_DEFAULT_LENGTHS = tuple(range(12, 46))  # 12-45


def _parse_serotypes(value: str) -> list[str]:
    """argparse type for --serotype: comma-separated serotype labels.

    Validates each token as a readable serotype here, at parse time, rather
    than leaving it to fail deep inside ``_filter_by_serotype`` after gene
    resolution and the full observations/evidence-index load have already
    run — mirrors ``--lengths``' ``_parse_lengths``.  ``_filter_by_serotype``
    re-derives the same keys to do the actual matching; this is a fail-fast
    duplicate of that same cheap check.
    """
    from .mhc import serotype_key

    tokens = _split_csv(value)
    unresolved = sorted(token for token in tokens if serotype_key(token) is None)
    if unresolved:
        raise argparse.ArgumentTypeError(
            f"--serotype could not read {', '.join(repr(t) for t in unresolved)} as a "
            "serotype. Use a serological designation such as A2, A24 or Bw4; pass "
            "molecular alleles to --allele instead."
        )
    return tokens


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
            "with `tsarina build observations` — it auto-builds on first use otherwise."
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
        "--alleles",
        "--hla",
        dest="allele",
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
        type=_parse_serotypes,
        default=[],
        help=(
            "Comma-separated serotype labels to keep (e.g. A2,A24,Bw4); case and "
            "the HLA- prefix are optional. Matches any restriction carrying that "
            "serotype, public epitopes and donor allele sets included; narrow to "
            "single-molecule restrictions with --min-resolution."
        ),
    )
    p.add_argument(
        "--min-resolution",
        choices=_SUPPORTED_RESOLUTIONS,
        default=None,
        help=(
            "Minimum allele resolution to keep, read from hitlist's stored "
            "annotation. Choices are ordered most specific first."
        ),
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
            "Include IEDB/CEDAR binding-assay rows in addition to positive MS "
            "observations (default: positive MS only).  On the cached-index path "
            "this routes through hitlist's load_all_evidence() union; rows are "
            "tagged with an evidence_kind column ('ms' vs 'binding'). Does not "
            "admit structural, unknown, or negative-MS assays."
        ),
    )
    p.add_argument(
        "--include-other-assays",
        action="store_true",
        help=(
            "Include structural, non-MS ligand, unknown, and negative-MS assays "
            "in addition to positive MS observations (kept separate from "
            "--include-binding-assays)."
        ),
    )
    p.add_argument(
        "--mono-allelic-only",
        action="store_true",
        help=(
            "Keep only hits with direct mono-allelic evidence, i.e. observed on a "
            "cell line expressing a single HLA allele.  Excludes multi-allelic "
            "studies where the allele was assigned by a predictor (NetMHCpan / "
            "MHCflurry) rather than observed directly."
        ),
    )
    p.add_argument(
        "--format",
        choices=_SUPPORTED_FORMATS,
        default="pmhc",
        help=(
            "Output aggregation, default pmhc. 'peptides' emits per-peptide "
            "summaries with deduplicated ms_alleles, ms_exact_alleles, and "
            "ms_restriction_sets. 'refs' adds per-pMHC ms_pmhc_ref_count, "
            "ms_pmhc_pmids, ms_pmhc_tissues, ms_pmhc_diseases, "
            "ms_pmhc_cell_lines, and ms_pmhc_in_cancer / "
            "ms_pmhc_in_healthy_tissue flags. 'table' renders a compact "
            "terminal table sorted by evidence strength."
        ),
    )
    p.add_argument(
        "--predict",
        action="store_true",
        help="Also score observed (peptide, allele) pairs via topiary.",
    )
    add_predictor_arg(p, context="--predict")
    add_iedb_cedar_args(p)
    p.add_argument(
        "--skip-ms-evidence",
        action="store_true",
        help=(
            "Mode switch: skip IEDB/CEDAR and output only the proteome-enumerated "
            "peptide menu. Filtering/aggregation flags (--allele, --serotype, "
            "--format, --predict, …) are ignored."
        ),
    )
    p.add_argument(
        "--healthy-tissue",
        action="store_true",
        help=(
            "Mode switch: output this gene's per-peptide healthy-SOMATIC-tissue "
            "MS hits (tissue, allele, allele_set, provenance, pmid, vital_organ "
            "flag) instead of the cancer-MS aggregation. Reproductive/thymic hits, "
            "expected for CTAs, are excluded unless --sex includes same-sex adult "
            "reproductive tract hits. Only --gene/--uniprot, --mhc-class, "
            "--species, --sex, --sex-eligibility, and --output apply; other flags "
            "are ignored. See cta_healthy_tissue_ms_hits()."
        ),
    )
    p.add_argument(
        "--sex",
        choices=("female", "male"),
        default=None,
        help="Optional patient sex for --healthy-tissue: include same-sex adult reproductive tract MS hits under all_adult.",
    )
    p.add_argument(
        "--sex-eligibility",
        choices=("all_adult", "sex_organ_only"),
        default="all_adult",
        help="Eligibility policy when --sex is used with --healthy-tissue, default all_adult.",
    )
    p.add_argument(
        "--ms-evidence-scope",
        choices=("all", "requested_mhc", "any_mhc"),
        default="all",
        help=(
            "Control MHC restriction attribution when --allele is provided: "
            "'all' (default), 'requested_mhc' (retain exact/donor-set matches "
            "to requested --allele), or 'any_mhc' (keep all exact-peptide MS "
            "observations while annotating restriction tier in raw output)."
        ),
    )
    p.add_argument(
        "--top",
        type=int,
        default=None,
        help="Limit output to the top N rows (after evidence-strength sorting on pmhc/refs/table).",
    )
    p.add_argument(
        "-o",
        "--output",
        default=None,
        help="Write CSV (or table text for --format table) to this path (default: stdout).",
    )
    return p


@functools.lru_cache(maxsize=1)
def _cached_proteome_index(ensembl_release: int, lengths: tuple[int, ...]):
    # Only one (release, lengths) is used per process. The index is
    # ~8-15 GB resident, so caching more than one entry would double
    # memory without improving hit rate.
    from hitlist.proteome import ProteomeIndex

    return ProteomeIndex.from_ensembl(release=ensembl_release, lengths=lengths, verbose=False)


def _enumerate_gene_peptides(
    gene: str,
    ensembl_release: int,
    lengths: tuple[int, ...],
) -> pd.DataFrame:
    """Build a DataFrame of (peptide, protein_id, gene_name, gene_id, position, n_flank, c_flank)."""
    idx = _cached_proteome_index(ensembl_release, lengths)

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
    import pandas as pd

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


def _require_annotation_column(hits: pd.DataFrame, column: str, flag: str) -> None:
    """Fail loudly when a filter needs an annotation column the index lacks.

    These columns are written by hitlist alongside ``mhc_restriction``, so a
    frame without them is a legacy artifact and the fix is a rebuild rather
    than a different query.
    """
    if column not in hits.columns:
        raise ValueError(
            f"{flag} needs hitlist's `{column}` column, which this observations "
            "index does not carry. Rebuild it with "
            "`tsarina build observations --force`."
        )


def _filter_by_allele(hits: pd.DataFrame, alleles: list[str]) -> pd.DataFrame:
    if not alleles:
        return hits
    from .mhc import mhc_restriction_matches_any, normalize_mhc_restriction_set

    wanted = normalize_mhc_restriction_set(alleles)
    return hits[
        hits["mhc_restriction"].map(lambda value: mhc_restriction_matches_any(value, wanted))
    ].copy()


def _filter_by_serotype(hits: pd.DataFrame, serotypes: list[str]) -> pd.DataFrame:
    """Filter observations to restrictions that name one of ``serotypes``.

    hitlist records every serotype a restriction belongs to in the
    semicolon-joined ``serotypes`` column — the locus serotype (``HLA-A23``)
    plus any public epitope it carries (``HLA-Bw4``) — written in the same
    annotation pass that produced ``mhc_restriction``.  Both the query and the
    stored tokens are reduced to a :func:`~tsarina.mhc.serotype_key`, so
    mhcgnomes decides what each designation means and the two sides agree
    across case, the ``HLA-`` prefix, and split serotypes.

    Matching is membership, exactly as ``--allele`` treats a semicolon-joined
    restriction: a donor set whose typed alleles carry the serotype matches,
    because membership is all the row supports.  Narrowing to restrictions that
    name one molecule is what ``--min-resolution four_digit`` already does, so
    this filter does not second-guess it.
    """
    if not serotypes:
        return hits
    from .mhc import serotype_key, serotype_keys

    requested = {token: serotype_key(token) for token in serotypes if token.strip()}
    if not requested:
        return hits
    unresolved = sorted(token for token, key in requested.items() if key is None)
    if unresolved:
        raise ValueError(
            f"--serotype could not read {', '.join(repr(t) for t in unresolved)} as a "
            "serotype. Use a serological designation such as A2, A24 or Bw4; pass "
            "molecular alleles to --allele instead."
        )
    wanted = set(requested.values())
    if hits.empty:
        return hits
    _require_annotation_column(hits, "serotypes", "--serotype")
    names = hits["serotypes"].astype("string").fillna("")
    # Map over the distinct cell values once, not once per row: the corpus
    # has far fewer distinct ``serotypes`` strings than rows, and computing
    # the intersection per unique value first is markedly faster than
    # re-running it under every row via a plain ``.map(lambda ...)``.
    match_by_cell = {cell: bool(wanted & serotype_keys(cell)) for cell in names.unique()}
    return hits[names.map(match_by_cell)].copy()


def _apply_min_resolution(hits: pd.DataFrame, min_resolution: str | None) -> pd.DataFrame:
    """Keep restrictions at least as specific as ``min_resolution``.

    Reads hitlist's stored ``allele_resolution`` instead of reclassifying
    ``mhc_restriction``.  The stored label comes from the annotation pass that
    produced the restriction, so it accounts for the study's MHC species
    context and for class-only rows promoted to a donor set; reclassifying the
    final restriction string alone sees neither.
    """
    if min_resolution is None or hits.empty:
        return hits
    _require_annotation_column(hits, "allele_resolution", "--min-resolution")
    from hitlist.curation import allele_resolution_rank

    min_rank = allele_resolution_rank(min_resolution)
    ranks = hits["allele_resolution"].astype("string").fillna("").map(allele_resolution_rank)
    return hits[ranks <= min_rank].copy()


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

    # Flags that only apply to the default cancer-MS aggregation path; warn if
    # the user combines them with a mode switch that ignores them.
    def _set_aggregation_flags() -> list[str]:
        flags = []
        if args.allele:
            flags.append("--allele")
        if args.serotype:
            flags.append("--serotype")
        if args.min_resolution is not None:
            flags.append("--min-resolution")
        if args.include_binding_assays:
            flags.append("--include-binding-assays")
        if args.mono_allelic_only:
            flags.append("--mono-allelic-only")
        if args.predict:
            flags.append("--predict")
        if args.format != "pmhc":
            flags.append("--format")
        return flags

    # ── 1b. Healthy-tissue safety screen ───────────────────────────────
    # Short-circuit: peptide x tissue x allele healthy-somatic MS hits for
    # the gene, with a vital-organ flag — the per-patient safety view.
    if args.healthy_tissue:
        ignored = _set_aggregation_flags()
        if args.lengths is not None:
            ignored.append("--lengths")
        if args.skip_ms_evidence:
            ignored.append("--skip-ms-evidence")
        if args.iedb_path is not None:
            ignored.append("--iedb")
        if args.cedar_path is not None:
            ignored.append("--cedar")
        _warn_ignored_flags("--healthy-tissue", ignored)
        from .ms_evidence import cta_healthy_tissue_ms_hits

        ht_kwargs: dict[str, object] = {"mhc_class": args.mhc_class, "mhc_species": mhc_species}
        if getattr(args, "sex", None) is not None:
            ht_kwargs["patient_sex"] = args.sex
            ht_kwargs["sex_eligibility"] = getattr(args, "sex_eligibility", "all_adult")
        hits = cta_healthy_tissue_ms_hits(gene, **ht_kwargs)
        if hits.empty:
            print(f"No healthy-somatic-tissue MS hits for {gene}.", file=sys.stderr)
        _write(args.output, hits)
        return

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
        _warn_ignored_flags("--skip-ms-evidence", _set_aggregation_flags())
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
        from .ms_evidence import filter_positive_ms_observations

        hits = filter_positive_ms_observations(
            hits,
            include_binding=args.include_binding_assays,
            include_other=getattr(args, "include_other_assays", False),
        )
    else:
        from .evidence_store import load_evidence_store_ms_hits, resolve_evidence_package

        resolved_package = resolve_evidence_package(getattr(args, "evidence_package", None))
        length_min = min(resolved_lengths) if resolved_lengths else None
        length_max = max(resolved_lengths) if resolved_lengths else None
        include_other = getattr(args, "include_other_assays", False)

        if resolved_package is not None:
            hits = load_evidence_store_ms_hits(
                resolved_package,
                gene_name=gene,
                mhc_class=args.mhc_class,
                mhc_species=mhc_species,
                drop_binding_assays=not args.include_binding_assays,
                include_binding=args.include_binding_assays,
                include_other=include_other,
                length_min=length_min,
                length_max=length_max,
            )
        else:
            # Fast path: direct gene_name pushdown through hitlist's mappings sidecar.
            from .indexing import ensure_index_built

            ensure_index_built()

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
            if include_other:
                import pandas as pd
                from hitlist.observations import load_other_assays

                other = load_other_assays(
                    gene_name=gene,
                    mhc_class=args.mhc_class,
                    species=mhc_species,
                    length_min=length_min,
                    length_max=length_max,
                )
                if not other.empty:
                    other = other.copy()
                    if "evidence_kind" not in other.columns:
                        other["evidence_kind"] = "other"
                    if "evidence_kind" not in hits.columns and not hits.empty:
                        hits = hits.copy()
                        hits["evidence_kind"] = "ms"
                    hits = pd.concat([hits, other], ignore_index=True)
            # hitlist 1.30.46 split gene_names / gene_ids / protein_ids out of
            # observations.parquet (#238) — they're now derived on demand from
            # the peptide_mappings sidecar.  No-projection loads don't get them
            # back automatically, so re-attach here using hitlist's own
            # aggregator so downstream gene_ident_cols logic keeps working on
            # both pre- and post-1.30.46 parquets.
            #
            # Test note: any test that mocks ``hitlist.observations.load_observations``
            # (or ``load_all_evidence``) with a non-empty frame lacking
            # ``gene_names`` must also mock ``hitlist.mappings.load_peptide_mappings``
            # — otherwise the real loader raises ``FileNotFoundError`` against
            # a CI environment without a built sidecar.
            if not hits.empty and "gene_names" not in hits.columns:
                from hitlist.mappings import (
                    annotate_observations_with_genes,
                    load_peptide_mappings,
                )

                mappings = load_peptide_mappings(
                    peptide=hits["peptide"].unique().tolist(),
                    columns=["peptide", "gene_name", "gene_id", "protein_id"],
                )
                hits = annotate_observations_with_genes(hits, mappings)
        if resolved_lengths and not hits.empty:
            hits = hits[hits["peptide"].str.len().isin(resolved_lengths)].copy()
        hits = _apply_min_resolution(hits, args.min_resolution)

    ms_evidence_scope = getattr(args, "ms_evidence_scope", "all")
    try:
        if ms_evidence_scope == "any_mhc" and args.allele:
            from .ms_attribution import attribute_ms_evidence

            hits = attribute_ms_evidence(
                hits[["peptide"]].drop_duplicates(),
                hits,
                requested_alleles=args.allele,
                evidence_scope="any_mhc",
                evidence_loaded=True,
            )
        elif ms_evidence_scope == "requested_mhc":
            if not args.allele:
                raise ValueError("--ms-evidence-scope requested_mhc requires --allele")
            from .ms_attribution import attribute_ms_evidence

            hits = attribute_ms_evidence(
                hits[["peptide"]].drop_duplicates(),
                hits,
                requested_alleles=args.allele,
                evidence_scope="requested_mhc",
                evidence_loaded=True,
            )
        else:
            hits = _filter_by_allele(hits, args.allele)
        hits = _filter_by_serotype(hits, args.serotype)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if args.mono_allelic_only and not hits.empty and "is_monoallelic" in hits.columns:
        hits = hits[hits["is_monoallelic"]].copy()

    if hits.empty:
        # An empty result is otherwise indistinguishable from "unknown gene" or
        # "command didn't run" — say so on stderr (output stays pipe-clean).
        print(
            f"No public MS hits for gene '{gene}' with the current filters "
            f"(check the spelling; try --include-binding-assays or --skip-ms-evidence).",
            file=sys.stderr,
        )
        import pandas as pd

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
        out = _annotate_pmhc_attribution(out, args.allele)
        out = _sort_pmhc_by_evidence(out)
    elif args.format in ("refs", "table"):
        out = aggregate_per_pmhc_with_refs(hits)
        if gene_ident_cols:
            ids = hits[["peptide", *gene_ident_cols]].drop_duplicates(subset="peptide")
            out = ids.merge(out, on="peptide", how="inner")
        elif pep_df is not None:
            legacy = pep_df[["peptide", "protein_id", "gene_name", "gene_id"]].drop_duplicates(
                subset="peptide"
            )
            out = legacy.merge(out, on="peptide", how="inner")
        out = _annotate_pmhc_attribution(out, args.allele)
        out = _sort_pmhc_by_evidence(out)
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

        if args.format in ("pmhc", "refs", "table"):
            peptide_allele = out[["peptide", "mhc_restriction"]].drop_duplicates()
            single_alleles = [
                a for a in peptide_allele["mhc_restriction"].unique().tolist() if ";" not in str(a)
            ]
            if single_alleles:
                scores = score_presentation(
                    peptides=peptide_allele["peptide"].tolist(),
                    alleles=single_alleles,
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

    top = getattr(args, "top", None)
    if top is not None:
        if top <= 0:
            print("Error: --top must be a positive integer.", file=sys.stderr)
            sys.exit(1)
        out = out.head(top)

    if args.format == "table":
        rendered = _format_hits_table(out) + "\n"
        if args.output:
            with open(args.output, "w") as f:
                f.write(rendered)
            print(f"Wrote {len(out)} rows to {args.output}", file=sys.stderr)
        else:
            sys.stdout.write(rendered)
        return

    _write(args.output, out)


def _sort_pmhc_by_evidence(df: pd.DataFrame) -> pd.DataFrame:
    """Order per-pMHC rows by evidence strength (#241): exact single-allele
    before multi-allelic donor sets, then monoallelic evidence, PMID/ref count,
    and hit count."""
    if df.empty or "mhc_restriction" not in df.columns:
        return df
    import pandas as pd

    out = df.copy()
    out["_is_donor_set"] = out["mhc_restriction"].astype(str).str.contains(";", regex=False)
    if "ms_pmhc_has_mono_evidence" in out.columns:
        mono = out["ms_pmhc_has_mono_evidence"].fillna(False).astype(bool)
    elif "ms_pmhc_mono_hit_count" in out.columns:
        mono = pd.to_numeric(out["ms_pmhc_mono_hit_count"], errors="coerce").fillna(0) > 0
    else:
        mono = pd.Series(False, index=out.index)
    out["_mono_rank"] = (~mono).astype(int)
    pmid_col = (
        "ms_pmhc_pmid_count"
        if "ms_pmhc_pmid_count" in out.columns
        else "ms_pmhc_ref_count"
        if "ms_pmhc_ref_count" in out.columns
        else None
    )
    out["_pmid_rank"] = (
        -pd.to_numeric(out[pmid_col], errors="coerce").fillna(0).astype(int)
        if pmid_col is not None
        else 0
    )
    out["_hit_rank"] = (
        -pd.to_numeric(out["ms_pmhc_hit_count"], errors="coerce").fillna(0).astype(int)
        if "ms_pmhc_hit_count" in out.columns
        else 0
    )
    sort_cols = ["_is_donor_set", "_mono_rank", "_pmid_rank", "_hit_rank"]
    if "peptide" in out.columns:
        sort_cols.append("peptide")
    sort_cols.append("mhc_restriction")
    out = out.sort_values(sort_cols, kind="mergesort").drop(
        columns=["_is_donor_set", "_mono_rank", "_pmid_rank", "_hit_rank"]
    )
    return out.reset_index(drop=True)


def _annotate_pmhc_attribution(df: pd.DataFrame, requested_alleles: list[str]) -> pd.DataFrame:
    """Tag per-pMHC rows with exact vs donor-set attribution when --allele is given (#241)."""
    if df.empty or not requested_alleles or "mhc_restriction" not in df.columns:
        return df
    from .mhc import normalize_mhc_restriction_set

    wanted = normalize_mhc_restriction_set(requested_alleles)
    out = df.copy()

    def _classify(restriction: object) -> str:
        text = str(restriction or "").strip()
        if not text:
            return "other_mhc_only"
        tokens = [tok.strip() for tok in text.split(";") if tok.strip()]
        normalized = normalize_mhc_restriction_set(tokens)
        if not (normalized & wanted):
            return "other_mhc_only"
        if len(tokens) == 1:
            return "exact_requested_mhc"
        return "donor_set_ambiguous"

    out["ms_attribution_tier"] = out["mhc_restriction"].map(_classify)
    return out


def _abbreviate_restriction(value: object, max_alleles: int = 3) -> str:
    text = str(value or "").strip()
    if ";" not in text:
        return text
    tokens = [tok.strip() for tok in text.split(";") if tok.strip()]
    if len(tokens) <= max_alleles:
        return ";".join(tokens)
    return f"{';'.join(tokens[:max_alleles])};... (+{len(tokens) - max_alleles})"


def _format_hits_table(df: pd.DataFrame) -> str:
    """Render a compact fixed-width terminal table for ``tsarina hits --format table``."""
    import pandas as pd

    if df.empty:
        return "(no hits)"

    attr_short = {
        "exact_requested_mhc": "exact",
        "donor_set_ambiguous": "donor_set",
        "other_mhc_only": "other_mhc",
    }

    def _safety(row: pd.Series) -> str:
        in_cancer = (
            bool(row.get("ms_pmhc_in_cancer")) if pd.notna(row.get("ms_pmhc_in_cancer")) else False
        )
        in_healthy = (
            bool(row.get("ms_pmhc_in_healthy_tissue"))
            if pd.notna(row.get("ms_pmhc_in_healthy_tissue"))
            else False
        )
        if in_healthy and in_cancer:
            return "cancer+healthy!"
        if in_healthy:
            return "healthy!"
        if in_cancer:
            return "cancer"
        return "unclassified"

    has_attr = "ms_attribution_tier" in df.columns
    has_safety = bool({"ms_pmhc_in_cancer", "ms_pmhc_in_healthy_tissue"} & set(df.columns))
    has_predict = "presentation_percentile" in df.columns

    headers = ["peptide", "genes", "mhc_restriction"]
    if has_attr:
        headers.append("ms_attr")
    headers.extend(["hits", "pmids", "mono"])
    if has_safety:
        headers.append("ms_safety")
    if has_predict:
        headers.extend(["pct_rank", "affinity_nM"])

    rows: list[list[str]] = []
    for _, row in df.iterrows():
        pmids_val = row.get("ms_pmhc_pmid_count")
        if pd.isna(pmids_val):
            pmids_val = row.get("ms_pmhc_ref_count")
        if pd.isna(pmids_val):
            pmid_str = str(row.get("ms_pmhc_pmids") or "").strip()
            pmids_val = len([p for p in pmid_str.split(";") if p.strip()]) if pmid_str else 0
        mono_hits = row.get("ms_pmhc_mono_hit_count")
        mono_str = (
            f"yes ({int(mono_hits)})"
            if pd.notna(mono_hits) and int(mono_hits) > 0
            else "yes"
            if bool(row.get("ms_pmhc_has_mono_evidence"))
            else "no"
        )
        entry = [
            str(row.get("peptide") or ""),
            str(row.get("gene_names") or row.get("gene_name") or ""),
            _abbreviate_restriction(row.get("mhc_restriction")),
        ]
        if has_attr:
            raw_attr = str(row.get("ms_attribution_tier") or "")
            entry.append(attr_short.get(raw_attr, raw_attr))
        entry.extend(
            [
                str(int(row.get("ms_pmhc_hit_count") or 0)),
                str(int(pmids_val or 0)),
                mono_str,
            ]
        )
        if has_safety:
            entry.append(_safety(row))
        if has_predict:
            pct = row.get("presentation_percentile")
            aff = row.get("affinity_nm")
            entry.extend(
                [
                    "" if pd.isna(pct) else f"{float(pct):.2f}",
                    "" if pd.isna(aff) else f"{float(aff):.1f}",
                ]
            )
        rows.append(entry)

    widths = [max(len(h), *(len(r[i]) for r in rows)) for i, h in enumerate(headers)]
    lines = [
        "  ".join(h.ljust(w) for h, w in zip(headers, widths)),
        "  ".join("-" * w for w in widths),
    ]
    lines.extend("  ".join(cell.ljust(w) for cell, w in zip(r, widths)) for r in rows)
    return "\n".join(lines)


def _write(output: str | None, df: pd.DataFrame) -> None:
    if output:
        df.to_csv(output, index=False)
        print(f"Wrote {len(df)} rows to {output}", file=sys.stderr)
    else:
        df.to_csv(sys.stdout, index=False)
