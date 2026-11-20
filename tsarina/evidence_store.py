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

"""Hitlist 2.0 verified ``EvidenceStore`` adapter for Tsarina MS queries.

A verified Hitlist 2.x evidence package binds typed tables
(``typed/mhc_observations.parquet``, ``typed/assay_context.parquet``) and
native lossless parquets (``native/observations.parquet``,
``native/binding.parquet``, ``native/peptide_mappings.parquet``) under one
SHA256-verified manifest.  Opening a package through
:class:`hitlist.EvidenceStore` verifies the complete manifest before any row
is served.
"""

from __future__ import annotations

import json
import os
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

if TYPE_CHECKING:
    from hitlist import EvidenceStore

#: Environment variables that configure a default verified Hitlist 2.x package.
EVIDENCE_PACKAGE_ENV_VARS: tuple[str, ...] = (
    "TSARINA_EVIDENCE_PACKAGE",
    "HITLIST_EVIDENCE_PACKAGE",
)

#: Typed ``mhc_positive_ms`` / ``mhc_binding`` / ``mhc_other`` fields selected
#: from ``EvidenceStore.query``; includes native row pointers for exact
#: lossless enrichment from the manifest-verified ``native/`` partition.
_TYPED_SELECT_FIELDS: tuple[str, ...] = (
    "observation_id",
    "peptide_sequence",
    "peptide_modifications_raw",
    "has_ptm",
    "peptide_extended",
    "is_non_peptide_ligand",
    "is_potential_contaminant",
    "reported_mhc",
    "tested_mhc",
    "candidate_mhc_set_raw",
    "peptide_mhc_restriction",
    "peptide_restriction_status",
    "restriction_evidence",
    "assay_method",
    "assay_readout",
    "assay_modality",
    "assay_modality_basis",
    "outcome_status",
    "source_dataset",
    "source_file",
    "source_assay_id",
    "source_reference_iri",
    "source_pmid",
    "native_partition",
    "native_row_index",
)

#: Mapping from typed ``EvidenceStore`` columns to canonical Tsarina/Hitlist
#: observation column names (used as fallback when a native column is absent).
_TYPED_TO_CANONICAL: dict[str, str] = {
    "observation_id": "provenance_id",
    "peptide_sequence": "peptide",
    "peptide_modifications_raw": "peptide_modifications",
    "reported_mhc": "mhc_restriction",
    "candidate_mhc_set_raw": "mhc_allele_set",
    "assay_readout": "response_measured",
    "assay_modality_basis": "assay_modality_source",
    "outcome_status": "ms_outcome_status",
    "source_dataset": "source",
    "source_file": "supplementary_file",
    "source_assay_id": "assay_iri",
    "source_reference_iri": "reference_iri",
    "source_pmid": "pmid",
}

_EVIDENCE_KIND_BY_TYPE: dict[str, str] = {
    "mhc_positive_ms": "ms",
    "mhc_binding": "binding",
    "mhc_other": "other",
}


def resolve_evidence_package(package_path: str | Path | None = None) -> Path | None:
    """Resolve an explicit or environment-configured Hitlist 2.x package path."""
    if package_path is not None:
        return Path(package_path)
    for var in EVIDENCE_PACKAGE_ENV_VARS:
        value = os.environ.get(var, "").strip()
        if value:
            return Path(value)
    return None


class _StandaloneMhcStore:
    """Thin verified wrapper for standalone MHC packages (schema_version == 2).

    Full multi-modality packages open through :class:`hitlist.EvidenceStore`;
    standalone MHC packages written by :func:`hitlist.write_mhc_evidence_package`
    are verified by :func:`hitlist.verify_mhc_evidence_package` and expose the
    exact same ``query`` / ``count`` / ``bindings`` / ``schemas`` interface for
    ``mhc_positive_ms``, ``mhc_binding``, and ``mhc_other``.
    """

    def __init__(self, directory: Path, manifest: dict):
        self.directory = directory
        self.manifest = manifest
        obs_path = directory / "typed" / "mhc_observations.parquet"
        obs_schema = pq.ParquetFile(obs_path).schema_arrow
        self.schemas = dict.fromkeys(("mhc_positive_ms", "mhc_binding", "mhc_other"), obs_schema)
        self.bindings = dict.fromkeys(self.schemas)
        self._obs_path = obs_path

    def _manifest_hash(self) -> str:
        from hitlist.provenance import file_digest

        return file_digest(self.directory / "manifest.json")["sha256"]

    def _filter_expr(self, evidence_type: str, where: dict | None):
        expr = ds.field("evidence_type") == evidence_type
        for col, cond in (where or {}).items():
            if isinstance(cond, dict):
                op, val = cond["op"], cond["value"]
            else:
                op, val = "eq", cond
            field = ds.field(col)
            if op == "eq":
                part = field == val
            elif op == "in":
                part = field.isin(list(val)) if val else ds.scalar(False)
            else:
                raise ValueError(f"Unsupported predicate op: {op}")
            expr = expr & part
        return expr

    def count(self, evidence_type: str, *, where: dict | None = None) -> int:
        return ds.dataset(self._obs_path, format="parquet").count_rows(
            filter=self._filter_expr(evidence_type, where)
        )

    def query(
        self,
        evidence_type: str,
        *,
        select=None,
        where: dict | None = None,
        joins=(),
        limit: int | None = 1000,
    ):
        from hitlist.evidence_store import QueryResult

        cols = list(select) if select is not None else list(self.schemas[evidence_type].names)
        table = ds.dataset(self._obs_path, format="parquet").to_table(
            columns=cols,
            filter=self._filter_expr(evidence_type, where),
        )
        if limit is not None:
            table = table.slice(0, limit)
        return QueryResult(
            evidence_type=evidence_type,
            grain="assay_result",
            source_population="mhc_v4",
            manifest_sha256=self._manifest_hash(),
            catalog_schema_version=int(self.manifest.get("catalog_schema_version", 3)),
            table=table,
        )


@lru_cache(maxsize=4)
def _open_verified_store(path_str: str, mtime_ns: int, size: int):
    from hitlist import EvidenceStore, verify_mhc_evidence_package

    directory = Path(path_str)
    envelope = json.loads((directory / "manifest.json").read_text())
    if envelope.get("schema_version") == 2 and "typed/mhc_observations.parquet" in envelope.get(
        "files", {}
    ):
        manifest = verify_mhc_evidence_package(directory)
        return _StandaloneMhcStore(directory, manifest)
    return EvidenceStore.open(path_str)


def clear_evidence_store_cache() -> None:
    """Clear the verified ``EvidenceStore`` cache (primarily for tests)."""
    _open_verified_store.cache_clear()


def open_evidence_store(package_path: str | Path) -> EvidenceStore:
    """Open and cryptographically verify a Hitlist 2.x ``EvidenceStore`` package.

    Cached by ``(resolved_path, manifest_mtime_ns, manifest_size)`` so repeated
    queries against the same verified package do not re-hash gigabyte files
    while any manifest change triggers re-verification.
    """
    resolved = Path(package_path).resolve()
    manifest_path = resolved / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Hitlist EvidenceStore package manifest not found: {manifest_path}"
        )
    stat = manifest_path.stat()
    return _open_verified_store(str(resolved), stat.st_mtime_ns, stat.st_size)


def evidence_store_provenance(package_path: str | Path) -> dict:
    """Return cryptographic provenance metadata for a verified Hitlist 2.x package."""
    import hitlist
    import mhcgnomes

    store = open_evidence_store(package_path)
    manifest_sha256 = store._manifest_hash()
    info: dict = {
        "evidence_backend": "hitlist_evidence_store",
        "evidence_package_path": str(store.directory),
        "evidence_package_manifest_sha256": manifest_sha256,
        "evidence_package_schema_version": store.manifest.get("schema_version"),
        "evidence_package_catalog_schema_version": store.manifest.get("catalog_schema_version"),
        "evidence_package_hitlist_version": store.manifest.get("hitlist_version"),
        "hitlist_version": hitlist.__version__,
        "mhcgnomes_version": mhcgnomes.__version__,
    }
    if "mhc_positive_ms" in store.bindings:
        info["mhc_positive_ms_rows"] = int(store.count("mhc_positive_ms"))
    if "mhc_binding" in store.bindings:
        info["mhc_binding_rows"] = int(store.count("mhc_binding"))
    return info


def _take_native_rows(
    parquet_path: Path,
    row_indices: list[int],
    columns: list[str] | None = None,
) -> pd.DataFrame:
    """Read specific 0-based rows from a parquet file by row-group take.

    Avoids loading the full native table into memory: only row groups containing
    at least one requested row index are read.
    """
    if not row_indices:
        return pd.DataFrame()
    pf = pq.ParquetFile(parquet_path)
    available = pf.schema_arrow.names
    selected_cols = (
        [c for c in columns if c in available] if columns is not None else list(available)
    )
    if not selected_cols:
        return pd.DataFrame(index=range(len(row_indices)))

    # Map each global row index to its output position(s).
    positions_by_row: dict[int, list[int]] = {}
    for pos, row_idx in enumerate(row_indices):
        positions_by_row.setdefault(int(row_idx), []).append(pos)

    tables: list[pa.Table] = []
    out_positions: list[int] = []
    rg_start = 0
    for rg_idx in range(pf.num_row_groups):
        rg_len = pf.metadata.row_group(rg_idx).num_rows
        rg_end = rg_start + rg_len
        local_indices: list[int] = []
        rg_positions: list[int] = []
        for global_idx in sorted(k for k in positions_by_row if rg_start <= k < rg_end):
            local_idx = global_idx - rg_start
            for pos in positions_by_row[global_idx]:
                local_indices.append(local_idx)
                rg_positions.append(pos)
        if local_indices:
            rg_table = pf.read_row_group(rg_idx, columns=selected_cols)
            tables.append(rg_table.take(local_indices))
            out_positions.extend(rg_positions)
        rg_start = rg_end

    if not tables:
        return pd.DataFrame(columns=selected_cols)

    combined = pa.concat_tables(tables) if len(tables) > 1 else tables[0]
    frame = combined.to_pandas()
    frame["_out_pos"] = out_positions
    return frame.sort_values("_out_pos").drop(columns=["_out_pos"]).reset_index(drop=True)


def _peptides_for_genes_in_package(
    package_dir: Path,
    gene_name: str | list[str],
) -> tuple[list[str], pd.DataFrame]:
    """Resolve ``gene_name`` to peptides using the package's verified mappings parquet."""
    mappings_path = package_dir / "native" / "peptide_mappings.parquet"
    if not mappings_path.exists():
        raise FileNotFoundError(
            f"Gene lookup requires native/peptide_mappings.parquet inside {package_dir}"
        )
    names = [gene_name] if isinstance(gene_name, str) else list(gene_name)
    if not names:
        return [], pd.DataFrame(columns=["peptide", "gene_name", "gene_id", "protein_id"])

    dataset = ds.dataset(mappings_path, format="parquet")
    available = set(dataset.schema.names)
    read_cols = [c for c in ("peptide", "gene_name", "gene_id", "protein_id") if c in available]
    table = dataset.to_table(columns=read_cols, filter=ds.field("gene_name").isin(names))
    mappings = table.to_pandas()
    if mappings.empty:
        return [], mappings
    peptides = sorted(mappings["peptide"].dropna().astype(str).unique().tolist())
    return peptides, mappings


def _query_one_type(
    store: EvidenceStore,
    evidence_type: str,
    where: dict,
    *,
    include_assay_context: bool,
) -> pd.DataFrame:
    schema_names = set(store.schemas[evidence_type].names)
    select = [f for f in _TYPED_SELECT_FIELDS if f in schema_names]
    joins: tuple[str, ...] = ()
    if include_assay_context and "assay_context" in store.bindings:
        joins = ("assay_context",)
        ctx_schema = set(store.schemas["assay_context"].names)
        for ctx_field in ("sample_attribution", "acquisition_mode", "instrument"):
            if ctx_field in ctx_schema:
                select.append(f"assay_context.{ctx_field}")

    result = store.query(
        evidence_type,
        select=select,
        where=where,
        joins=joins,
        limit=None,
    )
    typed_df = result.table.to_pandas()
    if typed_df.empty:
        return typed_df

    typed_df["evidence_type"] = evidence_type
    typed_df["evidence_kind"] = _EVIDENCE_KIND_BY_TYPE.get(evidence_type, evidence_type)

    # Enrich with lossless deferred native fields from manifest-verified native/
    # parquets via exact (native_partition, native_row_index) pointers.
    if {"native_partition", "native_row_index"} <= set(typed_df.columns):
        native_parts: list[pd.DataFrame] = []
        for partition, group in typed_df.groupby("native_partition", sort=False):
            native_path = store.directory / "native" / str(partition)
            if not native_path.exists():
                continue
            indices = [int(i) for i in group["native_row_index"].tolist()]
            native_rows = _take_native_rows(native_path, indices)
            native_rows.index = group.index
            native_parts.append(native_rows)
        if native_parts:
            native_df = pd.concat(native_parts).reindex(typed_df.index)
            # Keep every native column and overlay typed-only columns alongside.
            for col in typed_df.columns:
                if col not in native_df.columns:
                    native_df[col] = typed_df[col]
            merged = native_df
        else:
            merged = typed_df
    else:
        merged = typed_df

    # Ensure canonical Tsarina/Hitlist column names exist even when a package
    # omits optional native columns.
    for typed_col, canonical_col in _TYPED_TO_CANONICAL.items():
        if canonical_col not in merged.columns and typed_col in merged.columns:
            merged[canonical_col] = merged[typed_col]

    return merged


def load_evidence_store_ms_hits(
    package_path: str | Path,
    *,
    peptides: set[str] | list[str] | tuple[str, ...] | None = None,
    mhc_class: str | None = "I",
    mhc_species: str | None = "Homo sapiens",
    gene_name: str | list[str] | None = None,
    columns: list[str] | None = None,
    drop_binding_assays: bool = True,
    include_binding: bool = False,
    include_other: bool = False,
    exclude_non_peptide_ligand: bool = True,
    length_min: int | None = None,
    length_max: int | None = None,
) -> pd.DataFrame:
    """Query MHC MS evidence from a verified Hitlist 2.x ``EvidenceStore`` package.

    Opens and hash-verifies the package via :class:`hitlist.EvidenceStore`,
    queries ``mhc_positive_ms`` (and optionally ``mhc_binding`` / ``mhc_other``),
    enriches matched rows from the package's verified ``native/`` parquets,
    applies the same species/class/length filters as
    :func:`tsarina.indexing.load_ms_evidence`, and attaches cryptographic
    package provenance to ``df.attrs``.
    """
    store = open_evidence_store(package_path)
    if "mhc_positive_ms" not in store.bindings:
        raise ValueError(
            f"Hitlist EvidenceStore package at {store.directory} does not materialize "
            "'mhc_positive_ms'"
        )

    provenance = evidence_store_provenance(store.directory)

    gene_mappings: pd.DataFrame | None = None
    peptide_filter: list[str] | None = None
    if gene_name is not None:
        gene_peptides_list, gene_mappings = _peptides_for_genes_in_package(
            store.directory, gene_name
        )
        if not gene_peptides_list:
            empty = pd.DataFrame(columns=columns or [])
            empty.attrs.update(provenance)
            return empty
        peptide_filter = gene_peptides_list

    if peptides is not None:
        requested = sorted({str(p) for p in peptides})
        if peptide_filter is not None:
            requested = sorted(set(requested) & set(peptide_filter))
        if not requested:
            empty = pd.DataFrame(columns=columns or [])
            empty.attrs.update(provenance)
            return empty
        peptide_filter = requested

    where: dict = {}
    if peptide_filter is not None:
        where["peptide_sequence"] = {"op": "in", "value": peptide_filter}
    if exclude_non_peptide_ligand:
        where["is_non_peptide_ligand"] = False

    evidence_types = ["mhc_positive_ms"]
    if include_binding and "mhc_binding" in store.bindings:
        evidence_types.append("mhc_binding")
    if include_other and "mhc_other" in store.bindings:
        evidence_types.append("mhc_other")

    frames = [
        _query_one_type(
            store,
            ev_type,
            where,
            include_assay_context=True,
        )
        for ev_type in evidence_types
    ]
    non_empty = [f for f in frames if not f.empty]
    if not non_empty:
        empty = pd.DataFrame(columns=columns or [])
        empty.attrs.update(provenance)
        return empty

    df = pd.concat(non_empty, ignore_index=True) if len(non_empty) > 1 else non_empty[0]

    if mhc_class is not None and "mhc_class" in df.columns:
        df = df[df["mhc_class"].astype(str) == str(mhc_class)]
    if mhc_species is not None:
        species_col = (
            "mhc_species"
            if "mhc_species" in df.columns
            else ("species" if "species" in df.columns else None)
        )
        if species_col is not None:
            df = df[df[species_col].astype(str) == str(mhc_species)]
    if length_min is not None and "peptide" in df.columns:
        df = df[df["peptide"].astype(str).str.len() >= int(length_min)]
    if length_max is not None and "peptide" in df.columns:
        df = df[df["peptide"].astype(str).str.len() <= int(length_max)]
    if drop_binding_assays and not include_binding and "is_binding_assay" in df.columns:
        df = df[~df["is_binding_assay"].fillna(False).astype(bool)]

    if gene_mappings is not None and not df.empty and "gene_names" not in df.columns:
        from hitlist.mappings import annotate_observations_with_genes

        df = annotate_observations_with_genes(df, gene_mappings)

    df = df.reset_index(drop=True)
    if columns is not None:
        available = [c for c in columns if c in df.columns]
        df = df[available]

    df.attrs.update(provenance)
    return df
