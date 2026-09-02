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

"""Unified observations-index management, layered over hitlist.

The hitlist observations parquet is the canonical indexed MS evidence table.
First-time build takes a few minutes; afterwards queries are sub-second.
This module centralizes the "check-and-build" policy so every query path in
tsarina hits the same cache and emits the same progress message when a build
is actually triggered.
"""

from __future__ import annotations

import json
import sys
from contextlib import suppress
from pathlib import Path

import pandas as pd

_MAPPING_COMPATIBILITY_CONTRACT = "hitlist-1.55.2-length-independent-mappings"
_MAPPING_COMPATIBILITY_FILENAME = ".tsarina-peptide-mappings.json"
_MAPPING_PROBE_SIZE = 128


def _mapping_fingerprint(path: Path) -> dict[str, int]:
    stat = path.stat()
    return {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _mapping_compatibility_path(mapping_path: Path) -> Path:
    return mapping_path.with_name(_MAPPING_COMPATIBILITY_FILENAME)


def _mapping_compatibility_is_recorded(mapping_path: Path) -> bool:
    marker_path = _mapping_compatibility_path(mapping_path)
    try:
        marker = json.loads(marker_path.read_text())
        return marker.get("contract") == _MAPPING_COMPATIBILITY_CONTRACT and marker.get(
            "mapping"
        ) == _mapping_fingerprint(mapping_path)
    except (FileNotFoundError, OSError, TypeError, json.JSONDecodeError):
        return False


def _record_mapping_compatibility(mapping_path: Path) -> None:
    """Remember that this exact sidecar was built or verified under hitlist 1.55.2+."""
    if not mapping_path.exists():
        return
    marker_path = _mapping_compatibility_path(mapping_path)
    temporary = marker_path.with_suffix(marker_path.suffix + ".tmp")
    marker = {
        "contract": _MAPPING_COMPATIBILITY_CONTRACT,
        "mapping": _mapping_fingerprint(mapping_path),
    }
    try:
        temporary.write_text(json.dumps(marker, indent=2) + "\n")
        temporary.replace(marker_path)
    except OSError:
        # A read-only shared cache is still usable. It will simply be probed
        # again in the next process because tsarina cannot persist the result.
        with suppress(OSError):
            temporary.unlink()


def _probe_mapping_compatibility(mapping_path: Path) -> bool:
    """Detect the pre-1.55 sidecar using behavior rather than package metadata.

    hitlist 1.55.2 does not stamp its builder version into mapping metadata.
    Old artifacts map none of the class-II or length-7 observation peptides;
    current artifacts map almost all of both groups. A small observation-backed
    sample therefore distinguishes the artifacts without reading the full
    175 MB mappings table.
    """
    from hitlist.mappings import load_peptide_mappings
    from hitlist.observations import load_observations

    probes = (
        {"mhc_class": "II", "length_min": 12, "length_max": 45},
        {"mhc_class": None, "length_min": 7, "length_max": 7},
    )
    for filters in probes:
        observations = load_observations(columns=["peptide"], **filters)
        peptides = (
            observations["peptide"].dropna().astype(str).drop_duplicates().head(_MAPPING_PROBE_SIZE)
        )
        if peptides.empty:
            continue
        mappings = load_peptide_mappings(
            peptide=peptides.tolist(),
            columns=["peptide"],
        )
        mapped = set(mappings["peptide"].dropna().astype(str)) if not mappings.empty else set()
        if mapped.isdisjoint(set(peptides)):
            return False

    _record_mapping_compatibility(mapping_path)
    return True


def _mapping_cache_is_compatible() -> bool:
    from hitlist.mappings import mappings_path

    path = mappings_path()
    if _mapping_compatibility_is_recorded(path):
        return True
    return _probe_mapping_compatibility(path)


def ensure_index_built(force: bool = False, verbose: bool = True) -> Path:
    """Ensure both the hitlist observations parquet and its peptide_mappings
    sidecar exist; build if either is missing.

    tsarina depends on the sidecar for gene-identifier resolution
    (``annotate_observations_with_genes`` in the cached fast path,
    ``gene_name=`` filter pushdown in ``load_ms_evidence``). A missing or
    pre-hitlist-1.55 sidecar is rebuilt directly from the current observations,
    avoiding a full evidence rescan. Current sidecars are behavior-probed once
    and fingerprinted so later calls only need a metadata check.

    Parameters
    ----------
    force
        Rebuild even if a valid cached index exists.
    verbose
        Emit a stderr message when a build is actually performed (first use
        or forced rebuild).  No output when an existing cache is reused.

    Returns
    -------
    Path
        Path to ``observations.parquet``.
    """
    from hitlist.builder import build_observations
    from hitlist.mappings import build_peptide_mappings, is_mappings_built, mappings_path
    from hitlist.observations import is_built, observations_path

    observations_built = is_built()
    mappings_built = is_mappings_built()

    if force or not observations_built:
        if verbose:
            print(
                "Building hitlist observations index (one-time ~2-5 min; cached afterwards)...",
                file=sys.stderr,
            )
        build_observations(force=force)
        _record_mapping_compatibility(mappings_path())
    elif not mappings_built or not _mapping_cache_is_compatible():
        if verbose:
            reason = "missing" if not mappings_built else "predates hitlist 1.55.2"
            print(
                f"Building hitlist peptide mappings ({reason}; one-time migration)...",
                file=sys.stderr,
            )
        build_peptide_mappings(force=True)
        _record_mapping_compatibility(mappings_path())
    return observations_path()


def load_ms_evidence(
    peptides: set[str] | list[str] | None = None,
    mhc_class: str | None = "I",
    mhc_species: str | None = "Homo sapiens",
    gene_name: str | list[str] | None = None,
    columns: list[str] | None = None,
    auto_build: bool = True,
    drop_binding_assays: bool = True,
) -> pd.DataFrame:
    """Load MS evidence rows from the hitlist observations index.

    Ensures the index exists (building if needed) and runs a pushdown-filtered
    parquet read for ``mhc_class``, ``mhc_species``, ``gene_name``, and the
    peptide set.

    Parameters
    ----------
    peptides
        Peptide sequences to keep.  If None, the species+class filtered slice
        is returned unfiltered by peptide (profile mode).
    mhc_class
        ``"I"``, ``"II"``, or None.  Default ``"I"``.
    mhc_species
        Species filter (default ``"Homo sapiens"``; None disables).
    gene_name
        Filter to rows whose peptide maps to this gene (or list of genes).
        Resolved through hitlist's ``peptide_mappings`` sidecar, which is
        built alongside the observations index.
    columns
        Project the parquet read to these columns.
    auto_build
        If True (default), trigger :func:`ensure_index_built` when the
        parquet is missing.
    drop_binding_assays
        If True (default), drop rows flagged as binding-assay data.

    Returns
    -------
    pd.DataFrame
    """
    from hitlist.observations import is_built, load_observations

    # Peptide-only reads do not depend on the mappings sidecar. Gene-filtered
    # reads do, so validate/migrate it even when observations already exist.
    if auto_build and (gene_name is not None or not is_built()):
        ensure_index_built()

    load_kwargs: dict = {
        "mhc_class": mhc_class,
        "species": mhc_species,
        "gene_name": gene_name,
        "columns": columns,
    }
    if peptides is not None:
        load_kwargs["peptide"] = (
            sorted(peptides) if isinstance(peptides, (set, frozenset)) else list(peptides)
        )

    df = load_observations(**load_kwargs)

    if drop_binding_assays and "is_binding_assay" in df.columns:
        df = df[~df["is_binding_assay"]]

    return df.reset_index(drop=True)
