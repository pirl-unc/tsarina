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

Freshness is delegated to hitlist, which fingerprints its own curation inputs
and artifact contracts.  tsarina keeps no parallel notion of what a current
index looks like.
"""

from __future__ import annotations

import io
import sys
from contextlib import redirect_stdout
from pathlib import Path

import pandas as pd


def _artifact_state(path: Path) -> tuple[int, int] | None:
    """Size and nanosecond mtime of an artifact, or ``None`` when absent."""
    try:
        stat = path.stat()
    except OSError:
        return None
    return (stat.st_size, stat.st_mtime_ns)


def _artifact_states() -> tuple[tuple[int, int] | None, ...]:
    from hitlist.mappings import mappings_path
    from hitlist.observations import observations_path

    return tuple(_artifact_state(p) for p in (observations_path(), mappings_path()))


def _sources_registered() -> bool:
    """Whether hitlist can validate its artifacts against the raw exports.

    ``build_observations`` raises when no IEDB/CEDAR export is registered, so
    an index copied in without its sources can only be taken as it is found.
    """
    from .datasources import resolve_dataset_paths

    iedb_path, cedar_path = resolve_dataset_paths(require_iedb=False)
    return iedb_path is not None or cedar_path is not None


def ensure_index_built(force: bool = False, verbose: bool = True) -> Path:
    """Ensure the hitlist observations parquet and its mappings sidecar are current.

    Whether an artifact is stale is hitlist's judgement, not tsarina's:
    ``build_observations`` compares the stored observations artifact version
    against the fingerprints of the curation YAMLs (hitlist#429) and validates
    the mapping contract for the sidecar (hitlist#404), rebuilding either one
    that no longer matches.  Existence alone is not freshness — gating on it is
    how a legacy artifact keeps answering queries with retired curation.

    hitlist reports on stdout, which for ``tsarina hits`` carries the data, so
    its output goes to stderr instead: streamed while a known build runs, and
    otherwise captured and replayed only if an artifact actually changed.
    Announcing a stale rebuild before it starts needs a public validity
    predicate upstream (pirl-unc/hitlist#448).

    Parameters
    ----------
    force
        Rebuild even when both artifacts are current.
    verbose
        Report on stderr when a build is actually performed.  Silent when the
        cached artifacts are reused unchanged.

    Returns
    -------
    Path
        Path to ``observations.parquet``.
    """
    from hitlist.builder import build_observations
    from hitlist.observations import is_built, observations_path

    observations = observations_path()
    building = force or not is_built()
    if not building and not _sources_registered():
        return observations

    if building:
        if verbose:
            print(
                "Building hitlist observations index (one-time ~2-5 min; cached afterwards)...",
                file=sys.stderr,
            )
        with redirect_stdout(sys.stderr):
            build_observations(force=force)
        return observations

    before = _artifact_states()
    report = io.StringIO()
    with redirect_stdout(report):
        build_observations(force=False)
    if verbose and _artifact_states() != before:
        print(
            "hitlist artifacts were stale (curation or artifact version changed) "
            "and have been rebuilt:",
            file=sys.stderr,
        )
        print(report.getvalue().rstrip(), file=sys.stderr)
    return observations


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
        If True (default), route through :func:`ensure_index_built` so a
        missing or stale index is built before the read.
    drop_binding_assays
        If True (default), drop rows flagged as binding-assay data.

    Returns
    -------
    pd.DataFrame
    """
    from hitlist.observations import load_observations

    # One policy for every read: hitlist decides whether either artifact needs
    # rebuilding, so a peptide-only read no longer skips that check either.
    if auto_build:
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
