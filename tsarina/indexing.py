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

import sys
from pathlib import Path

import pandas as pd


def ensure_index_built(force: bool = False, verbose: bool = True) -> Path:
    """Ensure the hitlist observations parquet exists; build if not.

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
    from hitlist.observations import is_built, observations_path

    if force or not is_built():
        if verbose:
            print(
                "Building hitlist observations index (one-time ~2-5 min; cached afterwards)...",
                file=sys.stderr,
            )
        build_observations(force=force)
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

    Ensures the index exists (building if needed), then runs a pushdown-filtered
    parquet read for ``mhc_class`` / ``mhc_species`` (and ``gene_name`` if a
    flanking-built index is present), and finally filters to the supplied
    peptide set in memory.

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
        Only supported on flanking-built indices.  Use ``peptides`` instead
        when you already have a peptide list in hand.
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

    if auto_build and not is_built():
        ensure_index_built()

    # Push peptide set down to the parquet reader when possible — hitlist
    # 1.6.0+ accepts a ``peptide=`` filter.  Fall back to in-memory isin
    # for older installs.
    load_kwargs: dict = {
        "mhc_class": mhc_class,
        "species": mhc_species,
        "gene_name": gene_name,
        "columns": columns,
    }
    peptide_list: list[str] | None = None
    if peptides is not None:
        peptide_list = (
            sorted(peptides) if isinstance(peptides, (set, frozenset)) else list(peptides)
        )
        load_kwargs["peptide"] = peptide_list

    try:
        df = load_observations(**load_kwargs)
    except TypeError:
        # Older hitlist without ``peptide=`` pushdown.
        load_kwargs.pop("peptide", None)
        df = load_observations(**load_kwargs)
        if peptide_list is not None:
            df = df[df["peptide"].isin(set(peptide_list))]

    if drop_binding_assays and "is_binding_assay" in df.columns:
        df = df[~df["is_binding_assay"]]

    return df.reset_index(drop=True)
