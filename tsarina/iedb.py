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

"""Thin adapter over :func:`hitlist.scanner.scan`.

For default queries, prefer :func:`tsarina.indexing.load_ms_evidence` — it
reads the cached hitlist observations parquet and is sub-second after the
first build.  The functions here stay around for callers that need to point
at a *non-registered* IEDB/CEDAR snapshot (e.g. a newer export than the one
in the hitlist registry, or an ad-hoc subset for testing).

Historically this module reimplemented hitlist's scanner loop.  It no longer
does — every real behavior (two-header parsing, column resolution, species
filter via mhcgnomes, source classification, binding-assay detection,
allele-resolution filtering, assay-IRI dedup) lives in
``hitlist.scanner.scan``.  These wrappers exist to (a) keep a stable tsarina
import path and (b) translate tsarina's historical ``mhc_species`` argument
onto hitlist's ``human_only`` / ``mhc_species`` kwargs.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from hitlist.curation import normalize_species
from hitlist.scanner import scan as _hitlist_scan


def _species_kwargs(mhc_species: str | None) -> dict:
    """Translate tsarina's ``mhc_species`` argument onto hitlist scanner kwargs.

    For ``"Homo sapiens"`` (tsarina's historical default) use ``human_only=True``,
    which applies the host/species fallback when an MHC restriction string
    cannot be resolved via mhcgnomes (e.g. bare ``"HLA class I"``).  For any
    other species, use the strict ``mhc_species`` path.  For ``None``, disable
    species filtering entirely.
    """
    if mhc_species is None:
        return {"human_only": False}
    if normalize_species(mhc_species) == "Homo sapiens":
        return {"human_only": True}
    return {"mhc_species": mhc_species, "human_only": False}


def scan_public_ms(
    peptides: set[str],
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    mhc_species: str | None = "Homo sapiens",
    mhc_class: str | None = None,
    classify_source: bool = False,
    min_allele_resolution: str | None = None,
) -> pd.DataFrame:
    """Scan IEDB/CEDAR MHC ligand exports for matching peptides.

    Thin wrapper over :func:`hitlist.scanner.scan`.  Prefer
    :func:`tsarina.indexing.load_ms_evidence` when querying against the
    registered snapshot — it hits a cached parquet and is much faster.

    Parameters
    ----------
    peptides
        Peptide sequences to match against ``Epitope | Name``.
    iedb_path, cedar_path
        Paths to IEDB / CEDAR exports.  Missing paths are silently skipped;
        if both are None the function returns an empty DataFrame.
    mhc_species
        Species whose MHC molecules to keep.  Default ``"Homo sapiens"``;
        ``None`` disables filtering.
    mhc_class
        ``"I"`` or ``"II"`` to filter; ``None`` keeps both.
    classify_source
        Emit ``src_*`` source-context flags (cancer / healthy tissue / cell
        line / ex vivo / EBV-LCL) via ``hitlist.curation.classify_ms_row``.
    min_allele_resolution
        One of ``"four_digit"``, ``"two_digit"``, ``"serological"``,
        ``"class_only"``.  Rows with coarser resolution are dropped.

    Returns
    -------
    pd.DataFrame
        Deduplicated by assay IRI across both sources.
    """
    if iedb_path is None and cedar_path is None:
        return pd.DataFrame()
    return _hitlist_scan(
        peptides=peptides,
        iedb_path=iedb_path,
        cedar_path=cedar_path,
        mhc_class=mhc_class,
        classify_source=classify_source,
        min_allele_resolution=min_allele_resolution,
        **_species_kwargs(mhc_species),
    )


def profile_dataset(
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    mhc_species: str | None = "Homo sapiens",
) -> pd.DataFrame:
    """Profile an entire IEDB/CEDAR export with no peptide filter.

    Equivalent to :func:`scan_public_ms` with ``peptides=None``.  Always
    includes source classification columns (``src_cancer``, ``src_healthy_*``,
    ``src_cell_line``, ``src_ebv_lcl``, ``src_ex_vivo``), which
    :mod:`tsarina.negatives` consumes.
    """
    if iedb_path is None and cedar_path is None:
        return pd.DataFrame()
    return _hitlist_scan(
        peptides=None,
        iedb_path=iedb_path,
        cedar_path=cedar_path,
        classify_source=True,
        **_species_kwargs(mhc_species),
    )
