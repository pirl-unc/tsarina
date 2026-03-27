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

"""Scan IEDB and CEDAR MHC ligand exports for validated peptide-MHC observations.

Both IEDB (https://www.iedb.org/) and CEDAR (https://cedar.iedb.org/) provide
downloadable full MHC ligand assay exports.  This module provides utilities
to scan those exports for a set of peptides of interest, extract human-source
HLA-restricted entries, and deduplicate by assay IRI across both sources.

Typical usage::

    from ctabase.iedb import scan_public_ms

    hits = scan_public_ms(
        peptides={"SLYNTVATL", "GILGFVFTL"},
        iedb_path="mhc_ligand_full.csv",
        cedar_path="cedar-mhc-ligand-full.csv",
    )
"""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

# Column indices in the IEDB/CEDAR MHC ligand full export CSV.
# Row layout: assay IRI (0), reference IRI (1), PMID (3), title (8),
# epitope name (11), source organism (23), species (25),
# MHC restriction name (107).

_COL_ASSAY_IRI = 0
_COL_REF_IRI = 1
_COL_PMID = 3
_COL_REF_TITLE = 8
_COL_EPITOPE_NAME = 11
_COL_SOURCE_ORGANISM = 23
_COL_SPECIES = 25
_COL_MHC_RESTRICTION = 107


def scan_public_ms(
    peptides: set[str],
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    human_only: bool = True,
    hla_only: bool = True,
) -> pd.DataFrame:
    """Scan IEDB and/or CEDAR MHC ligand exports for matching peptides.

    Parameters
    ----------
    peptides
        Set of peptide sequences to search for (matched against the
        ``Epitope | Name`` column, index 11).
    iedb_path
        Path to the IEDB ``mhc_ligand_full.csv`` export.  Skipped if None.
    cedar_path
        Path to the CEDAR ``cedar-mhc-ligand-full.csv`` export.  Skipped if None.
    human_only
        If True (default), keep only rows where ``Epitope | Source Organism``
        or ``Epitope | Species`` is ``"Homo sapiens"``.
    hla_only
        If True (default), keep only rows where ``MHC Restriction | Name``
        starts with ``"HLA-"``.

    Returns
    -------
    pd.DataFrame
        Columns: ``reference_iri``, ``pmid``, ``reference_title``,
        ``peptide``, ``source_organism``, ``species``, ``mhc_restriction``.
        Deduplicated by assay IRI across both sources.
    """
    source_paths: list[Path] = []
    if iedb_path is not None:
        source_paths.append(Path(iedb_path))
    if cedar_path is not None:
        source_paths.append(Path(cedar_path))

    selected = set(peptides)
    rows: list[dict] = []
    seen_assay_iris: set[str] = set()

    for source_path in source_paths:
        if not source_path.exists():
            continue
        with open(source_path, newline="") as src:
            reader = csv.reader(src)
            next(reader)  # category header
            next(reader)  # field header
            for row in reader:
                peptide = row[_COL_EPITOPE_NAME] if len(row) > _COL_EPITOPE_NAME else ""
                if peptide not in selected:
                    continue

                assay_iri = row[_COL_ASSAY_IRI] if row else ""
                if assay_iri in seen_assay_iris:
                    continue
                seen_assay_iris.add(assay_iri)

                source_organism = (
                    row[_COL_SOURCE_ORGANISM] if len(row) > _COL_SOURCE_ORGANISM else ""
                )
                species = row[_COL_SPECIES] if len(row) > _COL_SPECIES else ""
                mhc_restriction = (
                    row[_COL_MHC_RESTRICTION] if len(row) > _COL_MHC_RESTRICTION else ""
                )

                if human_only and "Homo sapiens" not in (source_organism, species):
                    continue
                if hla_only and not mhc_restriction.startswith("HLA-"):
                    continue

                raw_pmid = row[_COL_PMID].strip() if len(row) > _COL_PMID else ""
                pmid: str | int = ""
                if raw_pmid:
                    try:
                        pmid = int(raw_pmid)
                    except ValueError:
                        pmid = raw_pmid

                rows.append(
                    {
                        "reference_iri": row[_COL_REF_IRI] if len(row) > _COL_REF_IRI else "",
                        "pmid": pmid,
                        "reference_title": row[_COL_REF_TITLE] if len(row) > _COL_REF_TITLE else "",
                        "peptide": peptide,
                        "source_organism": source_organism,
                        "species": species,
                        "mhc_restriction": mhc_restriction,
                    }
                )

    return pd.DataFrame(rows)


def peptide_ms_support(
    peptides: set[str],
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
) -> dict[str, set[str]]:
    """Return a mapping from peptide to the set of MHC restrictions observed.

    This is a convenience wrapper around :func:`scan_public_ms` that groups
    results by peptide.

    Parameters
    ----------
    peptides
        Set of peptide sequences to search for.
    iedb_path
        Path to the IEDB MHC ligand export.
    cedar_path
        Path to the CEDAR MHC ligand export.

    Returns
    -------
    dict[str, set[str]]
        Mapping from peptide sequence to set of MHC restriction names
        (e.g. ``{"SLYNTVATL": {"HLA-A*02:01", "HLA-B*08:01"}}``).
    """
    df = scan_public_ms(peptides, iedb_path=iedb_path, cedar_path=cedar_path)
    result: dict[str, set[str]] = {}
    for peptide, mhc in zip(df["peptide"], df["mhc_restriction"]):
        result.setdefault(peptide, set()).add(mhc)
    return result
