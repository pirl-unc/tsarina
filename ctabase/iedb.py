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

    # With disease/tissue context classification
    hits = scan_public_ms(
        peptides=my_peptides,
        iedb_path="mhc_ligand_full.csv",
        classify_source=True,   # adds src_cancer, src_healthy, etc.
        mhc_class="I",          # MHC class I only
    )
"""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

# Column indices in the IEDB/CEDAR MHC ligand full export CSV.

_COL_ASSAY_IRI = 0
_COL_REF_IRI = 1
_COL_PMID = 3
_COL_REF_TITLE = 8
_COL_EPITOPE_NAME = 11
_COL_SOURCE_ORGANISM = 23
_COL_SPECIES = 25
_COL_PROCESS_TYPE = 50
_COL_DISEASE = 51
_COL_SOURCE_TISSUE = 102
_COL_CELL_NAME = 104
_COL_CULTURE_CONDITION = 106
_COL_MHC_RESTRICTION = 107
_COL_MHC_CLASS = 111


# ── Source context classification ───────────────────────────────────────────


def classify_ms_row(
    process_type: str,
    disease: str,
    culture_condition: str,
) -> dict[str, bool]:
    """Classify a public-MS row into source-context flags.

    Uses structured IEDB/CEDAR columns (Process Type, Disease,
    Culture Condition) to determine the biological context in which
    a peptide was observed.

    Parameters
    ----------
    process_type
        IEDB "Process Type" field (e.g. "Occurrence of cancer",
        "No immunization").
    disease
        IEDB "Disease" field (e.g. "healthy", "melanoma").
    culture_condition
        IEDB "Culture Condition" field (e.g. "Cell Line / Clone",
        "Direct Ex Vivo").

    Returns
    -------
    dict[str, bool]
        Flags: ``src_cancer``, ``src_healthy``, ``src_cell_line``,
        ``src_ex_vivo``, ``src_ebv_lcl``.
    """
    process_type = str(process_type).strip() if pd.notna(process_type) else ""
    disease = str(disease).strip() if pd.notna(disease) else ""
    culture_condition = str(culture_condition).strip() if pd.notna(culture_condition) else ""

    return {
        "src_cancer": process_type == "Occurrence of cancer",
        "src_healthy": process_type == "No immunization" and disease in ("healthy", ""),
        "src_cell_line": culture_condition
        in (
            "Cell Line / Clone",
            "Cell Line / Clone (EBV transformed, B-LCL)",
        ),
        "src_ex_vivo": culture_condition == "Direct Ex Vivo",
        "src_ebv_lcl": culture_condition == "Cell Line / Clone (EBV transformed, B-LCL)",
    }


# ── Core scanning function ──────────────────────────────────────────────────


def _safe_col(row: list[str], idx: int) -> str:
    return row[idx] if len(row) > idx else ""


def scan_public_ms(
    peptides: set[str],
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    human_only: bool = True,
    hla_only: bool = True,
    mhc_class: str | None = None,
    classify_source: bool = False,
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
    mhc_class
        If set (e.g. ``"I"`` or ``"II"``), filter to rows matching that
        MHC class (column 111).  None means no class filtering.
    classify_source
        If True, add source-context columns (``src_cancer``,
        ``src_healthy``, ``src_cell_line``, ``src_ex_vivo``,
        ``src_ebv_lcl``) and additional columns (``disease``,
        ``source_tissue``, ``cell_name``).

    Returns
    -------
    pd.DataFrame
        Deduplicated by assay IRI across both sources.
        Base columns: ``reference_iri``, ``pmid``, ``reference_title``,
        ``peptide``, ``source_organism``, ``species``, ``mhc_restriction``.
        With ``classify_source=True``: also ``process_type``, ``disease``,
        ``source_tissue``, ``cell_name``, ``culture_condition``,
        ``src_cancer``, ``src_healthy``, ``src_cell_line``,
        ``src_ex_vivo``, ``src_ebv_lcl``.
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
                peptide = _safe_col(row, _COL_EPITOPE_NAME)
                if peptide not in selected:
                    continue

                assay_iri = row[_COL_ASSAY_IRI] if row else ""
                if assay_iri in seen_assay_iris:
                    continue
                seen_assay_iris.add(assay_iri)

                source_organism = _safe_col(row, _COL_SOURCE_ORGANISM)
                species = _safe_col(row, _COL_SPECIES)
                mhc_restriction = _safe_col(row, _COL_MHC_RESTRICTION)

                if human_only and "Homo sapiens" not in (source_organism, species):
                    continue
                if hla_only and not mhc_restriction.startswith("HLA-"):
                    continue
                if mhc_class is not None:
                    row_class = _safe_col(row, _COL_MHC_CLASS)
                    if row_class != mhc_class:
                        continue

                raw_pmid = _safe_col(row, _COL_PMID).strip()
                pmid: str | int = ""
                if raw_pmid:
                    try:
                        pmid = int(raw_pmid)
                    except ValueError:
                        pmid = raw_pmid

                record: dict = {
                    "reference_iri": _safe_col(row, _COL_REF_IRI),
                    "pmid": pmid,
                    "reference_title": _safe_col(row, _COL_REF_TITLE),
                    "peptide": peptide,
                    "source_organism": source_organism,
                    "species": species,
                    "mhc_restriction": mhc_restriction,
                }

                if classify_source:
                    process_type = _safe_col(row, _COL_PROCESS_TYPE)
                    disease = _safe_col(row, _COL_DISEASE)
                    culture_condition = _safe_col(row, _COL_CULTURE_CONDITION)
                    record["process_type"] = process_type
                    record["disease"] = disease
                    record["source_tissue"] = _safe_col(row, _COL_SOURCE_TISSUE)
                    record["cell_name"] = _safe_col(row, _COL_CELL_NAME)
                    record["culture_condition"] = culture_condition
                    record.update(classify_ms_row(process_type, disease, culture_condition))

                rows.append(record)

    return pd.DataFrame(rows)


def profile_dataset(
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    human_only: bool = True,
) -> pd.DataFrame:
    """Profile an entire IEDB/CEDAR MHC ligand export without peptide filtering.

    Scans every row and extracts key fields for summary statistics:
    MHC class distribution, tissue coverage, cancer type breakdown, etc.

    Parameters
    ----------
    iedb_path
        Path to the IEDB ``mhc_ligand_full.csv`` export.
    cedar_path
        Path to the CEDAR ``cedar-mhc-ligand-full.csv`` export.
    human_only
        If True (default), keep only human-source rows.

    Returns
    -------
    pd.DataFrame
        One row per assay (deduplicated by assay IRI).
        Columns: ``peptide``, ``mhc_restriction``, ``mhc_class``,
        ``process_type``, ``disease``, ``source_tissue``, ``cell_name``,
        ``culture_condition``, ``source_organism``, ``species``,
        ``src_cancer``, ``src_healthy``, ``src_cell_line``,
        ``src_ex_vivo``, ``src_ebv_lcl``.
    """
    source_paths: list[Path] = []
    if iedb_path is not None:
        source_paths.append(Path(iedb_path))
    if cedar_path is not None:
        source_paths.append(Path(cedar_path))

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
                assay_iri = row[_COL_ASSAY_IRI] if row else ""
                if assay_iri in seen_assay_iris:
                    continue
                seen_assay_iris.add(assay_iri)

                source_organism = _safe_col(row, _COL_SOURCE_ORGANISM)
                species = _safe_col(row, _COL_SPECIES)
                if human_only and "Homo sapiens" not in (source_organism, species):
                    continue

                process_type = _safe_col(row, _COL_PROCESS_TYPE)
                disease = _safe_col(row, _COL_DISEASE)
                culture_condition = _safe_col(row, _COL_CULTURE_CONDITION)

                record = {
                    "peptide": _safe_col(row, _COL_EPITOPE_NAME),
                    "mhc_restriction": _safe_col(row, _COL_MHC_RESTRICTION),
                    "mhc_class": _safe_col(row, _COL_MHC_CLASS),
                    "process_type": process_type,
                    "disease": disease,
                    "source_tissue": _safe_col(row, _COL_SOURCE_TISSUE),
                    "cell_name": _safe_col(row, _COL_CELL_NAME),
                    "culture_condition": culture_condition,
                    "source_organism": source_organism,
                    "species": species,
                }
                record.update(classify_ms_row(process_type, disease, culture_condition))
                rows.append(record)

    return pd.DataFrame(rows)


def peptide_ms_support(
    peptides: set[str],
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    mhc_class: str | None = None,
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
    mhc_class
        If set, filter to this MHC class (e.g. ``"I"``).

    Returns
    -------
    dict[str, set[str]]
        Mapping from peptide sequence to set of MHC restriction names
        (e.g. ``{"SLYNTVATL": {"HLA-A*02:01", "HLA-B*08:01"}}``).
    """
    df = scan_public_ms(peptides, iedb_path=iedb_path, cedar_path=cedar_path, mhc_class=mhc_class)
    result: dict[str, set[str]] = {}
    for peptide, mhc in zip(df["peptide"], df["mhc_restriction"]):
        result.setdefault(peptide, set()).add(mhc)
    return result
