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

    from perseus.iedb import scan_public_ms

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

# ── Column resolution ───────────────────────────────────────────────────────
#
# IEDB/CEDAR exports have two header rows: a category row and a field row.
# The field row contains pipe-separated names like "Epitope | Name".
# We parse the field header to discover column indices dynamically, falling
# back to known defaults if the header is missing or unparseable.

# Mapping from our internal key to the expected IEDB field-header substring.
# We match by substring so minor header changes (whitespace, case) are
# tolerated.
_COLUMN_NAMES: dict[str, list[str]] = {
    "assay_iri": ["Assay IRI", "Assay - Assay IRI"],
    "ref_iri": ["Reference IRI", "Reference - Reference IRI", "Reference | IEDB IRI"],
    "pmid": ["PMID", "Reference | PMID"],
    "ref_title": ["Title", "Reference | Title"],
    "epitope_name": ["Epitope | Name", "Epitope Name"],
    "source_organism": ["Source Organism", "Epitope | Source Organism"],
    "species": ["Epitope | Species", "Species"],
    "process_type": ["Process Type", "Host | Process Type"],
    "disease": ["Disease", "Host | Disease"],
    "source_tissue": ["Source Tissue", "Effector Cells | Source Tissue"],
    "cell_name": ["Cell Name", "Effector Cells | Cell Name"],
    "culture_condition": ["Culture Condition", "Assay | Culture Condition"],
    "mhc_restriction": ["MHC Restriction | Name", "MHC Restriction Name"],
    "mhc_class": ["MHC Allele Class", "Class"],
}

# Fallback hardcoded indices (IEDB mhc_ligand_full.csv as of 2024-2025).
_FALLBACK_INDICES: dict[str, int] = {
    "assay_iri": 0,
    "ref_iri": 1,
    "pmid": 3,
    "ref_title": 8,
    "epitope_name": 11,
    "source_organism": 23,
    "species": 25,
    "process_type": 50,
    "disease": 51,
    "source_tissue": 102,
    "cell_name": 104,
    "culture_condition": 106,
    "mhc_restriction": 107,
    "mhc_class": 111,
}


def _resolve_columns(category_header: list[str], field_header: list[str]) -> dict[str, int]:
    """Build column key -> index mapping from the IEDB field header row.

    Tries each candidate name for each key, checking both the field header
    and a combined "category | field" form.  Falls back to hardcoded indices
    for any key that cannot be resolved from the header.
    """
    indices: dict[str, int] = {}

    # Build combined header: "Category | Field" for each column
    combined = []
    for i in range(max(len(category_header), len(field_header))):
        cat = category_header[i].strip() if i < len(category_header) else ""
        fld = field_header[i].strip() if i < len(field_header) else ""
        if cat and fld:
            combined.append(f"{cat} | {fld}")
        else:
            combined.append(fld or cat)

    # Lowercase for matching
    combined_lower = [c.lower() for c in combined]
    field_lower = [f.strip().lower() for f in field_header]

    for key, candidates in _COLUMN_NAMES.items():
        found = False
        for candidate in candidates:
            cl = candidate.lower()
            # Try combined header first
            for i, header_val in enumerate(combined_lower):
                if cl in header_val or header_val.endswith(cl):
                    indices[key] = i
                    found = True
                    break
            if found:
                break
            # Try field header alone
            for i, header_val in enumerate(field_lower):
                if cl in header_val or header_val == cl:
                    indices[key] = i
                    found = True
                    break
            if found:
                break

    # Fill in any missing keys with fallback
    for key, fallback_idx in _FALLBACK_INDICES.items():
        if key not in indices:
            indices[key] = fallback_idx

    return indices


# ── Source context classification ───────────────────────────────────────────


# Tissue names in IEDB Source Tissue that map to reproductive or thymus.
_REPRODUCTIVE_TISSUES_IEDB = frozenset(
    {
        "testis",
        "ovary",
        "placenta",
        "cervix",
        "endometrium",
        "uterus",
        "fallopian tube",
        "prostate",
        "epididymis",
        "seminal vesicle",
        "vagina",
    }
)
_THYMUS_TISSUES_IEDB = frozenset({"thymus"})


def classify_ms_row(
    process_type: str,
    disease: str,
    culture_condition: str,
    source_tissue: str = "",
    cell_name: str = "",
) -> dict[str, bool | str]:
    """Classify a public-MS row into source-context flags.

    Uses structured IEDB/CEDAR columns to determine the biological context
    in which a peptide was observed.

    **Key rule**: all non-EBV cell lines are treated as cancer-derived, even
    when IEDB marks them "No immunization".  Many commercial cancer cell lines
    (HeLa, THP-1, A549, HCT 116, ...) appear in "healthy" studies but are
    cancer-derived.  Only direct ex vivo tissue from healthy donors qualifies
    as genuinely healthy.

    Parameters
    ----------
    process_type
        IEDB "Process Type" field.
    disease
        IEDB "Disease" field.
    culture_condition
        IEDB "Culture Condition" field.
    source_tissue
        IEDB "Source Tissue" field.
    cell_name
        IEDB "Cell Name" field (cell type or named cell line).

    Returns
    -------
    dict[str, bool | str]
        Flags:

        - ``src_cancer``: from a cancer sample OR a non-EBV cell line
        - ``src_healthy_tissue``: direct ex vivo, healthy, non-reproductive, non-thymus
        - ``src_healthy_thymus``: direct ex vivo, healthy, thymus tissue
        - ``src_healthy_reproductive``: direct ex vivo, healthy, reproductive tissue
        - ``src_cell_line``: cultured cell line (any kind)
        - ``src_ebv_lcl``: EBV-transformed B-LCL specifically
        - ``src_ex_vivo``: direct ex vivo (highest confidence)
        - ``cell_line_name``: named cell line (empty string if not a cell line)
    """
    process_type = str(process_type).strip() if pd.notna(process_type) else ""
    disease = str(disease).strip() if pd.notna(disease) else ""
    culture_condition = str(culture_condition).strip() if pd.notna(culture_condition) else ""
    source_tissue_lower = str(source_tissue).strip().lower() if pd.notna(source_tissue) else ""
    cell_name_str = str(cell_name).strip() if pd.notna(cell_name) else ""

    is_reproductive = source_tissue_lower in _REPRODUCTIVE_TISSUES_IEDB
    is_thymus = source_tissue_lower in _THYMUS_TISSUES_IEDB
    is_cell_line = culture_condition in (
        "Cell Line / Clone",
        "Cell Line / Clone (EBV transformed, B-LCL)",
    )
    is_ebv_lcl = culture_condition == "Cell Line / Clone (EBV transformed, B-LCL)"
    is_ex_vivo = culture_condition == "Direct Ex Vivo"

    # Key rule: non-EBV cell lines are cancer-derived even if marked "healthy"
    is_cancer = process_type == "Occurrence of cancer" or (is_cell_line and not is_ebv_lcl)

    # Healthy requires: direct ex vivo, no cancer, no disease
    is_healthy_donor = (
        is_ex_vivo and process_type == "No immunization" and disease in ("healthy", "")
    )

    # Extract cell line name only when it's actually a cell line
    cl_name = cell_name_str if (is_cell_line or is_ebv_lcl) else ""

    return {
        "src_cancer": is_cancer,
        "src_healthy_tissue": is_healthy_donor and not is_reproductive and not is_thymus,
        "src_healthy_thymus": is_healthy_donor and is_thymus,
        "src_healthy_reproductive": is_healthy_donor and is_reproductive,
        "src_cell_line": is_cell_line,
        "src_ebv_lcl": is_ebv_lcl,
        "src_ex_vivo": is_ex_vivo,
        "cell_line_name": cl_name,
    }


# ── Core scanning function ──────────────────────────────────────────────────


def _safe_col(row: list[str], idx: int) -> str:
    return row[idx] if len(row) > idx else ""


def _open_iedb_csv(path: Path) -> tuple[csv.reader, dict[str, int]]:
    """Open an IEDB/CEDAR CSV and return (reader, column_indices).

    Reads the two header rows, resolves column indices, and returns
    a reader positioned at the first data row.
    """
    fh = open(path, newline="")  # noqa: SIM115
    reader = csv.reader(fh)
    category_header = next(reader, [])
    field_header = next(reader, [])
    cols = _resolve_columns(category_header, field_header)
    return reader, cols


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
        reader, c = _open_iedb_csv(source_path)
        for row in reader:
            peptide = _safe_col(row, c["epitope_name"])
            if peptide not in selected:
                continue

            assay_iri = row[c["assay_iri"]] if row else ""
            if assay_iri in seen_assay_iris:
                continue
            seen_assay_iris.add(assay_iri)

            source_organism = _safe_col(row, c["source_organism"])
            species = _safe_col(row, c["species"])
            mhc_restriction = _safe_col(row, c["mhc_restriction"])

            if human_only and "Homo sapiens" not in (source_organism, species):
                continue
            if hla_only and not mhc_restriction.startswith("HLA-"):
                continue
            if mhc_class is not None:
                row_class = _safe_col(row, c["mhc_class"])
                if row_class != mhc_class:
                    continue

            raw_pmid = _safe_col(row, c["pmid"]).strip()
            pmid: str | int = ""
            if raw_pmid:
                try:
                    pmid = int(raw_pmid)
                except ValueError:
                    pmid = raw_pmid

            record: dict = {
                "reference_iri": _safe_col(row, c["ref_iri"]),
                "pmid": pmid,
                "reference_title": _safe_col(row, c["ref_title"]),
                "peptide": peptide,
                "source_organism": source_organism,
                "species": species,
                "mhc_restriction": mhc_restriction,
            }

            if classify_source:
                process_type = _safe_col(row, c["process_type"])
                disease = _safe_col(row, c["disease"])
                culture_condition = _safe_col(row, c["culture_condition"])
                record["process_type"] = process_type
                record["disease"] = disease
                source_tissue = _safe_col(row, c["source_tissue"])
                cell_name = _safe_col(row, c["cell_name"])
                record["source_tissue"] = source_tissue
                record["cell_name"] = cell_name
                record["culture_condition"] = culture_condition
                record.update(
                    classify_ms_row(
                        process_type, disease, culture_condition, source_tissue, cell_name
                    )
                )

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
        reader, c = _open_iedb_csv(source_path)
        for row in reader:
            assay_iri = row[c["assay_iri"]] if row else ""
            if assay_iri in seen_assay_iris:
                continue
            seen_assay_iris.add(assay_iri)

            source_organism = _safe_col(row, c["source_organism"])
            species = _safe_col(row, c["species"])
            if human_only and "Homo sapiens" not in (source_organism, species):
                continue

            process_type = _safe_col(row, c["process_type"])
            disease = _safe_col(row, c["disease"])
            culture_condition = _safe_col(row, c["culture_condition"])

            record = {
                "peptide": _safe_col(row, c["epitope_name"]),
                "mhc_restriction": _safe_col(row, c["mhc_restriction"]),
                "mhc_class": _safe_col(row, c["mhc_class"]),
                "process_type": process_type,
                "disease": disease,
                "source_tissue": (source_tissue := _safe_col(row, c["source_tissue"])),
                "cell_name": (cell_name := _safe_col(row, c["cell_name"])),
                "culture_condition": culture_condition,
                "source_organism": source_organism,
                "species": species,
            }
            record.update(
                classify_ms_row(process_type, disease, culture_condition, source_tissue, cell_name)
            )
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
