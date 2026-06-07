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

"""Human Protein Atlas tissue expression enrichment for CTA gene evidence.

Ingests raw HPA data files (proteinatlas.tsv, normal_tissue.tsv,
rna_tissue_consensus.tsv) and computes tissue restriction columns per
gene: per-tissue nTPM maps/sets, protein support rankings, and
core-restriction flags.

The bundled CTA table's filter-driving columns (``rna_deflated_reproductive_frac``,
``passes_filters``) are produced by the regeneration scripts in ``scripts/``
(see ``add_cta_gene.py``), not here -- this module only enriches ad-hoc
evidence frames via :func:`tsarina.evidence.CTA_detailed_evidence`.

Typical usage::

    from tsarina.hpa import enrich_hpa_evidence

    df = enrich_hpa_evidence(
        gene_df=my_gene_dataframe,
        hpa_bulk_path="proteinatlas.tsv.zip",
        hpa_rna_path="rna_tissue_consensus.tsv.zip",
    )
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .tissues import (
    CORE_REPRODUCTIVE_TISSUES,
    EXTENDED_REPRODUCTIVE_TISSUES,
    HPA_MARKER_STRICT_MAX_RNA_TISSUES,
    PERMISSIVE_REPRODUCTIVE_TISSUES,
)

# Protein reliability ordering (best to worst)
PROTEIN_SUPPORT_ORDER: dict[str, int] = {
    "Enhanced": 0,
    "Supported": 1,
    "Approved": 2,
    "Uncertain": 3,
    "Missing": 4,
}


# ── Parsing utilities ───────────────────────────────────────────────────────


def parse_tissues(value) -> set[str]:
    """Parse semicolon-separated tissue list to lowercase set."""
    if pd.isna(value):
        return set()
    return {t.strip().lower() for t in str(value).split(";") if t.strip() and t.strip() != "nan"}


def parse_ntpm_entries(value) -> dict[str, float]:
    """Parse 'tissue:nTPM;tissue:nTPM' string to dict."""
    if pd.isna(value):
        return {}
    entries: dict[str, float] = {}
    for token in str(value).split(";"):
        token = token.strip()
        if not token or ":" not in token:
            continue
        tissue, amount = token.rsplit(":", 1)
        try:
            entries[tissue.strip().lower()] = float(amount.strip())
        except ValueError:
            continue
    return entries


def best_protein_support_label(has_protein: bool, reliability: str) -> str:
    """Return the best HPA protein reliability tier from a semicolon-separated string."""
    if not has_protein:
        return "Missing"
    values = [v.strip() for v in str(reliability).split(";") if v.strip() and v.strip() != "nan"]
    ranked = [v for v in values if v in PROTEIN_SUPPORT_ORDER and v != "Missing"]
    if not ranked:
        return "Missing"
    return min(ranked, key=lambda v: PROTEIN_SUPPORT_ORDER[v])


# ── Main enrichment function ────────────────────────────────────────────────


def enrich_hpa_evidence(
    gene_df: pd.DataFrame,
    hpa_bulk_path: str | Path | None = None,
    hpa_rna_path: str | Path | None = None,
    hpa_protein_path: str | Path | None = None,
) -> pd.DataFrame:
    """Enrich a gene DataFrame with HPA tissue restriction columns.

    Computes per-tissue RNA nTPM maps/sets, protein support rankings, and
    core-restriction flags from Human Protein Atlas data.  The deflated
    reproductive-fraction filter that drives the bundled CTA table lives in
    ``scripts/`` (see module docstring), not here.

    Parameters
    ----------
    gene_df
        DataFrame with at least ``gene_name`` and ``gene_id_stripped`` columns.
        Typically from CTA partition frames.
    hpa_bulk_path
        Path to HPA ``proteinatlas.tsv`` (bulk summary with RNA tissue
        specificity, tissue distribution, nTPM strings).
    hpa_rna_path
        Path to HPA ``rna_tissue_consensus.tsv`` (per-tissue nTPM values).
        Optional -- if not provided, uses the nTPM strings from bulk.
    hpa_protein_path
        Path to HPA ``normal_tissue.tsv`` (IHC protein expression).
        Optional.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with additional HPA evidence columns.
    """
    df = gene_df.copy()

    # ── Merge HPA bulk annotations ──────────────────────────────────────
    if hpa_bulk_path is not None:
        bulk_cols = [
            "Gene",
            "Gene synonym",
            "Ensembl",
            "Evidence",
            "HPA evidence",
            "RNA tissue specificity",
            "RNA tissue distribution",
            "RNA tissue specific nTPM",
            "Tissue expression cluster",
        ]
        bulk_df = pd.read_csv(str(hpa_bulk_path), sep="\t", usecols=bulk_cols, low_memory=False)
        bulk_df["gene_id_stripped"] = bulk_df["Ensembl"].astype(str).str.split(".").str[0]
        bulk_by_id = bulk_df.drop_duplicates(subset=["gene_id_stripped"])
        merge_cols = [
            "gene_id_stripped",
            "Evidence",
            "HPA evidence",
            "RNA tissue specificity",
            "RNA tissue distribution",
            "RNA tissue specific nTPM",
            "Tissue expression cluster",
        ]
        df = df.merge(bulk_by_id[merge_cols], on="gene_id_stripped", how="left")

    # ── Compute tissue restriction columns ──────────────────────────────

    # Parse tissue lists if available from CTA evidence columns
    if "protein_strict_expression" in df.columns:
        df["protein_tissue_set"] = df["protein_strict_expression"].map(parse_tissues)
    else:
        df["protein_tissue_set"] = [set()] * len(df)

    # Parse RNA nTPM map
    ntpm_col = "RNA tissue specific nTPM" if "RNA tissue specific nTPM" in df.columns else None
    if ntpm_col:
        df["rna_tissue_ntpm_map"] = df[ntpm_col].map(parse_ntpm_entries)
    else:
        df["rna_tissue_ntpm_map"] = [{}] * len(df)

    df["rna_tissue_set"] = df["rna_tissue_ntpm_map"].map(
        lambda m: {t for t, v in m.items() if v >= 1.0}
    )

    # Protein restriction checks
    df["protein_core_restricted"] = df["protein_tissue_set"].map(
        lambda ts: bool(ts) and ts.issubset(CORE_REPRODUCTIVE_TISSUES)
    )

    # RNA restriction checks
    df["rna_core_restricted"] = df["rna_tissue_set"].map(
        lambda ts: bool(ts) and ts.issubset(CORE_REPRODUCTIVE_TISSUES)
    )

    # Protein support ranking
    has_protein = df.get("protein_reproductive", pd.Series([False] * len(df)))
    reliability = df.get("protein_reliability", pd.Series([""] * len(df)))
    df["hpa_best_protein_support"] = [
        best_protein_support_label(bool(hp), str(rel))
        for hp, rel in zip(has_protein.fillna(False), reliability.fillna(""))
    ]
    df["hpa_protein_reliable"] = df["hpa_best_protein_support"].isin(["Enhanced", "Supported"])
    df["hpa_protein_reliability_enhanced"] = df["hpa_best_protein_support"].eq("Enhanced")

    # RNA tissue count
    df["rna_detected_tissue_count"] = df["rna_tissue_set"].map(len)
    df["rna_detected_tissues_le_7"] = df["rna_detected_tissue_count"].le(
        HPA_MARKER_STRICT_MAX_RNA_TISSUES
    )

    # Non-reproductive tissues
    df["nonreproductive_rna_tissues__core"] = df["rna_tissue_set"].map(
        lambda ts: ";".join(sorted(ts - CORE_REPRODUCTIVE_TISSUES))
    )
    df["nonreproductive_rna_tissues__extended"] = df["rna_tissue_set"].map(
        lambda ts: ";".join(sorted(ts - EXTENDED_REPRODUCTIVE_TISSUES))
    )

    # Composite restriction flag: protein-reliable, core-restricted, few RNA tissues.
    df["hpa_marker_strict_restricted"] = (
        df["protein_core_restricted"]
        & df["hpa_protein_reliability_enhanced"]
        & df["rna_detected_tissues_le_7"]
    )

    return df


# ── Per-tissue detail extraction ───────────────────────────────────────────


#: Tissues excluded from "somatic" max calculation.
_NON_SOMATIC: frozenset[str] = PERMISSIVE_REPRODUCTIVE_TISSUES | frozenset({"thymus"})


def extract_per_tissue_detail(ntpm_map: dict[str, float]) -> dict[str, object]:
    """Extract named per-tissue RNA values from an nTPM map.

    Parameters
    ----------
    ntpm_map
        Dict mapping tissue name (lowercase) to nTPM value, e.g. from
        ``parse_ntpm_entries()``.

    Returns
    -------
    dict
        Keys: ``rna_testis_ntpm``, ``rna_ovary_ntpm``, ``rna_placenta_ntpm``,
        ``rna_max_somatic_tissue``, ``rna_max_somatic_ntpm``,
        ``rna_somatic_detected_count``.
    """
    rna_testis = ntpm_map.get("testis", 0.0)
    rna_ovary = ntpm_map.get("ovary", 0.0)
    rna_placenta = ntpm_map.get("placenta", 0.0)

    somatic = {t: v for t, v in ntpm_map.items() if t not in _NON_SOMATIC}
    somatic_detected = {t: v for t, v in somatic.items() if v >= 1.0}

    if somatic:
        max_tissue = max(somatic, key=somatic.get)  # type: ignore[arg-type]
        max_ntpm = somatic[max_tissue]
    else:
        max_tissue = ""
        max_ntpm = 0.0

    return {
        "rna_testis_ntpm": rna_testis,
        "rna_ovary_ntpm": rna_ovary,
        "rna_placenta_ntpm": rna_placenta,
        "rna_max_somatic_tissue": max_tissue,
        "rna_max_somatic_ntpm": max_ntpm,
        "rna_somatic_detected_count": len(somatic_detected),
    }
