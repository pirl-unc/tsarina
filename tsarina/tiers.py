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

"""Per-modality restriction assessment with cross-source synthesis.

Each data modality (protein IHC, RNA, MS) provides its own restriction
assessment and metadata.  A synthesis layer integrates them into a
unified ``restriction`` with ``restriction_confidence``.

**Protein** (IHC):

- ``protein_restriction``: TESTIS / PLACENTAL / REPRODUCTIVE / SOMATIC / empty.
  Based on which tissues have IHC-detected protein.

**RNA** (HPA consensus):

- ``rna_restriction``: TESTIS / PLACENTAL / REPRODUCTIVE / empty.
  Based on which reproductive tissues have nTPM >= 1.
- ``rna_restriction_level``: STRICT / MODERATE / PERMISSIVE / empty.
  Based on deflated reproductive fraction.

**MS** (IEDB/CEDAR via hitlist, runtime):

- ``ms_restriction``: CANCER_ONLY / EXPECTED_TISSUE / SINGLETON_HEALTHY /
  RECURRENT_HEALTHY / NO_MS_DATA.

**Synthesized**:

- ``restriction``: best available tissue category from protein > RNA.
- ``restriction_confidence``: HIGH / MODERATE / LOW based on modality
  agreement and quality.
"""

from __future__ import annotations

import pandas as pd

from .tissues import PERMISSIVE_REPRODUCTIVE_TISSUES

# ── Constants ──────────────────────────────────────────────────────────────

#: Tissue restriction categories (shared across protein, RNA, synthesis).
RESTRICTION_VALUES: list[str] = ["TESTIS", "PLACENTAL", "REPRODUCTIVE", "SOMATIC"]

#: RNA restriction quality levels, ordered strictest → loosest.
RNA_RESTRICTION_LEVELS: list[str] = ["STRICT", "MODERATE", "PERMISSIVE", "LEAKY"]

#: MS restriction classifications.
MS_RESTRICTION_VALUES: list[str] = [
    "CANCER_ONLY",
    "EXPECTED_TISSUE",
    "SINGLETON_HEALTHY",
    "RECURRENT_HEALTHY",
    "NO_MS_DATA",
]

#: Restriction confidence levels.
CONFIDENCE_VALUES: list[str] = ["HIGH", "MODERATE", "LOW"]

#: Tissues excluded from "somatic" calculations.
_NON_SOMATIC: frozenset[str] = PERMISSIVE_REPRODUCTIVE_TISSUES | frozenset({"thymus"})


# ── Rank functions ─────────────────────────────────────────────────────────

_RESTRICTION_RANK = {v: i for i, v in enumerate(RESTRICTION_VALUES)}
_MS_RANK = {v: i for i, v in enumerate(MS_RESTRICTION_VALUES)}
_CONFIDENCE_RANK = {v: i for i, v in enumerate(CONFIDENCE_VALUES)}


def restriction_rank(value: str | None) -> int:
    """Integer rank for a restriction value (lower = more restrictive)."""
    if not value:
        return len(RESTRICTION_VALUES)
    return _RESTRICTION_RANK.get(value, len(RESTRICTION_VALUES))


def ms_restriction_rank(value: str | None) -> int:
    """Integer rank for an MS restriction value (lower = safer)."""
    if not value:
        return len(MS_RESTRICTION_VALUES)
    return _MS_RANK.get(value, len(MS_RESTRICTION_VALUES))


def confidence_rank(value: str | None) -> int:
    """Integer rank for confidence (lower = higher confidence)."""
    if not value:
        return len(CONFIDENCE_VALUES)
    return _CONFIDENCE_RANK.get(value, len(CONFIDENCE_VALUES))


# ── Protein restriction ───────────────────────────────────────────────────


def _parse_protein_tissues(protein_strict_expression: str) -> set[str]:
    val = str(protein_strict_expression).strip().lower()
    if val in ("no data", "not detected", "", "nan"):
        return set()
    return {t.strip() for t in val.split(";") if t.strip()}


#: All tissues considered reproductive (core + extended + breast).
#: HPA uses "endometrium 1"/"endometrium 2" so include those.
_ALL_REPRODUCTIVE: frozenset[str] = _NON_SOMATIC | frozenset({"endometrium 1", "endometrium 2"})


def assign_protein_restriction(row: pd.Series) -> str:
    """Assign tissue restriction from IHC protein data.

    Returns one of: TESTIS, PLACENTAL, REPRODUCTIVE, SOMATIC, or empty.
    Each single-tissue value means *only* that core tissue detected.
    REPRODUCTIVE means multiple reproductive tissues.
    SOMATIC means non-reproductive tissue detected.
    """
    tissues = _parse_protein_tissues(str(row.get("protein_strict_expression", "")))
    if not tissues:
        return "NO_DATA"
    # Remove thymus (expected for CTAs)
    tissues = tissues - {"thymus"}
    if not tissues:
        return "NO_DATA"
    # Check for somatic (non-reproductive) tissues
    non_repro = tissues - _ALL_REPRODUCTIVE
    if non_repro:
        return "SOMATIC"
    # Which core reproductive tissues are present?
    core = tissues & {"testis", "ovary", "placenta"}
    if core == {"testis"}:
        return "TESTIS"
    if core == {"placenta"}:
        return "PLACENTAL"
    if core:
        return "REPRODUCTIVE"
    # Only extended reproductive tissues (epididymis, etc.), no core
    return "REPRODUCTIVE"


def _protein_tissue_flag(protein_strict_expression: str, tissue: str) -> str:
    val = str(protein_strict_expression).strip().lower()
    if val in ("no data", "not detected", "", "nan"):
        return ""
    tissues = {t.strip() for t in val.split(";") if t.strip()}
    return str(tissue in tissues)


# ── RNA restriction ───────────────────────────────────────────────────────


def assign_rna_restriction(row: pd.Series) -> str:
    """Assign tissue restriction from RNA per-tissue nTPM data.

    Requires ``rna_testis_ntpm``, ``rna_ovary_ntpm``, ``rna_placenta_ntpm``,
    ``rna_somatic_detected_count`` columns (from HPA rna_tissue_consensus).

    Always assigns based on RNA data, regardless of filter status.
    """
    testis = float(row.get("rna_testis_ntpm", 0) or 0)
    ovary = float(row.get("rna_ovary_ntpm", 0) or 0)
    placenta = float(row.get("rna_placenta_ntpm", 0) or 0)
    somatic_count = int(row.get("rna_somatic_detected_count", 0) or 0)

    has_testis = testis >= 1.0
    has_ovary = ovary >= 1.0
    has_placenta = placenta >= 1.0
    has_somatic = somatic_count > 0

    if has_somatic:
        return "SOMATIC"

    if not (has_testis or has_ovary or has_placenta) and not has_somatic:
        # Nothing detected anywhere at >= 1 nTPM
        return "NO_DATA"

    if has_testis and not has_ovary and not has_placenta:
        return "TESTIS"
    if has_placenta and not has_ovary and not has_testis:
        return "PLACENTAL"
    return "REPRODUCTIVE"


def assign_rna_restriction_level(row: pd.Series) -> str:
    """Assign RNA restriction quality level from deflated fraction.

    Uses ``rna_somatic_detected_count`` (from HPA per-tissue data) rather
    than the CSV's ``rna_reproductive`` to ensure consistency with
    ``rna_restriction``.
    """
    try:
        frac = float(row.get("rna_deflated_reproductive_frac", -1))
    except (ValueError, TypeError):
        frac = -1.0

    if frac < 0:
        return "NO_DATA"

    somatic_count = int(row.get("rna_somatic_detected_count", 0) or 0)
    rna_is_reproductive = somatic_count == 0

    if rna_is_reproductive and frac >= 0.99:
        return "STRICT"
    if frac >= 0.95:
        return "MODERATE"
    if frac >= 0.80:
        return "PERMISSIVE"
    return "LEAKY"


# ── Synthesized restriction ───────────────────────────────────────────────


def synthesize_restriction(row: pd.Series) -> tuple[str, str]:
    """Synthesize restriction + confidence from all modalities.

    Returns (restriction, restriction_confidence).
    """
    protein_r = str(row.get("protein_restriction", "") or "")
    protein_rel = str(row.get("protein_reliability", "") or "")
    rna_r = str(row.get("rna_restriction", "") or "")
    rna_level = str(row.get("rna_restriction_level", "") or "")
    ms_r = str(row.get("ms_restriction", "") or "")

    # Best available tissue category: protein > RNA
    # Skip NO_DATA values when choosing
    protein_has_data = protein_r and protein_r != "NO_DATA"
    rna_has_data = rna_r and rna_r != "NO_DATA"

    if protein_has_data:
        tissue = protein_r
    elif rna_has_data:
        tissue = rna_r
    else:
        tissue = "NO_DATA"

    # Confidence scoring
    score = 0.0
    sources = 0

    if protein_has_data:
        sources += 1
        score += 1.0
        if protein_rel in ("Enhanced", "Supported"):
            score += 0.5

    if rna_has_data:
        sources += 1
        rna_agrees = rna_r == tissue or (rna_r == "REPRODUCTIVE" and protein_has_data)
        if rna_agrees:
            score += 1.0
            if rna_level == "STRICT":
                score += 0.5

    if ms_r and ms_r not in ("NO_MS_DATA", "NO_DATA", ""):
        sources += 1
        if ms_r in ("CANCER_ONLY", "EXPECTED_TISSUE"):
            score += 1.0
        elif ms_r == "SINGLETON_HEALTHY":
            score += 0.5
        # RECURRENT_HEALTHY contributes 0

    if sources == 0:
        confidence = "NO_DATA"
    elif score / sources >= 1.2:
        confidence = "HIGH"
    elif score / sources >= 0.8:
        confidence = "MODERATE"
    else:
        confidence = "LOW"

    return (tissue, confidence)


# ── Main assignment function ──────────────────────────────────────────────


def assign_all_axes(df: pd.DataFrame) -> pd.DataFrame:
    """Add all per-modality and synthesized restriction columns.

    Parameters
    ----------
    df
        CTA evidence DataFrame. Must include per-tissue RNA nTPM columns
        (``rna_testis_ntpm``, ``rna_ovary_ntpm``, ``rna_placenta_ntpm``)
        for RNA restriction assignment.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with protein_restriction, rna_restriction,
        rna_restriction_level, ms_restriction, restriction,
        restriction_confidence, protein_testis/ovary/placenta columns.
    """
    out = df.copy()

    # Protein
    out["protein_restriction"] = out.apply(assign_protein_restriction, axis=1)
    pse = out.get("protein_strict_expression", pd.Series([""] * len(out)))
    out["protein_testis"] = pse.map(lambda v: _protein_tissue_flag(v, "testis"))
    out["protein_ovary"] = pse.map(lambda v: _protein_tissue_flag(v, "ovary"))
    out["protein_placenta"] = pse.map(lambda v: _protein_tissue_flag(v, "placenta"))

    # RNA
    out["rna_restriction"] = out.apply(assign_rna_restriction, axis=1)
    out["rna_restriction_level"] = out.apply(assign_rna_restriction_level, axis=1)

    # MS (defaults; populated at runtime when IEDB data available)
    if "ms_restriction" not in out.columns:
        out["ms_restriction"] = "NO_MS_DATA"

    # Synthesis
    synth = out.apply(synthesize_restriction, axis=1, result_type="expand")
    out["restriction"] = synth[0]
    out["restriction_confidence"] = synth[1]

    return out


# ── RNA per-tissue enrichment from HPA ────────────────────────────────────


def enrich_rna_per_tissue(
    df: pd.DataFrame,
    rna_tissue_path: str,
) -> pd.DataFrame:
    """Add per-tissue RNA nTPM columns from HPA rna_tissue_consensus.tsv.

    Parameters
    ----------
    df
        CTA evidence DataFrame with ``Ensembl_Gene_ID`` column.
    rna_tissue_path
        Path to HPA ``rna_tissue_consensus.tsv``.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with rna_testis_ntpm, rna_ovary_ntpm,
        rna_placenta_ntpm, rna_max_somatic_tissue, rna_max_somatic_ntpm,
        rna_somatic_detected_count columns added.
    """
    rna = pd.read_csv(rna_tissue_path, sep="\t")
    rna["tissue_lower"] = rna["Tissue"].str.strip().str.lower()

    non_somatic = {t.lower() for t in _NON_SOMATIC} | {"thymus"}

    out = df.copy()
    testis_vals = []
    ovary_vals = []
    placenta_vals = []
    max_somatic_tissues = []
    max_somatic_ntpms = []
    somatic_detected_counts = []

    for _, row in out.iterrows():
        gene_id = row.get("Ensembl_Gene_ID", "")
        gene_rna = rna[rna["Gene"] == gene_id]

        if gene_rna.empty:
            testis_vals.append(0.0)
            ovary_vals.append(0.0)
            placenta_vals.append(0.0)
            max_somatic_tissues.append("")
            max_somatic_ntpms.append(0.0)
            somatic_detected_counts.append(0)
            continue

        tissue_ntpm = dict(zip(gene_rna["tissue_lower"], gene_rna["nTPM"]))
        testis_vals.append(tissue_ntpm.get("testis", 0.0))
        ovary_vals.append(tissue_ntpm.get("ovary", 0.0))
        placenta_vals.append(tissue_ntpm.get("placenta", 0.0))

        somatic = {t: v for t, v in tissue_ntpm.items() if t not in non_somatic}
        somatic_detected = {t: v for t, v in somatic.items() if v >= 1.0}

        if somatic:
            max_t = max(somatic, key=somatic.get)  # type: ignore[arg-type]
            max_somatic_tissues.append(max_t)
            max_somatic_ntpms.append(somatic[max_t])
        else:
            max_somatic_tissues.append("")
            max_somatic_ntpms.append(0.0)

        somatic_detected_counts.append(len(somatic_detected))

    out["rna_testis_ntpm"] = testis_vals
    out["rna_ovary_ntpm"] = ovary_vals
    out["rna_placenta_ntpm"] = placenta_vals
    out["rna_max_somatic_tissue"] = max_somatic_tissues
    out["rna_max_somatic_ntpm"] = max_somatic_ntpms
    out["rna_somatic_detected_count"] = somatic_detected_counts

    return out


# ── MS safety aggregation ─────────────────────────────────────────────────


def aggregate_gene_ms_safety(
    classified_hits: pd.DataFrame,
    peptide_gene_map: pd.DataFrame,
) -> pd.DataFrame:
    """Aggregate per-peptide MS evidence to per-gene MS restriction.

    Parameters
    ----------
    classified_hits
        Raw classified MS hits with columns: ``peptide``,
        ``src_cancer``, ``src_healthy_tissue``, ``src_healthy_reproductive``,
        ``src_healthy_thymus``, ``source_tissue``.
    peptide_gene_map
        Peptide-to-gene mapping with columns: ``peptide``, ``gene_name``.

    Returns
    -------
    pd.DataFrame
        One row per gene with ms_restriction and metadata columns.
    """
    _MS_COLS = [
        "gene_name",
        "ms_restriction",
        "ms_peptide_count",
        "ms_cancer_peptide_count",
        "ms_ebv_lcl_peptide_count",
        "ms_healthy_somatic_peptide_count",
        "ms_healthy_somatic_tissue_count",
        "ms_healthy_somatic_tissues",
        "ms_healthy_reproductive_peptide_count",
        "ms_healthy_thymus_peptide_count",
        "ms_pmids",
    ]

    if classified_hits.empty or peptide_gene_map.empty:
        return pd.DataFrame(columns=_MS_COLS)

    merged = classified_hits.merge(peptide_gene_map[["peptide", "gene_name"]], on="peptide")
    if merged.empty:
        return pd.DataFrame(columns=_MS_COLS)

    # Per-gene, per-peptide: aggregate source flags
    agg_dict = {
        "has_cancer": ("src_cancer", "any"),
        "has_ebv_lcl": ("src_ebv_lcl", "any") if "src_ebv_lcl" in merged.columns else None,
        "has_healthy_somatic": ("src_healthy_tissue", "any"),
        "has_healthy_reproductive": ("src_healthy_reproductive", "any"),
        "has_healthy_thymus": ("src_healthy_thymus", "any"),
    }
    agg_dict = {k: v for k, v in agg_dict.items() if v is not None}
    pep_gene = merged.groupby(["gene_name", "peptide"]).agg(**agg_dict).reset_index()
    if "has_ebv_lcl" not in pep_gene.columns:
        pep_gene["has_ebv_lcl"] = False

    # Healthy somatic tissues per gene
    healthy_somatic_rows = merged[merged["src_healthy_tissue"].astype(bool)].copy()
    if "source_tissue" in healthy_somatic_rows.columns and not healthy_somatic_rows.empty:
        somatic_tissues_per_gene = healthy_somatic_rows.groupby("gene_name")["source_tissue"].agg(
            ms_healthy_somatic_tissue_count=lambda s: s.str.strip().str.lower().nunique(),
            ms_healthy_somatic_tissues=lambda s: ";".join(
                sorted(s.str.strip().str.lower().unique())
            ),
        )
    else:
        somatic_tissues_per_gene = None

    # PMIDs per gene
    if "pmid" in merged.columns and not merged["pmid"].isna().all():
        pmids_per_gene = (
            merged.groupby("gene_name")["pmid"]
            .apply(lambda s: ";".join(sorted(s.dropna().astype(str).unique())))
            .rename("ms_pmids")
            .reset_index()
        )
    else:
        pmids_per_gene = None

    # Per-gene aggregation
    gene_agg = (
        pep_gene.groupby("gene_name")
        .agg(
            ms_peptide_count=("peptide", "nunique"),
            ms_cancer_peptide_count=("has_cancer", "sum"),
            ms_ebv_lcl_peptide_count=("has_ebv_lcl", "sum"),
            ms_healthy_somatic_peptide_count=("has_healthy_somatic", "sum"),
            ms_healthy_reproductive_peptide_count=("has_healthy_reproductive", "sum"),
            ms_healthy_thymus_peptide_count=("has_healthy_thymus", "sum"),
        )
        .reset_index()
    )

    # Merge tissue detail and PMIDs
    if somatic_tissues_per_gene is not None:
        gene_agg = gene_agg.merge(somatic_tissues_per_gene, on="gene_name", how="left")
        gene_agg["ms_healthy_somatic_tissue_count"] = (
            gene_agg["ms_healthy_somatic_tissue_count"].fillna(0).astype(int)
        )
        gene_agg["ms_healthy_somatic_tissues"] = gene_agg["ms_healthy_somatic_tissues"].fillna("")
    else:
        gene_agg["ms_healthy_somatic_tissue_count"] = 0
        gene_agg["ms_healthy_somatic_tissues"] = ""
    if pmids_per_gene is not None:
        gene_agg = gene_agg.merge(pmids_per_gene, on="gene_name", how="left")
        gene_agg["ms_pmids"] = gene_agg["ms_pmids"].fillna("")
    else:
        gene_agg["ms_pmids"] = ""

    # Classify
    gene_agg["ms_restriction"] = gene_agg.apply(_classify_gene_ms_restriction, axis=1)

    return gene_agg


def _classify_gene_ms_restriction(row: pd.Series) -> str:
    """Classify a single gene's MS restriction level."""
    n_pep = int(row.get("ms_peptide_count", 0))
    if n_pep == 0:
        return "NO_MS_DATA"

    n_somatic = int(row.get("ms_healthy_somatic_peptide_count", 0))
    n_repro = int(row.get("ms_healthy_reproductive_peptide_count", 0))
    n_thymus = int(row.get("ms_healthy_thymus_peptide_count", 0))
    n_cancer = int(row.get("ms_cancer_peptide_count", 0))
    n_somatic_tissues = int(row.get("ms_healthy_somatic_tissue_count", 0))

    if n_somatic == 0:
        if n_repro > 0 or n_thymus > 0:
            return "EXPECTED_TISSUE"
        if n_cancer > 0:
            return "CANCER_ONLY"
        return "NO_MS_DATA"

    if n_somatic == 1 and n_somatic_tissues <= 1:
        return "SINGLETON_HEALTHY"
    return "RECURRENT_HEALTHY"
