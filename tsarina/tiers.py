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

"""Mass-spec restriction assessment over oncoref's HPA restriction axes.

oncoref owns the protein-IHC and RNA restriction calls. Tsarina classifies
public immunopeptidomics evidence and can synthesize an explicitly MS-aware
confidence value for target selection.

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
  RECURRENT_HEALTHY / UNCLASSIFIED_MS / NO_MS_DATA.

**Synthesized**:

- ``restriction``: best available tissue category from protein > RNA.
- ``restriction_confidence``: HIGH / MODERATE / LOW based on modality
  agreement and quality.
"""

from __future__ import annotations

import pandas as pd
from oncoref.cta_tissues import SAFETY_TISSUE_GROUPS as SAFETY_TISSUE_GROUPS

from .tissues import HPA_EXPRESSION_FLOOR_NTPM

# ── Constants ──────────────────────────────────────────────────────────────

#: Tissue restriction categories (shared across protein, RNA, synthesis).
RESTRICTION_VALUES: list[str] = ["TESTIS", "PLACENTAL", "REPRODUCTIVE", "SOMATIC"]

#: The reproductive-tissue restriction categories. RNA's coarse "REPRODUCTIVE"
#: call only *agrees* with a protein call that is itself reproductive.
_REPRODUCTIVE_RESTRICTIONS: frozenset[str] = frozenset({"TESTIS", "PLACENTAL", "REPRODUCTIVE"})

#: RNA restriction quality levels, ordered strictest → loosest.
RNA_RESTRICTION_LEVELS: list[str] = ["STRICT", "MODERATE", "PERMISSIVE", "LEAKY"]

#: MS restriction classifications.
MS_RESTRICTION_VALUES: list[str] = [
    "CANCER_ONLY",
    "EXPECTED_TISSUE",
    "SINGLETON_HEALTHY",
    "RECURRENT_HEALTHY",
    # Has MS evidence, but every observation is from an unclassified source
    # (e.g. cell-line-only). Distinct from NO_MS_DATA (no MS evidence at all).
    "UNCLASSIFIED_MS",
    "NO_MS_DATA",
]

#: Restriction confidence levels.
CONFIDENCE_VALUES: list[str] = ["HIGH", "MODERATE", "LOW"]

#: Vital-organ tissue names as they appear in MS (immunopeptidome) source-tissue
#: fields — the MS-side counterpart of :data:`SAFETY_TISSUE_GROUPS` (which uses
#: HPA tissue names).  A healthy-tissue MS hit in one of these is a vital-organ
#: off-target signal.  Public so downstream consumers (e.g. vaxrank safety
#: scoring) can screen against the same vocabulary.  See tsarina#76 / vaxrank#303.
VITAL_TISSUE_MS_NAMES: frozenset[str] = frozenset(
    {
        "brain",
        "central nervous system (cns)",
        "cerebellum",
        "heart",
        "lung",
        "liver",
        "pancreas",
    }
)


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
        # Credit RNA's coarse "REPRODUCTIVE" as agreeing with a finer protein
        # call only when that call is itself reproductive (TESTIS/PLACENTAL).
        # A SOMATIC protein call genuinely *disagrees* with reproductive RNA, so
        # it must not earn the agreement bonus (which would otherwise inflate a
        # least-safe SOMATIC call to HIGH confidence and pass the selection
        # filter). This is a defensive gate against a future upstream
        # evidence refresh producing one.
        rna_agrees = rna_r == tissue or (
            rna_r == "REPRODUCTIVE" and tissue in _REPRODUCTIVE_RESTRICTIONS
        )
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

    # tsarina#114: cap HIGH when the only evidence is RNA below the expression
    # floor. The scorer credits any STRICT RNA equally regardless of level, so a
    # gene expressed at ~1-2 nTPM (never_expressed, no protein, no MS) otherwise
    # earns HIGH from near-noise RNA -- over-stating a restriction call built on
    # weak evidence. Genes with protein or MS evidence are untouched: an MS
    # peptide hit means real expression despite low RNA.
    has_ms_evidence = bool(ms_r) and ms_r not in ("NO_MS_DATA", "NO_DATA", "")
    if confidence == "HIGH" and not protein_has_data and not has_ms_evidence:
        rna_max = pd.to_numeric(row.get("rna_max_ntpm"), errors="coerce")
        if pd.notna(rna_max) and rna_max < HPA_EXPRESSION_FLOOR_NTPM:
            confidence = "MODERATE"

    return (tissue, confidence)


# ── Main assignment function ──────────────────────────────────────────────


def assign_all_axes(df: pd.DataFrame) -> pd.DataFrame:
    """Re-synthesize restriction/confidence from oncoref HPA axes plus MS.

    Kept as a compatibility helper for callers that explicitly want an MS-aware
    synthesis. It no longer derives protein or RNA restriction columns; those
    definitions belong exclusively to oncoref.

    Parameters
    ----------
    df
        An oncoref CTA evidence frame carrying ``protein_restriction``,
        ``protein_reliability``, ``rna_restriction``, and
        ``rna_restriction_level``. ``ms_restriction`` is optional.

    Returns
    -------
    pd.DataFrame
        Copy with an MS-aware ``restriction`` and ``restriction_confidence``.
    """
    out = df.copy()
    required = {
        "protein_restriction",
        "protein_reliability",
        "rna_restriction",
        "rna_restriction_level",
    }
    missing = required - set(out.columns)
    if missing:
        raise ValueError(
            f"assign_all_axes requires oncoref HPA restriction column(s): {sorted(missing)}"
        )

    if "ms_restriction" not in out.columns:
        out["ms_restriction"] = "NO_MS_DATA"

    synth = out.apply(synthesize_restriction, axis=1, result_type="expand")
    out["restriction"] = synth[0]
    out["restriction_confidence"] = synth[1]
    return out


# ── MS safety aggregation ─────────────────────────────────────────────────


def aggregate_gene_ms_safety(
    classified_hits: pd.DataFrame,
    peptide_gene_map: pd.DataFrame,
    exclusive_peptide_gene_map: pd.DataFrame | None = None,
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
    exclusive_peptide_gene_map
        Optional CTA-exclusive peptide-to-gene mapping. When supplied, adds
        strict CTA-exclusive MS count columns alongside the all-CTA counts.

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
        "ms_cta_exclusive_peptide_count",
        "ms_cta_exclusive_cancer_peptide_count",
        "ms_cta_exclusive_healthy_somatic_peptide_count",
        "ms_cta_exclusive_healthy_reproductive_peptide_count",
        "ms_cta_exclusive_healthy_thymus_peptide_count",
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
    gene_agg = _attach_cta_exclusive_ms_counts(
        gene_agg,
        classified_hits=classified_hits,
        exclusive_peptide_gene_map=exclusive_peptide_gene_map,
    )

    return gene_agg


def _attach_cta_exclusive_ms_counts(
    gene_agg: pd.DataFrame,
    *,
    classified_hits: pd.DataFrame,
    exclusive_peptide_gene_map: pd.DataFrame | None,
) -> pd.DataFrame:
    count_columns = [
        "ms_cta_exclusive_peptide_count",
        "ms_cta_exclusive_cancer_peptide_count",
        "ms_cta_exclusive_healthy_somatic_peptide_count",
        "ms_cta_exclusive_healthy_reproductive_peptide_count",
        "ms_cta_exclusive_healthy_thymus_peptide_count",
    ]
    if gene_agg.empty:
        for col in count_columns:
            gene_agg[col] = pd.Series(dtype="int64")
        return gene_agg

    if exclusive_peptide_gene_map is None or exclusive_peptide_gene_map.empty:
        for col in count_columns:
            gene_agg[col] = 0
        return gene_agg

    required_hit_cols = {
        "peptide",
        "src_cancer",
        "src_healthy_tissue",
        "src_healthy_reproductive",
        "src_healthy_thymus",
    }
    required_map_cols = {"peptide", "gene_name"}
    if not required_hit_cols <= set(classified_hits.columns) or not required_map_cols <= set(
        exclusive_peptide_gene_map.columns
    ):
        for col in count_columns:
            gene_agg[col] = 0
        return gene_agg

    exclusive_merged = classified_hits.merge(
        exclusive_peptide_gene_map[["peptide", "gene_name"]],
        on="peptide",
    )
    if exclusive_merged.empty:
        for col in count_columns:
            gene_agg[col] = 0
        return gene_agg

    exclusive_pep_gene = (
        exclusive_merged.groupby(["gene_name", "peptide"])
        .agg(
            has_cancer=("src_cancer", "any"),
            has_healthy_somatic=("src_healthy_tissue", "any"),
            has_healthy_reproductive=("src_healthy_reproductive", "any"),
            has_healthy_thymus=("src_healthy_thymus", "any"),
        )
        .reset_index()
    )
    exclusive_agg = (
        exclusive_pep_gene.groupby("gene_name")
        .agg(
            ms_cta_exclusive_peptide_count=("peptide", "nunique"),
            ms_cta_exclusive_cancer_peptide_count=("has_cancer", "sum"),
            ms_cta_exclusive_healthy_somatic_peptide_count=("has_healthy_somatic", "sum"),
            ms_cta_exclusive_healthy_reproductive_peptide_count=(
                "has_healthy_reproductive",
                "sum",
            ),
            ms_cta_exclusive_healthy_thymus_peptide_count=("has_healthy_thymus", "sum"),
        )
        .reset_index()
    )
    gene_agg = gene_agg.merge(exclusive_agg, on="gene_name", how="left")
    for col in count_columns:
        gene_agg[col] = gene_agg[col].fillna(0).astype(int)
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
        # MS peptides exist (n_pep > 0) but none carry a classifiable source
        # flag (cancer/reproductive/thymus/somatic) -- e.g. cell-line-only
        # evidence. Not the same as having no MS data.
        return "UNCLASSIFIED_MS"

    if n_somatic == 1 and n_somatic_tissues <= 1:
        return "SINGLETON_HEALTHY"
    return "RECURRENT_HEALTHY"
