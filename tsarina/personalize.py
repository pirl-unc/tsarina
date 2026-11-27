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

"""Patient-level personalization: from tumor characteristics to prioritized targets.

The clinical entry point.  Given a patient's HLA type, tumor CTA expression,
detected mutations, and viral status, returns a ranked list of (peptide, HLA)
targets annotated with public MS evidence, predicted presentation, and an
explicit confidence tier.

Guardrails applied by default:

1. **CTA exclusivity** — only k-mers that do not appear in any non-CTA protein
   reach the candidate pool (via :func:`tsarina.peptides.cta_exclusive_peptides`).
2. **Viral human-exclusivity** — viral k-mers that match the human proteome
   are dropped (via :func:`tsarina.viral.human_exclusive_viral_peptides`).
3. **Current CTA admission** — default membership uses current OncoRef decisions.
   Historical HPA confidence is an optional filter. An explicitly requested
   historical clinical target outside the current set retains a flagged
   compatibility path with current decision reasons and specificity caveats.
4. **Tumor-specificity filter** — peptides observed by mass spec on healthy
   non-reproductive tissue are excluded entirely (not merely penalized).
5. **Mandatory MHC scoring** — when ``score_presentation=True`` the scorer
   runs; failures raise loudly rather than silently dropping the gate.
6. **Explicit tiers** — candidates are assigned to tier 1-4 based on
   presentation percentile + evidence category; tier 4 (weak) is dropped
   from the default output.
7. **Deterministic sort** — output order is stable across runs
   (tier, ms_hit_count, percentile, peptide, best_allele).

Optional extra filters:

- ``mtec_matrix_path`` — when given, CTAs are restricted to those with mean
  mTEC expression below ``mtec_max_tpm`` (thymic-tolerance blacklist).

Typical usage::

    from tsarina import personalized_targets

    targets = personalized_targets(
        hla_alleles=["HLA-A*02:01", "HLA-B*07:02"],
        cta_expression={"MAGEA4": 142.5, "PRAME": 87.3},
        mutations=["KRAS G12D"],
        viruses=["hpv16"],
    )
"""

from __future__ import annotations

import sys
import warnings
from collections.abc import Iterable
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from .scoring import PRESENTATION_PERCENTILE_THRESHOLDS

_TIER_LABELS: dict[int, str] = {
    1: "STRONG",
    2: "MODERATE",
    3: "CANDIDATE",
    4: "WEAK",
}

_OUTPUT_COLUMNS: tuple[str, ...] = (
    "peptide",
    "length",
    "category",
    "source",
    "source_detail",
    "source_tpm",
    "flag_reason",
    "unsafe_near_self_1mm_pct",
    "clinical_cleared_epitope",
    "cta_safety_note",
    "ms_evidence_status",
    "ms_attribution_tier",
    "ms_hit_count",
    "ms_exact_mhc_hits",
    "ms_donor_set_hits",
    "ms_alleles",
    "ms_allele_count",
    "ms_in_cancer",
    "ms_in_healthy_tissue",
    "ms_in_healthy_reproductive_female",
    "ms_in_healthy_reproductive_male",
    "ms_detectability_note",
    "best_allele",
    "presentation_percentile",
    "presentation_score",
    "affinity_nm",
    "tier",
    "tier_label",
    "tier_reason",
)

_NEAR_SELF_OUTPUT_COLUMNS: tuple[str, ...] = (
    "near_self_1mm_count",
    "near_self_1mm_genes",
    "near_self_1mm_peptides",
)


@lru_cache(maxsize=1)
def _proteoform_group_labels() -> dict[str, str]:
    """Member gene symbol -> the group's preferred display symbol.

    Both halves come from oncoref: current selected entities for group
    membership (the same decisions :mod:`tsarina.spanning` uses for the
    panel workflow) and ``proteoform_symbol`` for the name that survives
    the collapse -- a curated alias where one exists (``CTAG1A/CTAG1B``
    -> ``NY-ESO-1``), else the prefix-contracted members
    (``XAGE1A/XAGE1B`` -> ``XAGE1A/B``). Neither is restated here, so
    tsarina can't drift from the canonical naming.
    """
    from oncoref.proteoforms import proteoform_symbol

    from .cta_selection import selection_proteoform_groups

    return {
        member: proteoform_symbol(label)
        for label, members in selection_proteoform_groups().items()
        for member in members
    }


def _apply_proteoform_rollup(frame: pd.DataFrame) -> pd.DataFrame:
    """Collapse identical-protein CTA paralogs into one group-labeled row.

    NY-ESO-1 is CTAG1A and CTAG1B; XAGE1A/B is XAGE1A and XAGE1B; SSX2/B
    is SSX2 and SSX2B. Each pair translates to a byte-identical protein,
    so naming either member yields the same peptides -- reporting them as
    separate sources either double-counts the same finding (when the
    caller named both) or hides that the peptide isn't unique to the one
    they named (when they named one). Both are relabeled to the group's
    preferred symbol, and a (group, peptide, length) duplicate collapses
    to a single row.

    Only rows whose source is a known group member are touched; viral,
    mutant, and ungrouped CTA rows pass through untouched (a same-peptide
    duplicate there can be a genuinely distinct source and isn't ours to
    merge).
    """
    labels = _proteoform_group_labels()
    if "source" not in frame.columns:
        return frame
    grouped_mask = frame["source"].isin(labels)
    if not grouped_mask.any():
        return frame

    rest = frame[~grouped_mask]
    members = frame[grouped_mask].copy()
    members["source"] = members["source"].map(labels)

    agg: dict[str, object] = {}
    if "source_detail" in members.columns:
        agg["source_detail"] = lambda s: ";".join(sorted({str(v) for v in s if pd.notna(v)}))
    if "source_tpm" in members.columns:
        # Identical proteins, so an RNA quantifier's reads split between
        # the loci more or less arbitrarily: take the highest any member
        # reported rather than summing (which would inflate a single real
        # signal that happened to be split two ways). All-NaN -- every
        # member came in without a TPM -- stays NaN.
        agg["source_tpm"] = "max"
    for col in ("category", "flag_reason"):
        if col in members.columns:
            agg[col] = "first"

    rolled = members.groupby(["source", "peptide", "length"], as_index=False, sort=False).agg(agg)
    return pd.concat([rest, rolled], ignore_index=True)


_MS_ATTR_SHORT: dict[str, str] = {
    "exact_requested_mhc": "exact",
    "donor_set_ambiguous": "donor_set",
    "predicted_only": "pred_only",
    "other_mhc_only": "other_mhc",
    "no_match": "no_match",
    "evidence_not_loaded": "unloaded",
}


def _format_ms_safety(row: pd.Series) -> str:
    hits = row.get("ms_hit_count")
    if pd.isna(hits) or int(hits) <= 0:
        return "-"
    in_cancer = bool(row.get("ms_in_cancer")) if pd.notna(row.get("ms_in_cancer")) else False
    in_healthy = (
        bool(row.get("ms_in_healthy_tissue"))
        if pd.notna(row.get("ms_in_healthy_tissue"))
        else False
    )
    in_repro_f = (
        bool(row.get("ms_in_healthy_reproductive_female"))
        if pd.notna(row.get("ms_in_healthy_reproductive_female"))
        else False
    )
    in_repro_m = (
        bool(row.get("ms_in_healthy_reproductive_male"))
        if pd.notna(row.get("ms_in_healthy_reproductive_male"))
        else False
    )
    if in_healthy and in_cancer:
        return "cancer+healthy!"
    if in_healthy:
        return "healthy!"
    repro_tags: list[str] = []
    if in_repro_f:
        repro_tags.append("repro_f")
    if in_repro_m:
        repro_tags.append("repro_m")
    if repro_tags:
        suffix = "+".join(repro_tags)
        return f"cancer+{suffix}" if in_cancer else suffix
    if in_cancer:
        return "cancer"
    return "unclassified"


def format_table(df: pd.DataFrame) -> str:
    """Render a :func:`personalized_targets` result as a compact,
    fixed-width table -- one row per peptide, sorted as given (already
    tier-ranked by the caller), no box-drawing.

    Mirrors the plain, minimal-decoration style of
    :func:`hitlist.pmhc_query.format_table` (dashes under headers, no
    grid lines) since the two tools' output is often read side by side,
    but tsarina's own schema -- one row per peptide with a single
    already-chosen best allele, rather than hitlist's per-(gene, allele)
    rows -- doesn't map onto that function directly, so this is its own
    implementation rather than a shared one.

    A ``cta_flagged`` row's source gets a trailing ``*``, with a one-line
    footnote explaining it -- the full ``flag_reason`` sentence is long
    enough that inlining it into every such row would defeat the point of
    a compact table; use ``--format csv`` for that.
    """
    if df.empty:
        return "(no targets)"

    columns: list[tuple[str, str | None]] = [
        ("peptide", "peptide"),
        ("category", "category"),
        ("source", "source"),
        ("tier", "tier_label"),
        ("allele", "best_allele"),
        ("affinity_nM", "affinity_nm"),
        ("pct_rank", "presentation_percentile"),
        ("ms_hits", "ms_hit_count"),
    ]
    if "ms_attribution_tier" in df.columns:
        columns.append(("ms_attr", "ms_attribution_tier"))
    if {"ms_in_cancer", "ms_in_healthy_tissue"} & set(df.columns):
        columns.append(("ms_safety", None))
    if "near_self_1mm_count" in df.columns:
        columns.append(("1mm_self", "near_self_1mm_count"))
    if "tier_reason" in df.columns:
        columns.append(("reason", "tier_reason"))

    def _fmt(header: str, value: object, row: pd.Series) -> str:
        if header == "ms_safety":
            return _format_ms_safety(row)
        if pd.isna(value):
            return ""
        if header in ("ms_hits", "1mm_self"):
            return f"{int(value)}"
        if header in ("affinity_nM", "pct_rank"):
            return f"{float(value):.2f}"
        if header == "ms_attr":
            token = str(value).strip()
            return _MS_ATTR_SHORT.get(token, token)
        return str(value)

    has_flagged = (df["category"] == "cta_flagged").any() if "category" in df.columns else False
    rows: list[list[str]] = []
    for _, row in df.iterrows():
        cells = [
            _fmt(key, row.get(src_col) if src_col is not None else None, row)
            for key, src_col in columns
        ]
        if has_flagged and row.get("category") == "cta_flagged":
            cells[2] += "*"  # "source" is the third column
        rows.append(cells)

    headers = [key for key, _ in columns]
    widths = [max(len(h), *(len(r[i]) for r in rows)) for i, h in enumerate(headers)]
    lines = [
        "  ".join(h.ljust(w) for h, w in zip(headers, widths)),
        "  ".join("-" * w for w in widths),
    ]
    lines.extend("  ".join(cell.ljust(w) for cell, w in zip(row, widths)) for row in rows)
    if has_flagged:
        lines.append("")
        lines.append(
            "* flagged clinical-target CTA, outside the current CTA selection -- "
            "see flag_reason in --format csv for why."
        )
    if "ms_detectability_note" in df.columns:
        cys_notes = df["ms_detectability_note"].fillna("").astype(str).str.strip()
        n_cys = int(cys_notes.ne("").sum())
        if n_cys > 0:
            lines.append("")
            lines.append(
                f"Note: {n_cys} zero-MS peptide(s) contain Cys (unmodified MS under-detects "
                "Cys oxidation/cysteinylation; see ms_detectability_note in --format csv)."
            )
    if "cta_safety_note" in df.columns:
        safety_notes = df["cta_safety_note"].fillna("").astype(str).str.strip()
        n_safety = int(safety_notes.ne("").sum())
        if n_safety > 0:
            lines.append("")
            lines.append(
                f"Note: {n_safety} row(s) carry OncoRef CTA safety caveats (nuclear snRNA-seq, "
                "epitope-scoped clinical clearance, or >=40% 1-mm near-self overlap; "
                "see cta_safety_note in --format csv)."
            )
    if "near_self_1mm_count" in df.columns:
        mm1_counts = pd.to_numeric(df["near_self_1mm_count"], errors="coerce").fillna(0)
        n_mm1 = int(mm1_counts.gt(0).sum())
        if n_mm1 > 0:
            lines.append("")
            lines.append(
                f"Note: {n_mm1} peptide(s) have 1-aa near-self matches in the Unsafe proteome "
                "(see near_self_1mm_genes / near_self_1mm_peptides in --format csv)."
            )
    return "\n".join(lines)


def _cta_flagged_gene_peptides(
    gene_names: Iterable[str], *, ensembl_release: int, lengths: tuple[int, ...]
) -> pd.DataFrame:
    """Generate peptides for a clinical-target CTA that oncoref excludes
    from the current CTA gene universe (historical clinical-target
    exclusion, for example).

    :func:`tsarina.peptides.cta_peptides` (and the ``_cta_gene_ids_for_names``
    resolver it uses) gates gene-name resolution on the strict CTA gene-ID
    set, so a flagged gene resolves to no genes at all through that path --
    it needs its own resolution here, not a parameter override there.

    Returns the same ``gene_name``/``gene_id``/``peptide``/``length`` shape
    as ``cta_peptides``, minus the flank columns this caller doesn't need.
    No exclusivity screening against non-CTA proteins is applied (see the
    ``flag_reason`` caveat this gets merged into by the caller).
    """
    from pyensembl import EnsemblRelease

    from .gene_sets import is_coding_transcript
    from .loader import cta_dataframe
    from .peptides import AA20

    wanted = {str(g).strip() for g in gene_names if str(g).strip()}
    columns = ["gene_name", "gene_id", "peptide", "length"]
    if not wanted:
        return pd.DataFrame(columns=columns)

    df = cta_dataframe()
    gene_ids: dict[str, str] = {}
    if "Symbol" in df.columns and "Ensembl_Gene_ID" in df.columns:
        for symbol_cell, id_cell in df[["Symbol", "Ensembl_Gene_ID"]].itertuples(
            index=False, name=None
        ):
            hit = wanted & {s.strip() for s in str(symbol_cell).split(";") if s.strip()}
            if not hit:
                continue
            ids = [i.strip() for i in str(id_cell).split(";") if i.strip()]
            if ids:
                for symbol in hit:
                    gene_ids.setdefault(symbol, ids[0])

    ensembl = EnsemblRelease(ensembl_release)
    rows: list[dict] = []
    for gene_name, gene_id in gene_ids.items():
        try:
            gene = ensembl.gene_by_id(gene_id)
        except ValueError:
            continue
        best_transcript = None
        best_length = 0
        for t in gene.transcripts:
            if not is_coding_transcript(t):
                continue
            try:
                seq = t.protein_sequence
            except (ValueError, KeyError, TypeError):
                continue
            if seq and len(seq) > best_length:
                best_transcript = t
                best_length = len(seq)
        if best_transcript is None:
            continue
        protein = best_transcript.protein_sequence
        if not protein:
            continue
        for k in lengths:
            for i in range(len(protein) - k + 1):
                pep = protein[i : i + k]
                if set(pep).issubset(AA20):
                    rows.append(
                        {"gene_name": gene_name, "gene_id": gene_id, "peptide": pep, "length": k}
                    )
    return pd.DataFrame(rows, columns=columns)


def _cta_flag_rationale(genes: dict[str, float]) -> dict[str, str]:
    """Current decision reasons, with historical rationale only as a fallback."""
    if not genes:
        return {}
    from oncoref import cta as _oncoref_cta

    df = _oncoref_cta.cta_df()
    rows = df[df["Symbol"].isin(genes) & df["specificity_rationale"].notna()]
    rationale = dict(zip(rows["Symbol"], rows["specificity_rationale"]))
    from .cta_selection import selection_gene_table

    for row in selection_gene_table().itertuples(index=False):
        if row.Symbol in genes:
            classification = getattr(row, "selection_classification", None)
            if pd.isna(classification) or not classification:
                category = (
                    "unassigned" if pd.isna(row.selection_category) else row.selection_category
                )
            else:
                category = str(classification)
            reasons = "; ".join(row.selection_reasons)
            rationale[row.Symbol] = f"current OncoRef category: {category}" + (
                f" ({reasons})" if reasons else ""
            )
    return {
        gene: rationale.get(gene)
        or "not in OncoRef's current CTA selection; no curated rationale available"
        for gene in genes
    }


def personalized_targets(
    hla_alleles: list[str],
    cta_expression: dict[str, float | None] | None = None,
    mutations: list[str] | None = None,
    viruses: list[str] | None = None,
    lengths: tuple[int, ...] = (8, 9, 10, 11),
    ensembl_release: int | None = None,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    evidence_package: str | Path | None = None,
    mhc_class: str = "I",
    min_cta_tpm: float = 2.0,
    min_restriction_confidence: Iterable[str] | None = None,
    mtec_matrix_path: str | Path | None = None,
    mtec_max_tpm: float = 1.0,
    require_human_exclusive_viral: bool = True,
    enforce_tumor_specificity: bool = True,
    score_presentation: bool = True,
    skip_ms_evidence: bool = False,
    predictor: str = "mhcflurry",
    drop_weak_tier: bool = True,
    show_progress: bool = True,
    proteoform_rollup: bool = True,
    *,
    host_species: str = "human",
    mhc_species: str | None = None,
    patient_sex: str | None = None,
    sex_eligibility: str = "all_adult",
    host_context: object | None = None,
    pregnant: bool = False,
    lactating: bool = False,
    protect_tissues: Iterable[str] = (),
    accept_loss_tissues: Iterable[str] = (),
    absent_tissues: Iterable[str] = (),
    allow_incomplete_scope: bool = False,
    fusions: Iterable[dict] | pd.DataFrame | None = None,
    require_host_exclusive_fusions: bool = True,
    ms_evidence_scope: str = "all",
    ms_observations: pd.DataFrame | None = None,
    custom_targets: pd.DataFrame | None = None,
    host_proteins: dict[str, str] | list[str] | None = None,
    capability_table: pd.DataFrame | None = None,
    screen_near_self_1mm: bool = False,
    max_near_self_1mm_matches: int | None = None,
    near_self_background: str = "Unsafe",
    near_self_background_proteins: pd.DataFrame | dict[str, str] | Iterable[str] | None = None,
) -> pd.DataFrame:
    """Build a personalized, tier-ranked target list for a single patient.

    Parameters
    ----------
    hla_alleles
        Patient's HLA class I alleles (e.g. ``["HLA-A*02:01", "HLA-B*07:02"]``).
    cta_expression
        Dict mapping CTA gene symbol to tumor RNA expression in TPM.  Only
        genes passing all CTA gates (membership in the filter-passing CTA set,
        ``tpm >= min_cta_tpm``, ``restriction_confidence`` allowed, optional
        mTEC cutoff) contribute peptides.  Peptides are restricted to those
        exclusive to CTA proteins (not present in any non-CTA protein).
        Explicitly requested clinical targets outside the strict set use the
        separate ``cta_flagged`` path with visible specificity/overlap caveats;
        they bypass the HPA confidence gate but still obey tumor TPM and mTEC
        gates. A strict CTA rejected by a gate cannot switch to that path.
        Dropped source genes warn with the reason for exclusion.
        A ``None`` or ``NaN`` TPM means "no measured expression given" --
        the gene is included regardless of ``min_cta_tpm`` rather than
        excluded as if it measured zero.
    mutations
        List of hotspot labels (e.g. ``["KRAS G12D"]``).  Wildtype-identical
        k-mers and matches to any other reference-human coding protein are
        filtered in :mod:`tsarina.mutations` before scoring. Sequence
        exclusivity does not establish tumor presentation by itself.
    viruses
        List of virus keys (e.g. ``["hpv16", "ebv"]``).  By default only
        viral k-mers that do not appear in any human protein are kept; set
        ``require_human_exclusive_viral=False`` to relax this.
    lengths
        Peptide lengths (default 8-11).
    ensembl_release
        Ensembl release (default 112).
    iedb_path, cedar_path
        Explicit dataset paths.  Default auto-resolves via the hitlist
        registry.  When either is supplied, the direct scan path is used
        instead of the cached observations index.
    mhc_class
        MHC class filter for IEDB/CEDAR scanning (default ``"I"``).
    min_cta_tpm
        Minimum measured tumor CTA expression in TPM (default 2.0).
        This user-supplied tumor-expression gate is separate from admission.
    min_restriction_confidence
        Allowed ``restriction_confidence`` bins for CTAs.  Default
        ``None`` leaves current admission unchanged; pass explicit bins to filter.
    mtec_matrix_path
        Path to an mTEC gene TPM matrix (TSV).  When provided, CTAs are
        restricted to genes with mean mTEC TPM <= ``mtec_max_tpm``.
    mtec_max_tpm
        Maximum mean mTEC TPM allowed when ``mtec_matrix_path`` is set
        (default 1.0).
    require_human_exclusive_viral
        If True (default), use
        :func:`tsarina.viral.human_exclusive_viral_peptides`.
    enforce_tumor_specificity
        If True (default), drop peptides with any public MS evidence on
        healthy non-reproductive tissue.
    score_presentation
        If True (default), score peptide-allele pairs via topiary + mhctools.
        An ``ImportError`` from the scoring backend is propagated (no silent
        pass).
    skip_ms_evidence
        If True, do not look up IEDB/CEDAR evidence (dry-run path).
    predictor
        mhctools predictor key.  One of ``"mhcflurry"`` (default),
        ``"netmhcpan"``, ``"netmhcpan_el"``.  A warning is emitted when
        a non-mhcflurry predictor is chosen, since tier cutoffs are
        calibrated to mhcflurry's presentation percentile distribution.
    drop_weak_tier
        If True (default), tier-4 rows are dropped from the output.  Set
        False to retain them for diagnostics.
    show_progress
        If True (default), write a line to stderr before each stage
        (peptide generation per source, MS evidence lookup, presentation
        scoring) so a long-running call isn't silent. Set False for a
        quiet/library-embedded run.
    proteoform_rollup
        If True (default), CTAs that translate to a byte-identical
        protein are reported as one group under oncoref's preferred
        symbol for it (``CTAG1A`` + ``CTAG1B`` -> ``NY-ESO-1``,
        ``XAGE1A`` + ``XAGE1B`` -> ``XAGE1A/B``, ``SSX2`` + ``SSX2B``
        -> ``SSX2/B``).  Set False to keep one row per gene symbol.

    Returns
    -------
    pd.DataFrame
        Tier-ranked target list.  Columns:

        - ``peptide``, ``length``, ``category``, ``source``, ``flag_reason``,
          ``source_detail``, ``source_tpm``
        - ``ms_hit_count``, ``ms_alleles``, ``ms_allele_count``,
          ``ms_in_cancer``, ``ms_in_healthy_tissue``
        - ``best_allele``, ``presentation_percentile``,
          ``presentation_score``, ``affinity_nm``
        - ``tier`` (int 1-4, lower is better), ``tier_label`` (str),
          ``tier_reason`` (which gates were satisfied)

    Raises
    ------
    ImportError
        If ``score_presentation=True`` but the scoring backend (topiary +
        mhctools + the chosen predictor) is not installed.
    """
    from .ms_attribution import VALID_EVIDENCE_SCOPES
    from .species import (
        allele_prediction_capability,
        resolve_species_context,
        validate_alleles_for_species,
    )

    if ms_evidence_scope not in VALID_EVIDENCE_SCOPES:
        raise ValueError(
            f"ms_evidence_scope must be one of {sorted(VALID_EVIDENCE_SCOPES)}, "
            f"got {ms_evidence_scope!r}"
        )

    run_near_self_screen = bool(screen_near_self_1mm or max_near_self_1mm_matches is not None)
    if run_near_self_screen:
        from .near_self import VALID_NEAR_SELF_BACKGROUNDS

        if near_self_background not in VALID_NEAR_SELF_BACKGROUNDS:
            raise ValueError(
                f"near_self_background must be one of {VALID_NEAR_SELF_BACKGROUNDS}, "
                f"got {near_self_background!r}"
            )
        if max_near_self_1mm_matches is not None and max_near_self_1mm_matches < 0:
            raise ValueError(
                f"max_near_self_1mm_matches must be >= 0, got {max_near_self_1mm_matches!r}"
            )
    output_columns = list(_OUTPUT_COLUMNS) + (
        list(_NEAR_SELF_OUTPUT_COLUMNS) if run_near_self_screen else []
    )

    norm_patient_sex: str | None = None
    norm_sex_eligibility = "all_adult"
    if patient_sex is not None:
        from .cta_selection import _normalize_patient_sex, _normalize_sex_eligibility

        norm_patient_sex = _normalize_patient_sex(patient_sex)
        norm_sex_eligibility = _normalize_sex_eligibility(sex_eligibility)

    ctx_kwargs: dict[str, object] = {"mhc_species": mhc_species}
    if ensembl_release is not None:
        ctx_kwargs["ensembl_release"] = ensembl_release
    ctx = resolve_species_context(host_species, **ctx_kwargs)
    resolved_ensembl_release = int(ctx.ensembl_release) if ctx.ensembl_release is not None else 112

    has_host_flags = any(
        (
            bool(pregnant),
            bool(lactating),
            bool(protect_tissues),
            bool(accept_loss_tissues),
            bool(absent_tissues),
            bool(allow_incomplete_scope),
        )
    )
    if host_context is not None and has_host_flags:
        raise ValueError(
            "Pass either host_context or individual host flags "
            "(pregnant/lactating/protect_tissues/accept_loss_tissues/absent_tissues/"
            "allow_incomplete_scope), not both."
        )
    resolved_host = None
    if host_context is not None:
        from oncoref.host_context import HostContext as _HostContext
        from oncoref.species import resolve_species as _resolve_oncoref_species

        if not isinstance(host_context, _HostContext):
            raise ValueError("host_context must be an oncoref.host_context.HostContext instance")
        expected_species = _resolve_oncoref_species(ctx.host_species).name
        if host_context.species != expected_species:
            raise ValueError(
                f"host_context.species={host_context.species!r} does not match "
                f"host_species={expected_species!r}"
            )
        resolved_host = host_context
    elif has_host_flags:
        from .cta_selection import build_host_context

        resolved_host = build_host_context(
            patient_sex=norm_patient_sex,
            pregnant=pregnant,
            lactating=lactating,
            protect_tissues=protect_tissues,
            accept_loss_tissues=accept_loss_tissues,
            absent_tissues=absent_tissues,
            allow_incomplete_scope=allow_incomplete_scope,
            species=ctx.host_species,
        )
    if hla_alleles:
        validate_alleles_for_species(hla_alleles, ctx.mhc_species)
    allele_caps = allele_prediction_capability(
        hla_alleles,
        predictor=predictor,
        mhc_species=ctx.mhc_species,
        capability_table=capability_table,
    )

    nonhuman_cta_evaluated = False
    if not ctx.is_human_host:
        if cta_expression:
            from .cta_selection import CTA_selection_coverage

            cov = CTA_selection_coverage(species=ctx.host_species)
            if cov.get("selection_status") == "evaluated":
                nonhuman_cta_evaluated = True
            else:
                raise ValueError(
                    f"Human OncoRef CTA expression catalog cannot be silently transferred to "
                    f"host_species={ctx.host_species!r} "
                    f"(selection_status={cov.get('selection_status')!r}); "
                    "pass species-scoped custom_targets instead."
                )
        if mutations:
            raise ValueError(
                f"Human HOTSPOT_MUTATIONS cannot be silently transferred to "
                f"host_species={ctx.host_species!r}; pass species-scoped custom_targets instead."
            )

    frames: list[pd.DataFrame] = []

    if custom_targets is not None and not custom_targets.empty:
        req_cols = {"peptide", "category", "source"}
        missing = req_cols - set(custom_targets.columns)
        if missing:
            raise ValueError(f"custom_targets missing required columns: {sorted(missing)}")
        ct = custom_targets.copy()
        if "length" not in ct.columns:
            ct["length"] = ct["peptide"].astype(str).str.len()
        if "source_detail" not in ct.columns:
            ct["source_detail"] = ""
        if "source_tpm" not in ct.columns:
            ct["source_tpm"] = float("nan")
        frames.append(
            ct[["peptide", "length", "category", "source", "source_detail", "source_tpm"]]
        )

    def _report(message: str) -> None:
        if show_progress:
            print(message, file=sys.stderr)

    if score_presentation and predictor != "mhcflurry":
        warnings.warn(
            f"Predictor {predictor!r} is not mhcflurry; tier cutoffs "
            f"({PRESENTATION_PERCENTILE_THRESHOLDS}) are calibrated to mhcflurry's "
            "presentation-percentile distribution and may not translate 1:1 "
            "to other backends.",
            UserWarning,
            stacklevel=2,
        )

    # ── CTA peptides ────────────────────────────────────────────────────
    if cta_expression and nonhuman_cta_evaluated:
        from .cta_selection import CTA_selection, CTA_sex_restricted_antigens
        from .peptides import _reference_overlapping_peptides, gene_peptides
        from .viral import _proteome_kmers_from_sequences

        species_ctas = CTA_selection(category="CTA", species=ctx.host_species)
        all_species_cta_ids: set[str] = {
            str(gid)
            for members in species_ctas.get("gene_ids", ())
            if isinstance(members, (tuple, list))
            for gid in members
            if gid
        }
        host_passed_entity_ids: set[str] | None = None
        host_fail_by_key: dict[str, str] = {}
        if resolved_host is not None and resolved_host.protected:
            host_assessed = CTA_selection(
                category=("CTA", "Reproductive Antigen only"),
                species=ctx.host_species,
                host_context=resolved_host,
            )
            host_passed_entity_ids = set(
                host_assessed.loc[host_assessed["pass_host"].astype(bool), "entity_id"]
            )
            for row in host_assessed.loc[~host_assessed["pass_host"].astype(bool)].itertuples(
                index=False
            ):
                label = str(getattr(row, "label", "") or "").strip()
                status = str(getattr(row, "host_status", "excluded") or "excluded")
                reasons = str(getattr(row, "host_reasons", "") or "").strip()
                msg = f"HostContext {status}" + (f" ({reasons})" if reasons else "")
                if label:
                    host_fail_by_key[label.casefold()] = msg
        eligible_by_key: dict[str, tuple[str, tuple[str, ...]]] = {}
        sex_eligible_entity_ids: set[str] | None = None
        if norm_patient_sex is not None:
            sex_df = CTA_sex_restricted_antigens(
                sex=norm_patient_sex,
                eligibility=norm_sex_eligibility,
                include_reproductive_only=True,
                species=ctx.host_species,
            )
            sex_eligible_entity_ids = set(sex_df["entity_id"])
            for row in sex_df.itertuples(index=False):
                if (
                    host_passed_entity_ids is not None
                    and row.entity_id not in host_passed_entity_ids
                ):
                    continue
                label = str(getattr(row, "label", "") or "").strip()
                members = getattr(row, "gene_ids", ()) or ()
                if label and members:
                    eligible_by_key[label.casefold()] = (label, tuple(str(g) for g in members))
        for row in species_ctas.itertuples(index=False):
            if sex_eligible_entity_ids is not None and row.entity_id not in sex_eligible_entity_ids:
                continue
            if host_passed_entity_ids is not None and row.entity_id not in host_passed_entity_ids:
                continue
            label = str(getattr(row, "label", "") or "").strip()
            members = getattr(row, "gene_ids", ()) or ()
            if label and members:
                eligible_by_key[label.casefold()] = (label, tuple(str(g) for g in members))

        admitted_tpms: dict[str, float] = {}
        admitted_ids: set[str] = set()
        id_to_label: dict[str, str] = {}
        for gene, raw_tpm in cta_expression.items():
            tpm = float("nan") if raw_tpm is None else float(raw_tpm)
            gene_key = str(gene).strip().casefold()
            hit = eligible_by_key.get(gene_key)
            if hit is None:
                if gene_key in host_fail_by_key:
                    reason = host_fail_by_key[gene_key]
                else:
                    reason = (
                        f"not in OncoRef {ctx.host_species} CTA"
                        + (
                            f" / {norm_patient_sex} ({norm_sex_eligibility}) eligible"
                            if norm_patient_sex is not None
                            else ""
                        )
                        + " set"
                    )
                warnings.warn(f"--cta {gene} dropped: {reason}", UserWarning, stacklevel=2)
                continue
            if not pd.isna(tpm) and tpm < min_cta_tpm:
                warnings.warn(
                    f"--cta {gene} dropped: TPM {tpm:g} < min_cta_tpm {min_cta_tpm:g}",
                    UserWarning,
                    stacklevel=2,
                )
                continue
            canon_label, gids = hit
            admitted_tpms[canon_label] = tpm
            admitted_ids.update(gids)
            for gid in gids:
                id_to_label[gid] = canon_label

        if admitted_ids:
            _report(
                f"Generating {ctx.host_species} CTA peptides for {len(admitted_tpms)} gene(s): "
                f"{', '.join(sorted(admitted_tpms))}..."
            )
            sp_peps = gene_peptides(
                sorted(admitted_ids),
                ensembl_release=resolved_ensembl_release,
                lengths=lengths,
                species=ctx.pyensembl_species,
            )
            if not sp_peps.empty:
                if host_proteins is not None:
                    host_kmers = _proteome_kmers_from_sequences(host_proteins, tuple(lengths))
                    overlapping = set(sp_peps["peptide"]) & host_kmers
                else:
                    overlapping = _reference_overlapping_peptides(
                        set(sp_peps["peptide"]),
                        ensembl_release=resolved_ensembl_release,
                        lengths=lengths,
                        species=ctx.pyensembl_species,
                        exclude_gene_ids=all_species_cta_ids | admitted_ids,
                    )
                sp_peps = sp_peps[~sp_peps["peptide"].isin(overlapping)].copy()
            if not sp_peps.empty:
                sources = sp_peps["gene_id"].map(id_to_label).fillna(sp_peps["gene_name"])
                frames.append(
                    pd.DataFrame(
                        {
                            "peptide": sp_peps["peptide"],
                            "length": sp_peps["length"],
                            "category": "cta",
                            "source": sources,
                            "source_detail": sp_peps["gene_id"],
                            "source_tpm": sources.map(admitted_tpms),
                        }
                    )
                )

    elif cta_expression:
        from .gene_sets import (
            CTA_by_axes,
            CTA_clinical_target_gene_names,
            CTA_excluded_gene_names,
            CTA_gene_names,
            CTA_unfiltered_gene_names,
        )
        from .peptides import cta_exclusive_peptides

        # A caller who didn't supply a TPM for a gene (--cta GENE with no
        # =TPM on the CLI, or an explicit None here) means "include it
        # regardless of min_cta_tpm", not "measured zero expression, gate
        # it out." Normalize to NaN once so every filter below can use a
        # single float comparison.
        cta_expression = {
            gene: (float("nan") if tpm is None else tpm) for gene, tpm in cta_expression.items()
        }

        # Membership determines the generation path before optional filters.
        # A strict CTA rejected below must never return as an unscreened
        # clinical target merely because it also has a clinical reference.
        strict_ctas = CTA_gene_names()
        clinical_target_ctas = CTA_clinical_target_gene_names() - strict_ctas
        valid_ctas = strict_ctas.copy()
        sex_repro_ctas: set[str] = set()
        drop_reasons = {}

        if norm_patient_sex is not None:
            from .cta_selection import CTA_sex_eligible_gene_names, CTA_sex_restricted_antigens

            sex_eligible = CTA_sex_eligible_gene_names(
                norm_patient_sex,
                eligibility=norm_sex_eligibility,
                include_reproductive_only=True,
            )
            sex_table = CTA_sex_restricted_antigens(include_reproductive_only=True)
            sex_class_by_label = dict(
                zip(
                    sex_table["label"].astype(str),
                    sex_table["sex_restriction_class"].astype(str),
                )
            )
            sex_repro_ctas = sex_eligible - strict_ctas
            for gene in set(cta_expression) & (valid_ctas - sex_eligible):
                sex_cls = sex_class_by_label.get(gene, "not sex-eligible")
                drop_reasons[gene] = (
                    f"not {norm_patient_sex} {norm_sex_eligibility} eligible "
                    f"(OncoRef sex_restriction_class: {sex_cls})"
                )
            valid_ctas &= sex_eligible

        if resolved_host is not None and resolved_host.protected:
            from .cta_selection import CTA_selection, selection_gene_table

            host_assessed = CTA_selection(
                category=("CTA", "Reproductive Antigen only"),
                species="human",
                host_context=resolved_host,
            )
            passed_entities = set(
                host_assessed.loc[host_assessed["pass_host"].astype(bool), "entity_id"]
            )
            fail_by_entity: dict[str, str] = {}
            for row in host_assessed.loc[~host_assessed["pass_host"].astype(bool)].itertuples(
                index=False
            ):
                status = str(getattr(row, "host_status", "excluded") or "excluded")
                reasons = str(getattr(row, "host_reasons", "") or "").strip()
                fail_by_entity[str(row.entity_id)] = f"HostContext {status}" + (
                    f" ({reasons})" if reasons else ""
                )
            sel_genes = selection_gene_table()
            host_passed_symbols = set(
                sel_genes.loc[sel_genes["selection_entity_id"].isin(passed_entities), "Symbol"]
            )
            host_fail_by_symbol = {
                str(row.Symbol): fail_by_entity[str(row.selection_entity_id)]
                for row in sel_genes.itertuples(index=False)
                if str(row.selection_entity_id) in fail_by_entity
            }
            for gene in set(cta_expression) & ((valid_ctas | sex_repro_ctas) - host_passed_symbols):
                drop_reasons[gene] = host_fail_by_symbol.get(gene, "excluded by HostContext")
            valid_ctas &= host_passed_symbols
            sex_repro_ctas &= host_passed_symbols

        outside_panel = set(cta_expression) - strict_ctas - clinical_target_ctas - sex_repro_ctas
        if outside_panel:
            excluded = CTA_excluded_gene_names()
            known = CTA_unfiltered_gene_names()
            rationale = _cta_flag_rationale({gene: cta_expression[gene] for gene in outside_panel})
            for gene in outside_panel:
                if gene in excluded:
                    drop_reasons[gene] = (
                        "excluded from OncoRef's current CTA set with no clinical-target "
                        f"override: {rationale[gene]}"
                    )
                elif gene in known:
                    drop_reasons[gene] = rationale[gene]
                else:
                    drop_reasons[gene] = "not a recognized CTA symbol"

        if min_restriction_confidence is not None:
            confidence_set = {c.upper() for c in min_restriction_confidence}
            confident_ctas = CTA_by_axes(restriction_confidence=confidence_set)
            for gene in set(cta_expression) & (valid_ctas - confident_ctas):
                drop_reasons[gene] = (
                    f"restriction-confidence filter (allowed: {', '.join(sorted(confidence_set))})"
                )
            valid_ctas &= confident_ctas
        if mtec_matrix_path is not None:
            from .mtec import filter_by_mtec, load_mtec_gene_table

            mtec_df = load_mtec_gene_table(mtec_matrix_path)
            candidates = valid_ctas | sex_repro_ctas | clinical_target_ctas
            low_mtec = filter_by_mtec(candidates, mtec_df, threshold=mtec_max_tpm)
            for gene in set(cta_expression) & (candidates - low_mtec):
                drop_reasons[gene] = (
                    f"mTEC filter (requires measured mean TPM <= {mtec_max_tpm:g}; "
                    "above threshold or missing measurement)"
                )
            valid_ctas &= low_mtec
            sex_repro_ctas &= low_mtec
            clinical_target_ctas &= low_mtec

        for gene in set(cta_expression) & (valid_ctas | sex_repro_ctas | clinical_target_ctas):
            tpm = cta_expression[gene]
            if not pd.isna(tpm) and tpm < min_cta_tpm:
                drop_reasons[gene] = f"tumor TPM {tpm:g} below minimum {min_cta_tpm:g}"
        for gene, reason in sorted(drop_reasons.items()):
            warnings.warn(
                f"--cta {gene} dropped: {reason}",
                UserWarning,
                stacklevel=2,
            )

        expressed_ctas = {
            gene: tpm
            for gene, tpm in cta_expression.items()
            if gene in valid_ctas and (pd.isna(tpm) or tpm >= min_cta_tpm)
        }
        expressed_sex_repro = {
            gene: tpm
            for gene, tpm in cta_expression.items()
            if gene in sex_repro_ctas and (pd.isna(tpm) or tpm >= min_cta_tpm)
        }
        flagged_ctas = {
            gene: tpm
            for gene, tpm in cta_expression.items()
            if gene in clinical_target_ctas and (pd.isna(tpm) or tpm >= min_cta_tpm)
        }

        if expressed_ctas:
            _report(
                f"Generating exclusivity-screened CTA peptides for "
                f"{len(expressed_ctas)} gene(s): {', '.join(sorted(expressed_ctas))}..."
            )
            all_cta_peps = cta_exclusive_peptides(
                ensembl_release=resolved_ensembl_release,
                lengths=lengths,
                on_progress=_report if show_progress else None,
                progress_bar=show_progress,
            )
            cta_peps = all_cta_peps[all_cta_peps["gene_name"].isin(expressed_ctas)].copy()
            if not cta_peps.empty:
                cta_peps["source_tpm"] = cta_peps["gene_name"].map(expressed_ctas)
                frames.append(
                    pd.DataFrame(
                        {
                            "peptide": cta_peps["peptide"],
                            "length": cta_peps["length"],
                            "category": "cta",
                            "source": cta_peps["gene_name"],
                            "source_detail": cta_peps["gene_id"],
                            "source_tpm": cta_peps["source_tpm"],
                        }
                    )
                )

        if expressed_sex_repro:
            from .cta_selection import selection_gene_table
            from .partition import CTA_partition_gene_ids
            from .peptides import _reference_overlapping_peptides, gene_peptides

            _report(
                f"Generating exclusivity-screened sex-restricted reproductive antigen peptides "
                f"for {len(expressed_sex_repro)} gene(s): {', '.join(sorted(expressed_sex_repro))}..."
            )
            sel_table = selection_gene_table()
            admitted_ids = set(
                sel_table.loc[sel_table["Symbol"].isin(expressed_sex_repro), "Ensembl_Gene_ID"]
            )
            if admitted_ids:
                sex_peps = gene_peptides(
                    admitted_ids,
                    ensembl_release=resolved_ensembl_release,
                    lengths=lengths,
                )
                if not sex_peps.empty:
                    # Stream background coding genes (excluding strict CTAs and the
                    # admitted sex-restricted reproductive genes themselves).
                    partition = CTA_partition_gene_ids(ensembl_release=resolved_ensembl_release)
                    bg_gene_ids = sorted(partition.non_cta - admitted_ids)
                    overlapping = _reference_overlapping_peptides(
                        set(sex_peps["peptide"]),
                        ensembl_release=resolved_ensembl_release,
                        lengths=lengths,
                        gene_ids=bg_gene_ids,
                    )
                    sex_peps = sex_peps[~sex_peps["peptide"].isin(overlapping)].copy()
                if not sex_peps.empty:
                    sex_peps["source_tpm"] = sex_peps["gene_name"].map(expressed_sex_repro)
                    frames.append(
                        pd.DataFrame(
                            {
                                "peptide": sex_peps["peptide"],
                                "length": sex_peps["length"],
                                "category": "cta",
                                "source": sex_peps["gene_name"],
                                "source_detail": sex_peps["gene_id"],
                                "source_tpm": sex_peps["source_tpm"],
                            }
                        )
                    )

        if flagged_ctas:
            # A flagged gene is, by construction, part of the "non-CTA"
            # background that cta_exclusive_peptides() screens strict CTAs
            # against (CTA_partition_gene_ids puts anything outside the
            # strict set into non_cta, including a gene excluded only for
            # an expression-safety reason like CTAG2's heart signal), and
            # cta_peptides()'s own gene-name resolver gates on that same
            # strict set -- neither can be reused as-is for a flagged gene
            # without either zeroing out its own peptides against itself
            # or resolving to no genes at all. Generate directly instead,
            # and say so via flag_reason: a flagged row has not been
            # checked for sequence overlap with other proteins the way a
            # strict CTA has, on top of its own oncoref exclusion caveat.
            _report(
                f"Generating CTA peptides for {len(flagged_ctas)} flagged clinical-target "
                f"gene(s): {', '.join(sorted(flagged_ctas))}..."
            )
            flagged_peps = _cta_flagged_gene_peptides(
                flagged_ctas, ensembl_release=resolved_ensembl_release, lengths=lengths
            )
            if not flagged_peps.empty:
                flagged_peps["source_tpm"] = flagged_peps["gene_name"].map(flagged_ctas)
                rationale_by_gene = _cta_flag_rationale(flagged_ctas)
                frames.append(
                    pd.DataFrame(
                        {
                            "peptide": flagged_peps["peptide"],
                            "length": flagged_peps["length"],
                            "category": "cta_flagged",
                            "source": flagged_peps["gene_name"],
                            "source_detail": flagged_peps["gene_id"],
                            "source_tpm": flagged_peps["source_tpm"],
                            "flag_reason": flagged_peps["gene_name"].map(rationale_by_gene)
                            + " (not screened for peptide overlap with other proteins the "
                            "way a strict CTA is)",
                        }
                    )
                )

    # ── Mutant peptides ─────────────────────────────────────────────────
    if mutations:
        from .mutations import HOTSPOT_MUTATIONS
        from .mutations import mutant_peptides as _mutant_peptides

        mutation_labels = set(mutations)
        matched = [m for m in HOTSPOT_MUTATIONS if m["label"] in mutation_labels]
        if matched:
            _report(
                f"Generating mutant-spanning peptides for {len(matched)} mutation(s): "
                f"{', '.join(m['label'] for m in matched)}..."
            )
            mdf = _mutant_peptides(
                mutations=matched, lengths=lengths, ensembl_release=resolved_ensembl_release
            )
            if not mdf.empty:
                frames.append(
                    pd.DataFrame(
                        {
                            "peptide": mdf["peptide"],
                            "length": mdf["length"],
                            "category": "mutant",
                            "source": mdf["label"],
                            "source_detail": mdf["mutation"],
                            "source_tpm": float("nan"),
                        }
                    )
                )

    # ── Viral peptides ──────────────────────────────────────────────────
    if viruses:
        if require_human_exclusive_viral:
            if ctx.is_human_host and host_proteins is None:
                from .viral import human_exclusive_viral_peptides as _viral_peps
            else:
                from .viral import host_exclusive_viral_peptides as _viral_peps
        else:
            from .viral import viral_peptides as _viral_peps

        for vk in viruses:
            _report(f"Generating viral peptides for {vk}...")
            if require_human_exclusive_viral:
                if ctx.is_human_host and host_proteins is None:
                    vdf = _viral_peps(
                        virus=vk, lengths=lengths, ensembl_release=resolved_ensembl_release
                    )
                else:
                    vdf = _viral_peps(
                        virus=vk,
                        lengths=lengths,
                        ensembl_release=resolved_ensembl_release,
                        host_species=host_species,
                        host_proteins=host_proteins,
                    )
            else:
                vdf = _viral_peps(virus=vk, lengths=lengths)
            if not vdf.empty:
                frames.append(
                    pd.DataFrame(
                        {
                            "peptide": vdf["peptide"],
                            "length": vdf["length"],
                            "category": "viral",
                            "source": vdf["virus"],
                            "source_detail": vdf["protein_id"],
                            "source_tpm": float("nan"),
                        }
                    )
                )

    # ── Fusion junction peptides ────────────────────────────────────────
    if fusions is not None:
        from .fusions import fusion_junction_peptides

        _report("Generating junction-spanning peptides for supplied fusion breakpoints...")
        fdf = fusion_junction_peptides(
            fusions,
            lengths=lengths,
            ensembl_release=resolved_ensembl_release,
            require_host_exclusive=require_host_exclusive_fusions,
            host_species=host_species,
            host_proteins=host_proteins,
        )
        if not fdf.empty:
            frames.append(
                pd.DataFrame(
                    {
                        "peptide": fdf["peptide"],
                        "length": fdf["length"],
                        "category": "fusion",
                        "source": fdf["label"],
                        "source_detail": fdf["gene_5prime"].astype(str)
                        + "-"
                        + fdf["gene_3prime"].astype(str),
                        "source_tpm": float("nan"),
                    }
                )
            )

    if not frames:
        return _with_selection_provenance(
            pd.DataFrame(columns=output_columns),
            resolved_ensembl_release,
            species_context=ctx,
            allele_capabilities=allele_caps,
            ms_evidence_scope=ms_evidence_scope,
            patient_sex=norm_patient_sex,
            sex_eligibility=norm_sex_eligibility if norm_patient_sex is not None else None,
            host_context=resolved_host,
            near_self_meta=(
                {"near_self_1mm_screened": True, "near_self_background": near_self_background}
                if run_near_self_screen
                else None
            ),
        )

    combined = pd.concat(frames, ignore_index=True)
    _report(f"Generated {len(combined)} candidate peptide rows.")

    if proteoform_rollup and ctx.is_human_host:
        before = len(combined)
        combined = _apply_proteoform_rollup(combined)
        if len(combined) != before:
            _report(
                f"Rolled identical-protein CTA paralogs into groups: "
                f"{before} -> {len(combined)} rows."
            )

    # ── IEDB/CEDAR evidence ─────────────────────────────────────────────
    if not skip_ms_evidence:
        _report(
            f"Looking up public MS evidence for {combined['peptide'].nunique()} "
            "unique peptide(s)..."
        )
    combined, raw_ms_hits = _attach_ms_evidence(
        combined,
        iedb_path=iedb_path,
        cedar_path=cedar_path,
        evidence_package=evidence_package,
        mhc_class=mhc_class,
        skip_ms_evidence=skip_ms_evidence,
        requested_alleles=hla_alleles,
        species_context=ctx,
        ms_evidence_scope=ms_evidence_scope,
        observations=ms_observations,
    )

    if enforce_tumor_specificity and "ms_in_healthy_tissue" in combined.columns:
        safe_mask = ~combined["ms_in_healthy_tissue"].astype(bool)
        # Match Hitlist's sex-specific normal reproductive tract MS flags
        # (src_healthy_reproductive_female / src_healthy_reproductive_male)
        # against the patient's biological sex and eligibility policy:
        # under all_adult, exclude same-sex adult reproductive tract MS hits
        # while allowing opposite-sex hits; under sex_organ_only, retain
        # same-sex reproductive MS hits with explicit flags.
        if (
            norm_patient_sex == "female"
            and norm_sex_eligibility == "all_adult"
            and "ms_in_healthy_reproductive_female" in combined.columns
        ):
            safe_mask = safe_mask & ~combined["ms_in_healthy_reproductive_female"].astype(bool)
        elif (
            norm_patient_sex == "male"
            and norm_sex_eligibility == "all_adult"
            and "ms_in_healthy_reproductive_male" in combined.columns
        ):
            safe_mask = safe_mask & ~combined["ms_in_healthy_reproductive_male"].astype(bool)
        combined = combined[safe_mask].reset_index(drop=True)

    if combined.empty:
        return _with_selection_provenance(
            pd.DataFrame(columns=output_columns),
            resolved_ensembl_release,
            species_context=ctx,
            allele_capabilities=allele_caps,
            ms_evidence_scope=ms_evidence_scope,
            patient_sex=norm_patient_sex,
            sex_eligibility=norm_sex_eligibility if norm_patient_sex is not None else None,
            host_context=resolved_host,
            near_self_meta=(
                {"near_self_1mm_screened": True, "near_self_background": near_self_background}
                if run_near_self_screen
                else None
            ),
        )

    # ── Topiary presentation scoring ────────────────────────────────────
    if score_presentation and hla_alleles:
        from .scoring import score_presentation as _score

        unique_peps = combined["peptide"].unique().tolist()
        _report(
            f"Scoring {len(unique_peps)} peptide(s) against {len(hla_alleles)} "
            f"HLA allele(s) via {predictor}..."
        )
        scores = _score(peptides=unique_peps, alleles=hla_alleles, predictor=predictor)
        combined = _merge_best_allele(combined, scores)
        # Re-attribute predicted restriction explicitly only for positive MS rows
        # without measured exact or donor-set restriction, keeping prediction
        # separate from measured MHC restriction (#216).
        if not skip_ms_evidence and raw_ms_hits is not None and not raw_ms_hits.empty:
            combined = _apply_predicted_restriction_attribution(
                combined,
                raw_ms_hits,
                requested_alleles=hla_alleles,
                species_context=ctx,
                ms_evidence_scope=ms_evidence_scope,
            )
    else:
        for col in ("best_allele", "presentation_percentile", "presentation_score", "affinity_nm"):
            if col not in combined.columns:
                combined[col] = pd.NA

    # ── Tier assignment, detectability & OncoRef safety caveats ─────────
    from .ms_evidence import annotate_ms_detectability
    from .near_self import annotate_cta_safety_caveats

    combined = _assign_tiers(combined)
    combined = annotate_ms_detectability(combined)
    combined = annotate_cta_safety_caveats(
        combined,
        source_col="source",
        peptide_col="peptide",
        category_col="category",
        species=ctx.host_species,
    )
    if drop_weak_tier:
        combined = combined[combined["tier"] < 4].reset_index(drop=True)

    near_self_meta: dict[str, object] | None = None
    if run_near_self_screen:
        from .near_self import annotate_near_self_1mm

        _report(
            f"Screening {combined['peptide'].nunique()} peptide(s) for 1-aa near-self matches "
            f"against OncoRef {near_self_background} proteome..."
        )
        combined = annotate_near_self_1mm(
            combined,
            peptide_col="peptide",
            ensembl_release=resolved_ensembl_release,
            species=ctx.host_species,
            background=near_self_background,
            background_proteins=near_self_background_proteins,
            max_matches=max_near_self_1mm_matches,
        )
        near_self_meta = {
            "near_self_1mm_screened": True,
            "near_self_background": near_self_background,
            "near_self_max_matches": max_near_self_1mm_matches,
            "near_self_contract_id": combined.attrs.get("near_self_contract_id"),
        }

    # ── Deterministic sort ──────────────────────────────────────────────
    sort_cols = ["tier", "ms_hit_count", "presentation_percentile", "peptide", "best_allele"]
    sort_ascending = [True, False, True, True, True]
    available = [c for c in sort_cols if c in combined.columns]
    ascending = [a for c, a in zip(sort_cols, sort_ascending) if c in combined.columns]
    combined = combined.sort_values(
        available, ascending=ascending, na_position="last", kind="mergesort"
    ).reset_index(drop=True)

    for col in output_columns:
        if col not in combined.columns:
            combined[col] = pd.NA
    _report(f"Done: {len(combined)} target row(s).")
    return _with_selection_provenance(
        combined[output_columns],
        resolved_ensembl_release,
        species_context=ctx,
        allele_capabilities=allele_caps,
        ms_evidence_scope=ms_evidence_scope,
        patient_sex=norm_patient_sex,
        sex_eligibility=norm_sex_eligibility if norm_patient_sex is not None else None,
        host_context=resolved_host,
        near_self_meta=near_self_meta,
    )


def _with_selection_provenance(
    frame,
    ensembl_release,
    *,
    species_context=None,
    allele_capabilities=None,
    ms_evidence_scope: str = "all",
    patient_sex: str | None = None,
    sex_eligibility: str | None = None,
    host_context=None,
    near_self_meta: dict[str, object] | None = None,
):
    from .cta_selection import selection_provenance

    if species_context is None or species_context.is_human_host:
        frame.attrs.update(selection_provenance(ensembl_release=ensembl_release))
    else:
        frame.attrs["ensembl_release"] = ensembl_release
    if species_context is not None:
        frame.attrs["species_context"] = species_context.to_provenance()
    if allele_capabilities is not None:
        frame.attrs["allele_prediction_capability"] = allele_capabilities
    frame.attrs["ms_evidence_scope"] = ms_evidence_scope
    if patient_sex is not None:
        frame.attrs["patient_sex"] = patient_sex
        frame.attrs["sex_eligibility"] = sex_eligibility
    if host_context is not None:
        frame.attrs["host_context"] = host_context.as_dict()
    if near_self_meta:
        frame.attrs.update(near_self_meta)
    return frame


# Backward-compatible alias. ``personalized_targets`` is the canonical name (it's
# also exported from the top level, ``from tsarina import personalized_targets``);
# ``personalize`` is retained so existing ``from tsarina.personalize import
# personalize`` callers keep working.
personalize = personalized_targets


def _attach_ms_evidence(
    combined: pd.DataFrame,
    *,
    iedb_path: str | Path | None,
    cedar_path: str | Path | None,
    mhc_class: str,
    skip_ms_evidence: bool,
    evidence_package: str | Path | None = None,
    requested_alleles: list[str] | None = None,
    species_context=None,
    ms_evidence_scope: str = "all",
    observations: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Join restriction-aware MS attribution columns onto the candidate frame (#216)."""
    from .ms_attribution import (
        EVIDENCE_UNAVAILABLE,
        NO_MS_MATCH,
        attribute_ms_evidence,
        summarize_ms_attribution,
    )

    defaults: dict[str, object] = {
        "ms_evidence_status": "not_loaded" if skip_ms_evidence else "no_match",
        "ms_attribution_tier": EVIDENCE_UNAVAILABLE if skip_ms_evidence else NO_MS_MATCH,
        "ms_hit_count": 0,
        "ms_exact_mhc_hits": 0,
        "ms_donor_set_hits": 0,
        "ms_alleles": "",
        "ms_allele_count": 0,
        "ms_in_cancer": False,
        "ms_in_healthy_tissue": False,
        "ms_in_healthy_reproductive_female": False,
        "ms_in_healthy_reproductive_male": False,
    }

    if skip_ms_evidence:
        for col, default in defaults.items():
            combined[col] = default
        return combined, None

    if observations is not None:
        hits = observations.copy()
    else:
        from .ms_evidence import load_public_ms_hits

        hit_kwargs = {
            "peptides": set(combined["peptide"].unique()),
            "iedb_path": iedb_path,
            "cedar_path": cedar_path,
            "mhc_class": mhc_class,
            "classify_source": True,
            "drop_binding_assays": True,
        }
        if (species_context is not None and not species_context.is_human_host) or (
            species_context is not None and species_context.is_cross_species_mhc
        ):
            hit_kwargs["mhc_species"] = species_context.mhc_species
        if evidence_package is not None:
            hit_kwargs["evidence_package"] = evidence_package
        hits = load_public_ms_hits(**hit_kwargs)

    attributed = attribute_ms_evidence(
        combined[["peptide"]].drop_duplicates(),
        hits,
        requested_alleles=requested_alleles,
        species_context=species_context,
        evidence_scope="all",
        evidence_loaded=True,
    )
    summary = summarize_ms_attribution(attributed, group_cols=("peptide",))
    merge_cols = [
        "peptide",
        "ms_evidence_status",
        "ms_attribution_tier",
        "ms_hit_count",
        "ms_exact_mhc_hits",
        "ms_donor_set_hits",
        "ms_alleles",
        "ms_allele_count",
        "ms_in_cancer",
        "ms_in_healthy_tissue",
        "ms_in_healthy_reproductive_female",
        "ms_in_healthy_reproductive_male",
    ]
    # validate="m:1": summary is one row per peptide; fail loudly rather than
    # silently multiply candidate rows if that invariant ever breaks.
    combined = combined.merge(summary[merge_cols], on="peptide", how="left", validate="m:1")

    int_cols = {"ms_hit_count", "ms_exact_mhc_hits", "ms_donor_set_hits", "ms_allele_count"}
    bool_cols = {
        "ms_in_cancer",
        "ms_in_healthy_tissue",
        "ms_in_healthy_reproductive_female",
        "ms_in_healthy_reproductive_male",
    }
    for col, default in defaults.items():
        if col not in combined.columns:
            combined[col] = default
            continue
        if col in int_cols:
            combined[col] = (
                pd.to_numeric(combined[col], errors="coerce").fillna(default).astype(int)
            )
        elif col in bool_cols:
            combined[col] = combined[col].where(combined[col].notna(), default).astype(bool)
        else:
            combined[col] = combined[col].where(combined[col].notna(), default).astype(str)

    if ms_evidence_scope == "any_mhc":
        combined = combined[combined["ms_hit_count"] > 0].reset_index(drop=True)
    elif ms_evidence_scope == "requested_mhc":
        combined = combined[
            combined["ms_attribution_tier"].isin(
                {"exact_requested_mhc", "donor_set_ambiguous", "predicted_restriction"}
            )
        ].reset_index(drop=True)
    return combined, hits


def _apply_predicted_restriction_attribution(
    combined: pd.DataFrame,
    raw_hits: pd.DataFrame,
    *,
    requested_alleles: list[str],
    species_context,
    ms_evidence_scope: str,
) -> pd.DataFrame:
    """Upgrade unresolved positive-MS rows to ``predicted_restriction`` when scored."""
    from .ms_attribution import (
        ANY_MHC_PEPTIDE,
        CLASS_OR_SEROTYPE_ONLY,
        attribute_ms_evidence,
        summarize_ms_attribution,
    )

    t3_cut = PRESENTATION_PERCENTILE_THRESHOLDS[2]
    pct = pd.to_numeric(combined.get("presentation_percentile"), errors="coerce")
    pred_rows = combined.loc[
        pct.le(t3_cut)
        & combined["best_allele"].notna()
        & combined["ms_attribution_tier"].isin({ANY_MHC_PEPTIDE, CLASS_OR_SEROTYPE_ONLY}),
        ["peptide", "best_allele"],
    ].drop_duplicates()
    if pred_rows.empty:
        return combined

    attributed = attribute_ms_evidence(
        pred_rows[["peptide"]].drop_duplicates(),
        raw_hits,
        requested_alleles=requested_alleles,
        predicted_assignments=pred_rows,
        species_context=species_context,
        evidence_scope="all",
        evidence_loaded=True,
    )
    summary = summarize_ms_attribution(attributed, group_cols=("peptide",))
    tier_map = dict(zip(summary["peptide"], summary["ms_attribution_tier"]))
    combined = combined.copy()
    mask = combined["peptide"].isin(tier_map) & combined["ms_attribution_tier"].isin(
        {ANY_MHC_PEPTIDE, CLASS_OR_SEROTYPE_ONLY}
    )
    combined.loc[mask, "ms_attribution_tier"] = combined.loc[mask, "peptide"].map(tier_map)
    if ms_evidence_scope == "requested_mhc":
        combined = combined[
            combined["ms_attribution_tier"].isin(
                {"exact_requested_mhc", "donor_set_ambiguous", "predicted_restriction"}
            )
        ].reset_index(drop=True)
    return combined


def _merge_best_allele(combined: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    """Attach each peptide's lowest-percentile allele + associated score columns."""
    for col in ("best_allele", "presentation_percentile", "presentation_score", "affinity_nm"):
        if col not in combined.columns:
            combined[col] = pd.NA

    if scores.empty or "presentation_percentile" not in scores.columns:
        return combined

    best = (
        scores.sort_values(["presentation_percentile", "peptide", "allele"], kind="mergesort")
        .groupby("peptide", as_index=False)
        .first()
        .rename(columns={"allele": "best_allele"})
    )
    cols = ["peptide", "best_allele", "presentation_percentile"]
    for extra in ("presentation_score", "affinity_nm"):
        if extra in best.columns:
            cols.append(extra)
    combined = combined.drop(
        columns=[
            c
            for c in ("best_allele", "presentation_percentile", "presentation_score", "affinity_nm")
            if c in combined.columns
        ],
        errors="ignore",
    ).merge(best[cols], on="peptide", how="left")
    return combined


def _assign_tiers(combined: pd.DataFrame) -> pd.DataFrame:
    """Assign tier (1-4) + tier_label + tier_reason using explicit gates.

    Uses :data:`tsarina.scoring.PRESENTATION_PERCENTILE_THRESHOLDS` as
    ``(t1_cut, t2_cut, t3_cut)``:

    - **T1 (STRONG)**: ``percentile <= t1_cut`` AND
      (``ms_in_cancer`` OR category in ``{mutant, viral}``).
    - **T2 (MODERATE)**: ``percentile <= t2_cut`` AND
      (``ms_hit_count > 0`` OR category in ``{mutant, viral}``).
    - **T3 (CANDIDATE)**: ``percentile <= t3_cut``.
    - **T4 (WEAK)**: otherwise (including rows with no percentile).
    """
    t1_cut, t2_cut, t3_cut = PRESENTATION_PERCENTILE_THRESHOLDS
    idx = combined.index

    pct = pd.to_numeric(combined.get("presentation_percentile"), errors="coerce")
    ms_cancer = (
        combined["ms_in_cancer"].astype(bool)
        if "ms_in_cancer" in combined.columns
        else pd.Series(False, index=idx)
    )
    ms_hits = pd.to_numeric(combined.get("ms_hit_count", 0), errors="coerce").fillna(0)
    category = combined["category"].astype(str)
    non_self_category = category.isin(("mutant", "viral", "fusion"))

    t1_mask = pct.le(t1_cut) & (ms_cancer | non_self_category)
    t2_mask = pct.le(t2_cut) & ((ms_hits > 0) | non_self_category)
    t3_mask = pct.le(t3_cut)

    tier = pd.Series(4, index=idx, dtype="int64")
    tier = tier.mask(t3_mask, 3)
    tier = tier.mask(t2_mask, 2)
    tier = tier.mask(t1_mask, 1)

    # Vectorized tier_reason assembly — prefix + per-tier driver, chained .mask
    t1_driver = pd.Series(
        np.where(ms_cancer.to_numpy(), "cancer_ms", category.to_numpy()), index=idx
    )
    t2_driver = pd.Series(
        np.where(non_self_category.to_numpy(), category.to_numpy(), "any_ms"), index=idx
    )
    t4_driver = pd.Series(np.where(pct.isna().to_numpy(), "unscored", "below_threshold"), index=idx)
    reasons = t4_driver.copy()
    reasons = reasons.mask(tier == 3, "presentation_only")
    reasons = reasons.mask(tier == 2, "moderate_presentation+" + t2_driver)
    reasons = reasons.mask(tier == 1, "strong_presentation+" + t1_driver)

    combined = combined.copy()
    combined["tier"] = tier
    combined["tier_label"] = combined["tier"].map(_TIER_LABELS)
    combined["tier_reason"] = reasons.to_numpy()
    return combined
