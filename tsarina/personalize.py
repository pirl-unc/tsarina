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
3. **CTA restriction confidence** — CTAs are gated to the HIGH/MODERATE bins
   of ``restriction_confidence``; LOW-confidence CTAs never contribute
   peptides unless the caller explicitly opts in.
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

    from tsarina.personalize import personalize

    targets = personalize(
        hla_alleles=["HLA-A*02:01", "HLA-B*07:02"],
        cta_expression={"MAGEA4": 142.5, "PRAME": 87.3},
        mutations=["KRAS G12D"],
        viruses=["hpv16"],
    )
"""

from __future__ import annotations

import warnings
from collections.abc import Iterable
from pathlib import Path

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
    "ms_hit_count",
    "ms_alleles",
    "ms_allele_count",
    "ms_in_cancer",
    "ms_in_healthy_tissue",
    "best_allele",
    "presentation_percentile",
    "presentation_score",
    "affinity_nm",
    "tier",
    "tier_label",
    "tier_reason",
)


def personalize(
    hla_alleles: list[str],
    cta_expression: dict[str, float] | None = None,
    mutations: list[str] | None = None,
    viruses: list[str] | None = None,
    lengths: tuple[int, ...] = (8, 9, 10, 11),
    ensembl_release: int = 112,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    mhc_class: str = "I",
    min_cta_tpm: float = 2.0,
    min_restriction_confidence: Iterable[str] | None = ("HIGH", "MODERATE"),
    mtec_matrix_path: str | Path | None = None,
    mtec_max_tpm: float = 1.0,
    require_human_exclusive_viral: bool = True,
    enforce_tumor_specificity: bool = True,
    score_presentation: bool = True,
    skip_ms_evidence: bool = False,
    predictor: str = "mhcflurry",
    drop_weak_tier: bool = True,
) -> pd.DataFrame:
    """Build a personalized, tier-ranked target list for a single patient.

    Parameters
    ----------
    hla_alleles
        Patient's HLA class I alleles (e.g. ``["HLA-A*02:01", "HLA-B*07:02"]``).
    cta_expression
        Dict mapping CTA gene symbol to tumor RNA expression in TPM.  Only
        genes passing all CTA gates (membership in the filtered CTA set,
        ``tpm >= min_cta_tpm``, ``restriction_confidence`` allowed, optional
        mTEC cutoff) contribute peptides.  Peptides are restricted to those
        exclusive to CTA proteins (not present in any non-CTA protein).
    mutations
        List of hotspot labels (e.g. ``["KRAS G12D"]``).  Wildtype-identical
        k-mers are filtered in :mod:`tsarina.mutations`; they are already
        tumor-specific by construction.
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
        Minimum CTA expression in TPM (default 2.0 — matches the
        ``never_expressed`` HPA cutoff).
    min_restriction_confidence
        Allowed ``restriction_confidence`` bins for CTAs.  Default
        ``("HIGH", "MODERATE")``; pass ``None`` to disable this gate.
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

    Returns
    -------
    pd.DataFrame
        Tier-ranked target list.  Columns:

        - ``peptide``, ``length``, ``category``, ``source``,
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
    frames: list[pd.DataFrame] = []

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
    if cta_expression:
        from .gene_sets import CTA_by_axes, CTA_gene_names
        from .peptides import cta_exclusive_peptides

        valid_ctas = CTA_gene_names()
        if min_restriction_confidence is not None:
            confidence_set = {c.upper() for c in min_restriction_confidence}
            confident_ctas = CTA_by_axes(restriction_confidence=confidence_set)
            valid_ctas = valid_ctas & confident_ctas
        if mtec_matrix_path is not None:
            from .mtec import filter_by_mtec, load_mtec_gene_table

            mtec_df = load_mtec_gene_table(mtec_matrix_path)
            valid_ctas = filter_by_mtec(valid_ctas, mtec_df, threshold=mtec_max_tpm)

        expressed_ctas = {
            gene: tpm
            for gene, tpm in cta_expression.items()
            if gene in valid_ctas and tpm >= min_cta_tpm
        }
        if expressed_ctas:
            all_cta_peps = cta_exclusive_peptides(ensembl_release=ensembl_release, lengths=lengths)
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

    # ── Mutant peptides ─────────────────────────────────────────────────
    if mutations:
        from .mutations import HOTSPOT_MUTATIONS
        from .mutations import mutant_peptides as _mutant_peptides

        mutation_labels = set(mutations)
        matched = [m for m in HOTSPOT_MUTATIONS if m["label"] in mutation_labels]
        if matched:
            mdf = _mutant_peptides(
                mutations=matched, lengths=lengths, ensembl_release=ensembl_release
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
            from .viral import human_exclusive_viral_peptides as _viral_peps
        else:
            from .viral import viral_peptides as _viral_peps

        for vk in viruses:
            if require_human_exclusive_viral:
                vdf = _viral_peps(virus=vk, lengths=lengths, ensembl_release=ensembl_release)
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

    if not frames:
        return pd.DataFrame(columns=list(_OUTPUT_COLUMNS))

    combined = pd.concat(frames, ignore_index=True)

    # ── IEDB/CEDAR evidence ─────────────────────────────────────────────
    combined = _attach_ms_evidence(
        combined,
        iedb_path=iedb_path,
        cedar_path=cedar_path,
        mhc_class=mhc_class,
        skip_ms_evidence=skip_ms_evidence,
    )

    if enforce_tumor_specificity and "ms_in_healthy_tissue" in combined.columns:
        safe_mask = ~combined["ms_in_healthy_tissue"].astype(bool)
        combined = combined[safe_mask].reset_index(drop=True)

    if combined.empty:
        return pd.DataFrame(columns=list(_OUTPUT_COLUMNS))

    # ── Topiary presentation scoring ────────────────────────────────────
    if score_presentation and hla_alleles:
        from .scoring import score_presentation as _score

        unique_peps = combined["peptide"].unique().tolist()
        scores = _score(peptides=unique_peps, alleles=hla_alleles, predictor=predictor)
        combined = _merge_best_allele(combined, scores)
    else:
        for col in ("best_allele", "presentation_percentile", "presentation_score", "affinity_nm"):
            if col not in combined.columns:
                combined[col] = pd.NA

    # ── Tier assignment ─────────────────────────────────────────────────
    combined = _assign_tiers(combined)
    if drop_weak_tier:
        combined = combined[combined["tier"] < 4].reset_index(drop=True)

    # ── Deterministic sort ──────────────────────────────────────────────
    sort_cols = ["tier", "ms_hit_count", "presentation_percentile", "peptide", "best_allele"]
    sort_ascending = [True, False, True, True, True]
    available = [c for c in sort_cols if c in combined.columns]
    ascending = [a for c, a in zip(sort_cols, sort_ascending) if c in combined.columns]
    combined = combined.sort_values(
        available, ascending=ascending, na_position="last", kind="mergesort"
    ).reset_index(drop=True)

    for col in _OUTPUT_COLUMNS:
        if col not in combined.columns:
            combined[col] = pd.NA
    return combined[list(_OUTPUT_COLUMNS)]


def _attach_ms_evidence(
    combined: pd.DataFrame,
    *,
    iedb_path: str | Path | None,
    cedar_path: str | Path | None,
    mhc_class: str,
    skip_ms_evidence: bool,
) -> pd.DataFrame:
    """Join MS-aggregated columns onto the candidate frame."""
    defaults: dict[str, object] = {
        "ms_hit_count": 0,
        "ms_alleles": "",
        "ms_allele_count": 0,
        "ms_in_cancer": False,
        "ms_in_healthy_tissue": False,
    }

    if skip_ms_evidence:
        for col, default in defaults.items():
            combined[col] = default
        return combined

    if iedb_path is not None or cedar_path is not None:
        from .datasources import resolve_dataset_paths
        from .iedb import scan_public_ms

        resolved_iedb, resolved_cedar = resolve_dataset_paths(iedb_path, cedar_path)
        hits = scan_public_ms(
            peptides=set(combined["peptide"].unique()),
            iedb_path=resolved_iedb,
            cedar_path=resolved_cedar,
            mhc_class=mhc_class,
            classify_source=True,
        )
    else:
        from .indexing import load_ms_evidence

        hits = load_ms_evidence(
            peptides=set(combined["peptide"].unique()),
            mhc_class=mhc_class,
            mhc_species="Homo sapiens",
        )

    if hits.empty:
        for col, default in defaults.items():
            combined[col] = default
        return combined

    agg_cols: dict[str, tuple] = {
        "ms_hit_count": ("peptide", "size"),
        "ms_alleles": ("mhc_restriction", lambda x: ";".join(sorted(set(x)))),
        "ms_allele_count": ("mhc_restriction", lambda x: len(set(x))),
    }
    if "src_cancer" in hits.columns:
        agg_cols["ms_in_cancer"] = ("src_cancer", "any")
    if "src_healthy_tissue" in hits.columns:
        agg_cols["ms_in_healthy_tissue"] = ("src_healthy_tissue", "any")

    hit_agg = hits.groupby("peptide", as_index=False).agg(**agg_cols)
    combined = combined.merge(hit_agg, on="peptide", how="left")

    int_cols = {"ms_hit_count", "ms_allele_count"}
    bool_cols = {"ms_in_cancer", "ms_in_healthy_tissue"}
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

    pct = pd.to_numeric(combined.get("presentation_percentile"), errors="coerce")
    ms_cancer = combined.get("ms_in_cancer", False)
    ms_hits = pd.to_numeric(combined.get("ms_hit_count", 0), errors="coerce").fillna(0)
    category = combined["category"]
    non_self_category = category.isin(("mutant", "viral"))

    t1_mask = pct.le(t1_cut) & (ms_cancer.astype(bool) | non_self_category)
    t2_mask = pct.le(t2_cut) & ((ms_hits > 0) | non_self_category)
    t3_mask = pct.le(t3_cut)

    tier = pd.Series(4, index=combined.index, dtype="int64")
    tier = tier.mask(t3_mask, 3)
    tier = tier.mask(t2_mask, 2)
    tier = tier.mask(t1_mask, 1)

    reasons = []
    for idx in combined.index:
        t = int(tier.loc[idx])
        cat = category.loc[idx]
        if t == 1:
            driver = (
                "cancer_ms"
                if bool(ms_cancer.loc[idx] if hasattr(ms_cancer, "loc") else False)
                else cat
            )
            reasons.append(f"strong_presentation+{driver}")
        elif t == 2:
            if cat in ("mutant", "viral"):
                reasons.append(f"moderate_presentation+{cat}")
            else:
                reasons.append("moderate_presentation+any_ms")
        elif t == 3:
            reasons.append("presentation_only")
        else:
            reasons.append("unscored" if pd.isna(pct.loc[idx]) else "below_threshold")

    combined = combined.copy()
    combined["tier"] = tier
    combined["tier_label"] = combined["tier"].map(_TIER_LABELS)
    combined["tier_reason"] = reasons
    return combined
