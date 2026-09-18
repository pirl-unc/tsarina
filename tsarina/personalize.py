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
   peptides unless the caller explicitly opts in. A gene oncoref excludes
   from the strict set entirely but still tracks as a known clinical
   target (e.g. CTAG2/LAGE-1, excluded for a low-level HPA heart RNA
   signal) is not simply dropped: it's included under
   ``category="cta_flagged"`` with the exclusion reason in ``flag_reason``,
   so a caller who named the gene explicitly sees it and its caveat rather
   than silence. A ``--cta`` gene that is neither a recognized CTA nor a
   known clinical target (a typo, or a gene with no CTA evidence at all)
   is dropped with a warning naming it, rather than the same silence.
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

import warnings
from collections.abc import Iterable
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


def _cta_flagged_gene_peptides(
    gene_names: Iterable[str], *, ensembl_release: int, lengths: tuple[int, ...]
) -> pd.DataFrame:
    """Generate peptides for a clinical-target CTA that oncoref excludes
    from the strict CTA gene universe (CTAG2/LAGE-1's heart-signal
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
    """oncoref's ``specificity_rationale`` for each flagged gene, for the
    ``flag_reason`` column -- the actual reason a caller should read before
    trusting a "cta_flagged" row, not just that one exists."""
    if not genes:
        return {}
    from oncoref import cta as _oncoref_cta

    df = _oncoref_cta.cta_df()
    rows = df[df["Symbol"].isin(genes) & df["specificity_rationale"].notna()]
    rationale = dict(zip(rows["Symbol"], rows["specificity_rationale"]))
    return {
        gene: rationale.get(gene, "excluded from the strict default CTA set by oncoref")
        for gene in genes
    }


def personalized_targets(
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
        genes passing all CTA gates (membership in the filter-passing CTA set,
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
        from .gene_sets import CTA_by_axes, CTA_clinical_target_gene_names, CTA_gene_names
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

        # Genes with real CTA-source evidence that oncoref nonetheless
        # excludes from the strict default set (e.g. CTAG2/LAGE-1: a
        # low-level HPA heart RNA signal) but keeps as a known clinical
        # target. A caller who explicitly named the gene should see it and
        # its caveat, not have it silently vanish the way a typo would.
        clinical_target_ctas = CTA_clinical_target_gene_names() - valid_ctas

        unrecognized = sorted(
            gene
            for gene in cta_expression
            if gene not in valid_ctas and gene not in clinical_target_ctas
        )
        if unrecognized:
            warnings.warn(
                "--cta gene(s) not in the current CTA panel, dropped: "
                + ", ".join(unrecognized)
                + " (not a recognized CTA symbol, or excluded with no clinical-target "
                "override -- see tsarina.gene_sets.CTA_excluded_gene_names for why)",
                UserWarning,
                stacklevel=2,
            )

        expressed_ctas = {
            gene: tpm
            for gene, tpm in cta_expression.items()
            if gene in valid_ctas and tpm >= min_cta_tpm
        }
        flagged_ctas = {
            gene: tpm
            for gene, tpm in cta_expression.items()
            if gene in clinical_target_ctas and tpm >= min_cta_tpm
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
            flagged_peps = _cta_flagged_gene_peptides(
                flagged_ctas, ensembl_release=ensembl_release, lengths=lengths
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

    from .ms_evidence import aggregate_ms_hits_by_peptide, load_public_ms_hits

    hits = load_public_ms_hits(
        peptides=set(combined["peptide"].unique()),
        iedb_path=iedb_path,
        cedar_path=cedar_path,
        mhc_class=mhc_class,
        classify_source=True,
        drop_binding_assays=True,
    )

    if hits.empty:
        for col, default in defaults.items():
            combined[col] = default
        return combined

    hit_agg = aggregate_ms_hits_by_peptide(
        hits,
        source_flag_outputs={
            "src_cancer": "ms_in_cancer",
            "src_healthy_tissue": "ms_in_healthy_tissue",
        },
        cell_line_output=None,
    )
    # validate="m:1": hit_agg is one row per peptide; fail loudly rather than
    # silently multiply candidate rows if that invariant ever breaks.
    combined = combined.merge(hit_agg, on="peptide", how="left", validate="m:1")

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
    idx = combined.index

    pct = pd.to_numeric(combined.get("presentation_percentile"), errors="coerce")
    ms_cancer = (
        combined["ms_in_cancer"].astype(bool)
        if "ms_in_cancer" in combined.columns
        else pd.Series(False, index=idx)
    )
    ms_hits = pd.to_numeric(combined.get("ms_hit_count", 0), errors="coerce").fillna(0)
    category = combined["category"].astype(str)
    non_self_category = category.isin(("mutant", "viral"))

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
