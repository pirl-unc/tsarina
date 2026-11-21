"""Deterministic marginal-gain allocation, target-specific HLA support, and tradeoff audits."""

from __future__ import annotations

from dataclasses import replace

import pandas as pd

from .vaccine_construct import construct_fits


def _expression_tables(cancer):
    prevalence, weights = {}, {}
    for row in cancer.itertuples(index=False):
        if row.complete_measurement:
            prevalence.setdefault(row.proteoform_key, {})[row.burden_category] = row.prevalence_p95
        weights[row.burden_category] = (
            row.world_mortality_pct / 100,
            row.world_incidence_pct / 100,
        )
    return prevalence, weights


def _segment_tier_alleles(segments, ligand_rows=None):
    """Map each segment_id to per-evidence-tier allele sets without conflating tiers (#201/#202)."""
    by_seg = {
        s.segment_id: {
            "monoallelic_ms": set(),
            "sample_allele_ms": set(),
            "unrestricted_ms": set(),
        }
        for s in segments
    }
    if ligand_rows:
        for row in ligand_rows:
            sid = row.get("segment_id")
            tier = row.get("evidence_tier")
            allele = row.get("allele")
            if sid in by_seg and tier in by_seg[sid] and allele:
                by_seg[sid][tier].add(allele)
    return by_seg


def _evaluate_segment_set(selected, prevalence, weights, tier_by_seg, min_padding):
    """Compute conservative expression bounds and measured/inferred evidence totals."""
    current = {}
    global_alleles = set()
    target_alleles = set()
    mono_target_alleles = set()
    peptides = set()
    proteins = set()
    native_aa = 0
    for segment in selected:
        proteins.add(segment.proteoform_key)
        native_aa += len(segment.sequence(min_padding, min_padding))
        global_alleles.update(segment.alleles)
        target_alleles.update((segment.proteoform_key, a) for a in segment.alleles)
        mono_target_alleles.update(
            (segment.proteoform_key, a)
            for a in tier_by_seg.get(segment.segment_id, {}).get("monoallelic_ms", ())
        )
        peptides.update(segment.peptides)
        for category, value in prevalence.get(segment.proteoform_key, {}).items():
            current[category] = max(value, current.get(category, 0.0))
    mortality = sum(weights[c][0] * v for c, v in current.items())
    incidence = sum(weights[c][1] * v for c, v in current.items())
    return {
        "n_proteins": len(proteins),
        "n_segments": len(selected),
        "native_aa": native_aa,
        "mortality_lower_pp": mortality * 100.0,
        "incidence_lower_pp": incidence * 100.0,
        "distinct_ms_peptides": len(peptides),
        "global_panel_alleles": len(global_alleles),
        "target_allele_pairs": len(target_alleles),
        "monoallelic_target_alleles": len(mono_target_alleles),
    }


def select_budget_segments(segments, cancer, config, *, ligand_rows=None):
    """Greedy whole-piece allocation under a total sequence budget.

    Supports three explicit, documented objectives via
    ``config.allocation_objective``:

    - ``"global_union"`` (default): historical lexicographic order — mortality-
      weighted lower expression bound gain per aa, incidence gain per aa,
      new global panel alleles per aa, distinct observed MS peptides per aa,
      already-selected protein, longer contiguous piece, rank.
    - ``"target_specific"`` (#201): credits new target-specific
      ``(proteoform_key, allele)`` pairs (separating measured monoallelic MS
      from inferred alleles) so additional segments of an already-selected CTA
      can earn HLA support on alleles not yet presenting *that* CTA.
    - ``"evidence_tradeoff"`` (#202): applies ``config.min_mortality_gain_pp``
      so marginal mortality-expression gains below the configured percentage-
      point floor do not automatically outrank MS-rich and target-HLA-rich
      regions of common CTAs.

    ``max(p)`` is the conservative union bound; no patient independence or
    clinical coverage product is assumed. Missing measurements contribute no
    *known* gain. Additional same-protein pieces can add evidence but never
    count expression twice.
    """
    prevalence, weights = _expression_tables(cancer)
    tier_by_seg = _segment_tier_alleles(segments, ligand_rows)
    objective = getattr(config, "allocation_objective", "global_union")
    min_gain_frac = float(getattr(config, "min_mortality_gain_pp", 0.0)) / 100.0

    chosen, remaining, history = [], list(segments), []
    current, alleles, target_alleles, mono_target_alleles, peptides, proteins = (
        {},
        set(),
        set(),
        set(),
        set(),
        set(),
    )
    while remaining:
        options = []
        for segment in remaining:
            placement = (segment, config.min_padding, config.min_padding, "")
            if not construct_fits([*chosen, placement], config):
                continue
            gains = [0.0, 0.0]
            for category, value in prevalence.get(segment.proteoform_key, {}).items():
                delta = max(0, value - current.get(category, 0))
                for i in range(2):
                    gains[i] += weights[category][i] * delta
            new_alleles = len(set(segment.alleles) - alleles)
            seg_target_alleles = {(segment.proteoform_key, a) for a in segment.alleles}
            new_target_alleles = len(seg_target_alleles - target_alleles)
            seg_mono = {
                (segment.proteoform_key, a)
                for a in tier_by_seg.get(segment.segment_id, {}).get("monoallelic_ms", ())
            }
            new_mono_target_alleles = len(seg_mono - mono_target_alleles)
            new_peptides = len(set(segment.peptides) - peptides)
            length = len(segment.sequence(config.min_padding, config.min_padding))
            if objective == "global_union":
                if not any(gains) and not new_alleles and not new_peptides:
                    continue
            elif not any(gains) and not new_alleles and not new_target_alleles and not new_peptides:
                continue
            if objective == "target_specific":
                score = (
                    gains[0] / length,
                    gains[1] / length,
                    new_mono_target_alleles / length,
                    new_target_alleles / length,
                    new_peptides / length,
                    new_alleles / length,
                    segment.proteoform_key in proteins,
                    length,
                    -segment.rank,
                )
            elif objective == "evidence_tradeoff":
                effective_mortality = gains[0] if gains[0] >= min_gain_frac else 0.0
                score = (
                    effective_mortality / length,
                    new_mono_target_alleles / length,
                    new_target_alleles / length,
                    new_peptides / length,
                    gains[0] / length,
                    gains[1] / length,
                    segment.proteoform_key in proteins,
                    length,
                    -segment.rank,
                )
            else:
                score = (
                    *[gain / length for gain in gains],
                    new_alleles / length,
                    new_peptides / length,
                    segment.proteoform_key in proteins,
                    length,
                    -segment.rank,
                )
            options.append(
                (
                    score,
                    segment.segment_id,
                    placement,
                    gains,
                    new_alleles,
                    new_target_alleles,
                    new_mono_target_alleles,
                    new_peptides,
                )
            )
        if not options:
            break
        (
            score,
            _,
            placement,
            gains,
            new_alleles,
            new_target_alleles,
            new_mono_target_alleles,
            new_peptides,
        ) = max(options, key=lambda x: (x[0], x[1]))
        segment = placement[0]
        new_proteoform = segment.proteoform_key not in proteins
        chosen.append(placement)
        remaining.remove(segment)
        proteins.add(segment.proteoform_key)
        alleles.update(segment.alleles)
        target_alleles.update((segment.proteoform_key, a) for a in segment.alleles)
        mono_target_alleles.update(
            (segment.proteoform_key, a)
            for a in tier_by_seg.get(segment.segment_id, {}).get("monoallelic_ms", ())
        )
        peptides.update(segment.peptides)
        for category, value in prevalence.get(segment.proteoform_key, {}).items():
            current[category] = max(value, current.get(category, 0))
        history.append(
            {
                "step": len(chosen),
                "segment_id": segment.segment_id,
                "name": segment.name,
                "allocation_objective": objective,
                "mortality_gain": gains[0],
                "incidence_gain": gains[1],
                "new_alleles": new_alleles,
                "new_target_alleles": new_target_alleles,
                "new_monoallelic_target_alleles": new_mono_target_alleles,
                "new_ms_peptides": new_peptides,
                "new_proteoform": new_proteoform,
                "native_length_aa": len(segment.sequence(config.min_padding, config.min_padding)),
                "gain_per_aa": score[0],
            }
        )
    return [p[0] for p in chosen], history


def target_hla_support_table(segments, cancer, alleles=None, *, ligand_rows=None) -> pd.DataFrame:
    """Report per-proteoform HLA support, tier breakdown, and peptide redundancy (#201).

    Keeps measured monoallelic MS, sample-allele inference, and untyped panel
    inference distinct without multiplying carrier probabilities by expression
    prevalence.
    """
    columns = [
        "proteoform_key",
        "name",
        "rank",
        "n_segments",
        "native_aa",
        "distinct_ms_peptides",
        "supported_panel_alleles",
        "monoallelic_alleles",
        "sample_allele_alleles",
        "unrestricted_alleles",
        "mean_peptides_per_allele",
        "mortality_weighted_score",
        "measured_cancers",
        "missing_cancers",
    ]
    if not segments:
        return pd.DataFrame(columns=columns)

    by_key = {}
    for seg in segments:
        by_key.setdefault(seg.proteoform_key, []).append(seg)

    lig_by_key = {}
    if ligand_rows:
        for row in ligand_rows:
            lig_by_key.setdefault(row["proteoform_key"], []).append(row)

    cancer_by_key = {}
    if cancer is not None and not cancer.empty:
        for row in cancer.itertuples(index=False):
            cancer_by_key.setdefault(row.proteoform_key, []).append(row)

    rows = []
    for key, segs in sorted(by_key.items(), key=lambda kv: (kv[1][0].rank, kv[0])):
        peptides = {p for s in segs for p in s.peptides}
        supported_alleles = {a for s in segs for a in s.alleles}
        if alleles is not None:
            allowed = set(alleles)
            supported_alleles &= allowed

        ligs = lig_by_key.get(key, [])
        if alleles is not None:
            ligs = [r for r in ligs if r.get("allele") in allowed]
        tier_alleles = {
            tier: {r["allele"] for r in ligs if r.get("evidence_tier") == tier}
            for tier in ("monoallelic_ms", "sample_allele_ms", "unrestricted_ms")
        }
        peps_by_allele = {}
        for r in ligs:
            peps_by_allele.setdefault(r["allele"], set()).add(r["peptide"])
        if peps_by_allele:
            mean_redundancy = sum(len(v) for v in peps_by_allele.values()) / len(peps_by_allele)
        elif supported_alleles:
            mean_redundancy = len(peptides) / len(supported_alleles)
        else:
            mean_redundancy = 0.0

        crows = cancer_by_key.get(key, [])
        measured = sum(1 for r in crows if r.complete_measurement)
        missing = sum(1 for r in crows if not r.complete_measurement)
        rows.append(
            {
                "proteoform_key": key,
                "name": segs[0].name,
                "rank": segs[0].rank,
                "n_segments": len(segs),
                "native_aa": sum(s.specific_end - s.specific_start for s in segs),
                "distinct_ms_peptides": len(peptides),
                "supported_panel_alleles": len(supported_alleles),
                "monoallelic_alleles": len(tier_alleles["monoallelic_ms"]),
                "sample_allele_alleles": len(tier_alleles["sample_allele_ms"]),
                "unrestricted_alleles": len(tier_alleles["unrestricted_ms"]),
                "mean_peptides_per_allele": round(mean_redundancy, 4),
                "mortality_weighted_score": segs[0].score,
                "measured_cancers": measured,
                "missing_cancers": missing,
            }
        )
    return pd.DataFrame(rows, columns=columns)


def allocation_justification_table(
    selected_segments,
    all_segments,
    cancer,
    config,
    *,
    ligand_rows=None,
) -> pd.DataFrame:
    """Compute final-set leave-one-out removal and same-budget replacement effects (#202)."""
    columns = [
        "proteoform_key",
        "name",
        "retained_segments",
        "retained_aa",
        "delta_mortality_lower_pp",
        "delta_incidence_lower_pp",
        "delta_distinct_ms_peptides",
        "delta_global_alleles",
        "delta_target_alleles",
        "delta_monoallelic_target_alleles",
        "replacement_segments",
        "replacement_names",
        "replacement_aa",
        "replacement_delta_mortality_pp",
        "replacement_delta_ms_peptides",
        "replacement_delta_target_alleles",
    ]
    if not selected_segments:
        return pd.DataFrame(columns=columns)

    prevalence, weights = _expression_tables(cancer)
    tier_by_seg = _segment_tier_alleles(all_segments, ligand_rows)
    base = _evaluate_segment_set(
        selected_segments, prevalence, weights, tier_by_seg, config.min_padding
    )
    selected_ids = {s.segment_id for s in selected_segments}
    unselected = [s for s in all_segments if s.segment_id not in selected_ids]

    by_key = {}
    for seg in selected_segments:
        by_key.setdefault(seg.proteoform_key, []).append(seg)

    rows = []
    for key, segs in sorted(by_key.items(), key=lambda kv: (kv[1][0].rank, kv[0])):
        without = [s for s in selected_segments if s.proteoform_key != key]
        without_metrics = _evaluate_segment_set(
            without, prevalence, weights, tier_by_seg, config.min_padding
        )
        # Greedy same-budget replacement among unselected segments (excluding the removed protein).
        repl_pool = [s for s in unselected if s.proteoform_key != key]
        repl_placements = [(s, config.min_padding, config.min_padding, "") for s in without]
        repl_chosen = []
        used_peps = {p for s in without for p in s.peptides}
        used_target_alleles = {(s.proteoform_key, a) for s in without for a in s.alleles}
        while repl_pool:
            candidates = []
            for cand in repl_pool:
                pl = (cand, config.min_padding, config.min_padding, "")
                if not construct_fits([*repl_placements, pl], config):
                    continue
                new_p = len(set(cand.peptides) - used_peps)
                new_ta = len({(cand.proteoform_key, a) for a in cand.alleles} - used_target_alleles)
                if new_p == 0 and new_ta == 0:
                    continue
                cand_len = len(cand.sequence(config.min_padding, config.min_padding))
                candidates.append(((new_p / cand_len, new_ta / cand_len, -cand.rank), cand, pl))
            if not candidates:
                break
            _, best_cand, best_pl = max(candidates, key=lambda x: (x[0], x[1].segment_id))
            repl_placements.append(best_pl)
            repl_chosen.append(best_cand)
            repl_pool.remove(best_cand)
            used_peps.update(best_cand.peptides)
            used_target_alleles.update((best_cand.proteoform_key, a) for a in best_cand.alleles)

        repl_metrics = _evaluate_segment_set(
            [*without, *repl_chosen], prevalence, weights, tier_by_seg, config.min_padding
        )
        rows.append(
            {
                "proteoform_key": key,
                "name": segs[0].name,
                "retained_segments": len(segs),
                "retained_aa": sum(
                    len(s.sequence(config.min_padding, config.min_padding)) for s in segs
                ),
                "delta_mortality_lower_pp": round(
                    base["mortality_lower_pp"] - without_metrics["mortality_lower_pp"], 6
                ),
                "delta_incidence_lower_pp": round(
                    base["incidence_lower_pp"] - without_metrics["incidence_lower_pp"], 6
                ),
                "delta_distinct_ms_peptides": (
                    base["distinct_ms_peptides"] - without_metrics["distinct_ms_peptides"]
                ),
                "delta_global_alleles": (
                    base["global_panel_alleles"] - without_metrics["global_panel_alleles"]
                ),
                "delta_target_alleles": (
                    base["target_allele_pairs"] - without_metrics["target_allele_pairs"]
                ),
                "delta_monoallelic_target_alleles": (
                    base["monoallelic_target_alleles"]
                    - without_metrics["monoallelic_target_alleles"]
                ),
                "replacement_segments": ";".join(s.segment_id for s in repl_chosen),
                "replacement_names": ";".join(dict.fromkeys(s.name for s in repl_chosen)),
                "replacement_aa": sum(
                    len(s.sequence(config.min_padding, config.min_padding)) for s in repl_chosen
                ),
                "replacement_delta_mortality_pp": round(
                    repl_metrics["mortality_lower_pp"] - base["mortality_lower_pp"], 6
                ),
                "replacement_delta_ms_peptides": (
                    repl_metrics["distinct_ms_peptides"] - base["distinct_ms_peptides"]
                ),
                "replacement_delta_target_alleles": (
                    repl_metrics["target_allele_pairs"] - base["target_allele_pairs"]
                ),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def compare_allocation_tradeoffs(
    segments,
    cancer,
    config,
    *,
    tolerances_pp=(0.0, 0.05, 0.1, 0.25),
    ligand_rows=None,
) -> pd.DataFrame:
    """Compare global-union, target-specific, and explicit tolerance tradeoffs under one size cap (#201/#202)."""
    prevalence, weights = _expression_tables(cancer)
    tier_by_seg = _segment_tier_alleles(segments, ligand_rows)

    runs = [
        (
            "global_union",
            0.0,
            replace(config, allocation_objective="global_union", min_mortality_gain_pp=0.0),
        ),
        (
            "target_specific",
            0.0,
            replace(config, allocation_objective="target_specific", min_mortality_gain_pp=0.0),
        ),
        *(
            (
                "evidence_tradeoff",
                float(tol),
                replace(
                    config,
                    allocation_objective="evidence_tradeoff",
                    min_mortality_gain_pp=float(tol),
                ),
            )
            for tol in tolerances_pp
        ),
    ]
    rows = []
    for policy, tol, cfg in runs:
        selected, _ = select_budget_segments(segments, cancer, cfg, ligand_rows=ligand_rows)
        metrics = _evaluate_segment_set(selected, prevalence, weights, tier_by_seg, cfg.min_padding)
        rows.append(
            {
                "policy": policy,
                "min_mortality_gain_pp": tol,
                **metrics,
            }
        )
    frame = pd.DataFrame(rows)
    if not frame.empty:
        max_mort = float(frame["mortality_lower_pp"].max())
        frame["mortality_loss_vs_max_pp"] = (max_mort - frame["mortality_lower_pp"]).round(6)
    return frame
