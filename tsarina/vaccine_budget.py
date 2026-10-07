"""Deterministic marginal-gain allocation under a total sequence budget."""

from __future__ import annotations

from .vaccine_construct import construct_fits


def select_budget_segments(segments, cancer, config):
    """Greedy whole-piece allocation, without a protein count cap.

    Primary gain is the increase in mortality-weighted max marginal p95
    prevalence per aa, then incidence gain per aa, newly supported panel
    alleles per aa and distinct observed peptides per aa. max(p) is the
    conservative union bound; no patient independence is assumed. Missing
    measurements contribute no *known* gain. Additional same-protein pieces
    can add evidence but cannot count expression twice.
    """
    prevalence, weights = {}, {}
    for row in cancer.itertuples(index=False):
        if row.complete_measurement:
            prevalence.setdefault(row.proteoform_key, {})[row.burden_category] = row.prevalence_p95
        weights[row.burden_category] = (
            row.world_mortality_pct / 100,
            row.world_incidence_pct / 100,
        )
    chosen, remaining, history = [], list(segments), []
    current, alleles, peptides = {}, set(), set()
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
            new_peptides = len(set(segment.peptides) - peptides)
            length = len(segment.sequence(config.min_padding, config.min_padding))
            if not any(gains) and not new_alleles and not new_peptides:
                continue
            score = (
                *[gain / length for gain in gains],
                new_alleles / length,
                new_peptides / length,
                -segment.rank,
            )
            options.append((score, segment.segment_id, placement, gains, new_alleles, new_peptides))
        if not options:
            break
        score, _, placement, gains, new_alleles, new_peptides = max(
            options, key=lambda x: (x[0], x[1])
        )
        segment = placement[0]
        chosen.append(placement)
        remaining.remove(segment)
        alleles.update(segment.alleles)
        peptides.update(segment.peptides)
        for category, value in prevalence.get(segment.proteoform_key, {}).items():
            current[category] = max(value, current.get(category, 0))
        history.append(
            {
                "step": len(chosen),
                "segment_id": segment.segment_id,
                "name": segment.name,
                "mortality_gain": gains[0],
                "incidence_gain": gains[1],
                "new_alleles": new_alleles,
                "new_ms_peptides": new_peptides,
                "gain_per_aa": score[0],
            }
        )
    return [p[0] for p in chosen], history
