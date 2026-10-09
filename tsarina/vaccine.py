"""Mortality-prioritized, proteoform-level CTA vaccine design."""

from __future__ import annotations

from dataclasses import asdict
from fnmatch import fnmatchcase
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version

import pandas as pd

from .vaccine_construct import (
    NoFeasibleConstruct,
    VaccineConfig,
    encode_construct,
    optimize_construct,
    reserve_target_segments,
)
from .vaccine_inputs import (
    VaccineInputs,
    load_vaccine_inputs,
    panel_ms_support,
    rank_proteoforms,
    resolve_background_cta_ids,
)
from .vaccine_sequences import shared_kmers, specific_intervals, supported_segments


def design_vaccine(
    config: VaccineConfig | None = None,
    inputs: VaccineInputs | None = None,
    *,
    cohorts=None,
    auto_fetch=False,
    output_dir=None,
    on_progress=None,
    affinity_fn=None,
    cleavage_fn=None,
):
    """Design one vaccine and optionally write the complete evidence bundle.

    Prediction callbacks are explicit dependency injection for offline research
    and tests; their provenance is labeled in the manifest. Production defaults
    use live Hitlist observations, MHCflurry and Pepsickle.
    """
    from .spanning import _resolve_alleles

    config = config or VaccineConfig()
    config.validate()
    if config.species == "canine":
        if inputs is None or inputs.canine_evidence is None:
            raise ValueError("Canine design requires a frozen canine input bundle")
        if cohorts is not None or auto_fetch:
            raise ValueError(
                "Canine design uses the frozen bundle, not human cohort/download defaults"
            )
        from .vaccine_canine import design_canine_vaccine

        return design_canine_vaccine(
            config,
            inputs,
            output_dir=output_dir,
            on_progress=on_progress,
            affinity_fn=affinity_fn,
            cleavage_fn=cleavage_fn,
        )
    if inputs is not None and (
        inputs.canine_evidence is not None or inputs.provenance.get("taxon") == 9615
    ):
        raise ValueError("Canine inputs require species='canine'")
    if output_dir is not None:
        from pathlib import Path

        out = Path(output_dir)
        if out.exists() and any(out.iterdir()):
            raise ValueError(f"Output directory already contains files: {out}")
    alleles = list(dict.fromkeys(_resolve_alleles(config.alleles, config.panel)))
    if not alleles:
        raise ValueError("Specify a nonempty HLA class-I allele panel")
    if any(not a.startswith(("HLA-A*", "HLA-B*", "HLA-C*")) for a in alleles):
        raise ValueError("Vaccine design requires exact HLA-A/B/C alleles")
    inputs = inputs or load_vaccine_inputs(
        config.definition, config.ensembl_release, cohorts, auto_fetch, on_progress
    )
    if inputs.provenance.get("definition", config.definition) != config.definition:
        raise ValueError("Input CTA definition disagrees with design configuration")
    ranking, cancer_cohorts, cancer_summary = rank_proteoforms(inputs)
    exceptions = {s.strip().upper() for s in config.allow_genes}
    patterns = [s.strip().upper() for s in config.exclude_gene_patterns]
    ranking["eligible"] = ranking.name.map(
        lambda name: (
            not any(
                gene.upper() not in exceptions
                and any(fnmatchcase(gene.upper(), pattern) for pattern in patterns)
                for gene in name.split("/")
            )
        )
    )
    ranking["selection_reason"] = ranking.eligible.map(
        {True: "not_screened", False: "excluded_gene_pattern"}
    )
    eligible = ranking[ranking.eligible].copy()
    if config.selection_mode in {"supported", "budget"}:
        eligible = eligible[eligible.mortality_weighted_score > 0]
        ranking.loc[
            ranking.eligible & (ranking.mortality_weighted_score <= 0), "selection_reason"
        ] = "non_positive_score"
    if eligible.empty or eligible.mortality_weighted_score.max() <= 0:
        raise ValueError("No positive mortality-weighted CTA prevalence is available")
    background_cta_ids, aliases = resolve_background_cta_ids(inputs.proteins, inputs.cta_gene_ids)
    from .vaccine_inputs import filter_ms_modality
    from .vaccine_tissues import tissue_map, verified_normal_ms

    full_peptides = {
        row.sequence[i : i + k]
        for row in eligible.itertuples(index=False)
        for k in range(8, 16)
        for i in range(len(row.sequence) - k + 1)
    }
    if inputs.ms_hits is None:
        from .indexing import load_ms_evidence

        full_hits = load_ms_evidence(
            peptides=full_peptides, drop_binding_assays=False, include_binding=True
        )
    else:
        full_hits = inputs.ms_hits.copy()
    raw_hits = full_hits
    full_hits, rejected_full_hits = filter_ms_modality(raw_hits)
    normal_exclusions, verified_hits, normal_summary, normal_audit, normal_provenance = (
        pd.DataFrame(),
        pd.DataFrame(),
        pd.DataFrame(),
        pd.DataFrame(),
        {},
    )
    if config.normal_ms_atlas_dir:
        normal_exclusions, verified_hits, normal_summary, normal_audit, normal_provenance = (
            verified_normal_ms(
                config.normal_ms_atlas_dir, full_peptides, config.normal_ms_min_donors
            )
        )
    normal_kmers = (
        {p[i : i + 8] for p in normal_exclusions.get("peptide", []) for i in range(len(p) - 7)}
        if config.normal_ms_policy == "exclude"
        else set()
    )
    batch = max(10, 2 * config.top_k) if config.selection_mode == "supported" else config.top_k
    inspected = len(eligible) if config.selection_mode == "budget" else min(batch, len(eligible))
    allocation_history = []
    allocation_ids = None
    candidate_funnel = []
    before_normal_rows = []
    background_intervals = {}
    while True:
        candidates = eligible.head(inspected)
        if on_progress:
            on_progress(
                f"Screening {len(candidates)} eligible proteoforms; subtracting non-CTA 8-mers"
            )
        forbidden = shared_kmers(
            inputs.proteins, background_cta_ids, candidates.sequence, config.shared_k
        )
        background_intervals = {
            r.proteoform_key: specific_intervals(r.sequence, forbidden, config.shared_k)
            for r in candidates.itertuples(index=False)
        }
        before_normal_rows = [
            {
                "proteoform_key": r.proteoform_key,
                "name": r.name,
                "start": a,
                "end": b,
                "length_aa": b - a,
            }
            for r in candidates.itertuples(index=False)
            for a, b in background_intervals[r.proteoform_key]
        ]
        forbidden.update(normal_kmers)
        intervals = {
            r.proteoform_key: specific_intervals(r.sequence, forbidden, config.shared_k)
            for r in candidates.itertuples(index=False)
        }
        interval_rows, peptides = [], set()
        for row in candidates.itertuples(index=False):
            for start, end in intervals[row.proteoform_key]:
                interval_rows.append(
                    {
                        "proteoform_key": row.proteoform_key,
                        "name": row.name,
                        "start": start,
                        "end": end,
                        "length_aa": end - start,
                        "sequence": row.sequence[start:end],
                    }
                )
                for k in config.lengths:
                    peptides.update(row.sequence[a : a + k] for a in range(start, end - k + 1))
        if on_progress:
            on_progress(f"Reading panel MS evidence for {len(peptides)} CTA-specific peptides")
        support = panel_ms_support(
            peptides,
            alleles,
            raw_hits,
            config.predictor,
            mode=config.ms_support_mode,
            affinity_nm=config.ms_affinity_nm,
            allow_untyped=config.allow_untyped_ms,
        )
        segments, ligand_rows = supported_segments(candidates, intervals, support)
        supported_keys = {s.proteoform_key for s in segments}
        candidate_funnel = [
            {
                "proteoform_key": row.proteoform_key,
                "name": row.name,
                "raw_aa": len(row.sequence),
                "specific_aa": sum(b - a for a, b in background_intervals[row.proteoform_key]),
                "specific_pieces": len(background_intervals[row.proteoform_key]),
                "normal_ms_filtered_aa": sum(b - a for a, b in intervals[row.proteoform_key]),
                "normal_ms_filtered_pieces": len(intervals[row.proteoform_key]),
                "ms_supported_pieces": sum(
                    s.proteoform_key == row.proteoform_key for s in segments
                ),
                "ms_peptides": len(
                    {r["peptide"] for r in ligand_rows if r["proteoform_key"] == row.proteoform_key}
                ),
            }
            for row in candidates.itertuples(index=False)
        ]
        if config.selection_mode == "budget":
            from .vaccine_budget import select_budget_segments

            allocated, allocation_history = select_budget_segments(segments, cancer_summary, config)
            allocation_ids = {s.segment_id for s in allocated}
            chosen = {s.proteoform_key for s in allocated}
            length_rejected = []
            break
        if config.selection_mode == "ranked":
            chosen = set(candidates.proteoform_key)
            length_rejected = []
            break
        reserved, length_rejected = reserve_target_segments(segments, config)
        chosen = {placement[0].proteoform_key for placement in reserved}
        if len(chosen) >= config.top_k or inspected == len(eligible):
            break
        inspected = min(inspected + batch, len(eligible))
    ranking["selected"] = ranking.proteoform_key.isin(chosen)
    for row in candidates.itertuples(index=False):
        reason = "selected"
        if row.proteoform_key not in chosen:
            reason = (
                "no_cta_specific_sequence"
                if not background_intervals[row.proteoform_key]
                else "no_sequence_after_normal_ms"
                if not intervals[row.proteoform_key]
                else "no_qualified_panel_ms_ligands"
                if row.proteoform_key not in supported_keys
                else "minimum_segment_exceeds_length"
                if row.proteoform_key in length_rejected
                else "budget_not_selected"
                if config.selection_mode == "budget"
                else "lower_rank_supported"
            )
        ranking.loc[ranking.proteoform_key.eq(row.proteoform_key), "selection_reason"] = reason
    selected = ranking[ranking.selected].copy()
    segments = [s for s in segments if s.proteoform_key in chosen]
    ligand_rows = [r for r in ligand_rows if r["proteoform_key"] in chosen]
    target_count_error = (
        f"Only {len(selected)} supported proteoforms fit; requested {config.top_k}"
        if config.selection_mode == "supported" and len(selected) < config.top_k
        else None
    )
    if on_progress:
        on_progress(
            f"Selected {len(selected)} proteoforms after screening {len(candidates)} candidates"
        )
    design = None
    constraint_error = None
    if segments:
        try:
            design = optimize_construct(
                [s for s in segments if allocation_ids is None or s.segment_id in allocation_ids],
                alleles,
                config,
                affinity_fn,
                cleavage_fn,
                on_progress,
            )
        except NoFeasibleConstruct as error:
            constraint_error = str(error)
    if design:
        if allocation_ids is not None:
            design["excluded_segments"].extend(
                {"segment_id": s.segment_id, "reason": "budget_allocation"}
                for s in segments
                if s.segment_id not in allocation_ids
            )
        cds, full = encode_construct(design, config)
        design.update(
            {
                "cds_nt": cds,
                "full_nt": full,
                "length_aa": len(design["protein"]),
                "length_nt": len(full),
            }
        )
    # Map every native ligand to its actual final position, retaining exclusions.
    placements = (
        {layer["segment_id"]: layer for layer in design["layers"] if layer["kind"] == "cta_segment"}
        if design
        else {}
    )
    for ligand in ligand_rows:
        layer = placements.get(ligand["segment_id"])
        ligand["assembled"] = layer is not None
        ligand["construct_start"] = (
            layer["start_aa"] + ligand["start"] - layer["native_start"] if layer else None
        )
        ligand["construct_end"] = (
            layer["start_aa"] + ligand["end"] - layer["native_start"] if layer else None
        )
    funnel = []
    for row in selected.itertuples(index=False):
        pieces = intervals[row.proteoform_key]
        supported = [s for s in segments if s.proteoform_key == row.proteoform_key]
        assembled = (
            [
                layer
                for layer in design["layers"]
                if layer.get("proteoform_key") == row.proteoform_key
            ]
            if design
            else []
        )
        specific_aa = sum(b - a for a, b in background_intervals[row.proteoform_key])
        normal_filtered_aa = sum(b - a for a, b in pieces)
        ms_aa = sum(s.specific_end - s.specific_start for s in supported)
        padded_aa = sum(len(s.sequence(config.max_padding, config.max_padding)) for s in supported)
        allocated = [s for s in supported if s.segment_id in placements]
        allocated_aa = sum(
            len(s.sequence(config.max_padding, config.max_padding)) for s in allocated
        )
        assembled_aa = sum(len(layer["sequence"]) for layer in assembled)
        status = "retained"
        if not row.sequence:
            status = "missing_protein_sequence"
        elif not specific_aa:
            status = "no_cta_specific_sequence"
        elif not normal_filtered_aa:
            status = "no_sequence_after_normal_ms"
        elif not supported:
            status = "no_qualified_panel_ms_ligands"
        elif not assembled:
            status = "construct_length_limit"
        funnel.append(
            {
                "rank": row.rank,
                "proteoform_key": row.proteoform_key,
                "name": row.name,
                "status": status,
                "raw_aa": row.length_aa,
                "raw_pieces": int(bool(row.length_aa)),
                "specific_aa": specific_aa,
                "specific_pieces": len(background_intervals[row.proteoform_key]),
                "normal_ms_filtered_aa": normal_filtered_aa,
                "normal_ms_filtered_pieces": len(pieces),
                "ms_supported_aa": ms_aa,
                "ms_supported_pieces": len(supported),
                "padded_aa": padded_aa,
                "padded_pieces": len(supported),
                "allocated_aa": allocated_aa,
                "allocated_pieces": len(allocated),
                "budget_excluded_aa": padded_aa - allocated_aa,
                "terminal_trimmed_aa": allocated_aa - assembled_aa,
                "assembled_aa": assembled_aa,
                "assembled_pieces": len(assembled),
                "retained_fraction": assembled_aa / row.length_aa if row.length_aa else 0,
                "panel_alleles": ";".join(sorted({a for s in supported for a in s.alleles})),
                "assembled_panel_alleles": ";".join(
                    sorted(
                        {
                            a
                            for s in supported
                            if any(layer.get("segment_id") == s.segment_id for layer in assembled)
                            for a in s.alleles
                        }
                    )
                ),
                "ms_ligand_count": len({p for s in supported for p in s.peptides}),
                "assembled_ms_ligand_count": len(
                    {
                        r["peptide"]
                        for r in ligand_rows
                        if r["proteoform_key"] == row.proteoform_key and r["assembled"]
                    }
                ),
                "assembled_pmhc_count": len(
                    {
                        (r["peptide"], r["allele"])
                        for r in ligand_rows
                        if r["proteoform_key"] == row.proteoform_key and r["assembled"]
                    }
                ),
            }
        )
    versions = {}
    for package in (
        "tsarina",
        "oncoref",
        "hitlist",
        "pyensembl",
        "mhcflurry",
        "mhctools",
        "pepsickle",
        "vaxrank",
    ):
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = None
    from hitlist.version import __version__ as hitlist_source_version

    from .version import __version__

    versions["tsarina_source"] = __version__
    versions["hitlist_source"] = hitlist_source_version
    source_tables = dict(inputs.source_tables)
    source_tables["cta_specific_before_normal_ms"] = pd.DataFrame(
        before_normal_rows, columns=["proteoform_key", "name", "start", "end", "length_aa"]
    )
    source_tables["candidate_funnel"] = pd.DataFrame(candidate_funnel)
    if config.selection_mode == "budget":
        source_tables["budget_allocation"] = pd.DataFrame(allocation_history)
    source_tables["background_cta_aliases"] = aliases
    source_tables["selection_screen"] = ranking[
        ["rank", "proteoform_key", "name", "eligible", "selected", "selection_reason"]
    ].copy()
    for key in ("queried_ms_observations", "rejected_ms_observations", "ms_assignments"):
        if key in support.attrs:
            source_tables[key] = support.attrs[key]
    source_tables["normal_ms_exclusions"] = normal_exclusions
    source_tables["normal_ms_summary"] = normal_summary
    source_tables["normal_ms_audit"] = normal_audit
    source_tables["protein_ms_map"] = tissue_map(
        candidates,
        pd.concat([full_hits, verified_hits], ignore_index=True),
        pd.DataFrame(
            interval_rows,
            columns=["proteoform_key", "name", "start", "end", "length_aa", "sequence"],
        ),
        design["layers"] if design else [],
        source_tables.get("ms_assignments"),
    )
    source_tables["full_protein_ms_observations"] = full_hits
    source_tables["rejected_full_protein_ms_observations"] = rejected_full_hits
    allele_counts = []
    for allele in alleles:
        rows = [r for r in ligand_rows if r["assembled"] and r["allele"] == allele]
        allele_counts.append(
            {
                "allele": allele,
                "retained_peptides": len({r["peptide"] for r in rows}),
                "retained_proteoforms": len({r["proteoform_key"] for r in rows}),
                **{
                    tier: len({r["peptide"] for r in rows if r["evidence_tier"] == tier})
                    for tier in ("monoallelic_ms", "sample_allele_ms", "unrestricted_ms")
                },
            }
        )
    source_tables["hla_support_counts"] = pd.DataFrame(allele_counts)
    provenance = {
        **inputs.provenance,
        "ms_input_kind": inputs.provenance.get(
            "ms_input_kind",
            "supplied_observations" if inputs.ms_hits is not None else "hitlist_observations_index",
        ),
        "prevalence_input_sha256": sha256(
            inputs.prevalence.to_csv(index=False).encode()
        ).hexdigest(),
        "cta_gene_ids_sha256": sha256("\n".join(sorted(inputs.cta_gene_ids)).encode()).hexdigest(),
        "background_cta_gene_ids_sha256": sha256(
            "\n".join(sorted(background_cta_ids)).encode()
        ).hexdigest(),
        "background_cta_alias_count": len(aliases),
        "background_cta_identity": "curated IDs plus exact-symbol HSCHR alternate-haplotype annotations of the same gene",
        "ms_input_sha256": support.attrs.get("ms_input_sha256"),
        "ms_modality_policy": "positive MS assay method, curated MS-only supplement, or declared supplied modality; explicit non-MS rejected",
        "inspected_proteoforms": len(candidates),
        "selected_proteoforms": len(selected),
        "background_translated_occurrences": sum(
            p.gene_id not in background_cta_ids for p in inputs.proteins
        ),
        "world_mortality_share_represented_pct": float(
            inputs.burden.set_index("burden_category")
            .loc[list(inputs.cohorts), "world_mortality_pct"]
            .sum()
        ),
        "tissue_map_scope": inputs.provenance.get("ms_source_provenance", {}).get(
            "query_scope",
            "Full-protein query of current index"
            if inputs.ms_hits is None
            else "Supplied observation set; completeness outside the supplied peptide scope is unknown",
        ),
        "budget_stop_reason": "no remaining whole segment fits with positive known expression or new allele/peptide gain"
        if config.selection_mode == "budget"
        else None,
        "verified_normal_ms": normal_provenance,
    }
    status = (
        "insufficient_supported_proteoforms"
        if target_count_error
        else "junction_review_required"
        if design and not design["clean_junctions"]
        else "designed"
        if design
        else "no_feasible_construct"
        if constraint_error
        else "no_ms_supported_construct"
    )
    funnel_columns = [
        "rank",
        "proteoform_key",
        "name",
        "status",
        "raw_aa",
        "raw_pieces",
        "specific_aa",
        "specific_pieces",
        "normal_ms_filtered_aa",
        "normal_ms_filtered_pieces",
        "ms_supported_aa",
        "ms_supported_pieces",
        "padded_aa",
        "padded_pieces",
        "allocated_aa",
        "allocated_pieces",
        "budget_excluded_aa",
        "terminal_trimmed_aa",
        "assembled_aa",
        "assembled_pieces",
        "retained_fraction",
        "panel_alleles",
        "assembled_panel_alleles",
        "ms_ligand_count",
        "assembled_ms_ligand_count",
        "assembled_pmhc_count",
    ]
    result = {
        "config": asdict(config),
        "alleles": alleles,
        "provenance": provenance,
        "versions": versions,
        "cohort_mapping": inputs.cohorts,
        "coordinate_system": "zero-based, half-open",
        "objective": "sum_category(world_mortality_pct / 100 * sample_weighted_prevalence_p95)",
        "prediction_callbacks": {
            "affinity": affinity_fn is not None,
            "cleavage": cleavage_fn is not None,
        },
        "ranking": ranking,
        "cancer_cohorts": cancer_cohorts,
        "cancer_summary": cancer_summary,
        "specific_intervals": pd.DataFrame(
            interval_rows,
            columns=["proteoform_key", "name", "start", "end", "length_aa", "sequence"],
        ),
        "ligands": pd.DataFrame(ligand_rows)
        if ligand_rows
        else pd.DataFrame(
            columns=[
                "proteoform_key",
                "segment_id",
                "peptide",
                "allele",
                "start",
                "end",
                "evidence_tier",
                "assembled",
                "construct_start",
                "construct_end",
            ]
        ),
        "funnel": pd.DataFrame(funnel, columns=funnel_columns),
        "source_tables": source_tables,
        "design": design,
        "status": status,
    }
    if output_dir is not None:
        from .vaccine_report import write_vaccine_report

        write_vaccine_report(result, output_dir)
    if target_count_error:
        raise NoFeasibleConstruct(target_count_error + "; see selection_screen.csv and report.md")
    if design is None:
        if constraint_error:
            raise NoFeasibleConstruct(constraint_error + "; see the evidence report")
        raise ValueError(
            "No selected CTA retains panel MS-supported sequence; see the evidence report"
        )
    if config.require_clean_junctions and not design["clean_junctions"]:
        raise ValueError(
            "Final construct has unresolved junction binders; see junctions.csv and report.md"
        )
    return result
