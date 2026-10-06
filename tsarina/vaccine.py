"""Mortality-prioritized, proteoform-level CTA vaccine design."""

from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version

import pandas as pd

from .vaccine_construct import (
    NoFeasibleConstruct,
    VaccineConfig,
    encode_construct,
    optimize_construct,
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
    ranking["selected"] = ranking["rank"] <= config.top_k
    selected = ranking[ranking.selected].copy()
    if selected.empty or selected.mortality_weighted_score.max() <= 0:
        raise ValueError("No positive mortality-weighted CTA prevalence is available")
    if on_progress:
        on_progress(f"Selected {len(selected)} proteoforms; subtracting non-CTA 8-mer intervals")
    background_cta_ids, aliases = resolve_background_cta_ids(inputs.proteins, inputs.cta_gene_ids)
    forbidden = shared_kmers(
        inputs.proteins, background_cta_ids, selected.sequence, config.shared_k
    )
    intervals = {
        r.proteoform_key: specific_intervals(r.sequence, forbidden, config.shared_k)
        for r in selected.itertuples(index=False)
    }
    interval_rows, peptides = [], set()
    for row in selected.itertuples(index=False):
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
    support = panel_ms_support(peptides, alleles, inputs.ms_hits, config.predictor)
    segments, ligand_rows = supported_segments(selected, intervals, support)
    design = None
    constraint_error = None
    if segments:
        try:
            design = optimize_construct(
                segments, alleles, config, affinity_fn, cleavage_fn, on_progress
            )
        except NoFeasibleConstruct as error:
            constraint_error = str(error)
    if design:
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
        specific_aa = sum(b - a for a, b in pieces)
        ms_aa = sum(s.specific_end - s.specific_start for s in supported)
        padded_aa = sum(len(s.sequence(config.max_padding, config.max_padding)) for s in supported)
        assembled_aa = sum(len(layer["sequence"]) for layer in assembled)
        status = "retained"
        if not row.sequence:
            status = "missing_protein_sequence"
        elif not specific_aa:
            status = "no_cta_specific_sequence"
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
                "specific_pieces": len(pieces),
                "ms_supported_aa": ms_aa,
                "ms_supported_pieces": len(supported),
                "padded_aa": padded_aa,
                "padded_pieces": len(supported),
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
    source_tables["background_cta_aliases"] = aliases
    if "queried_ms_observations" in support.attrs:
        source_tables["queried_ms_observations"] = support.attrs["queried_ms_observations"]
    provenance = {
        **inputs.provenance,
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
        "background_translated_occurrences": sum(
            p.gene_id not in background_cta_ids for p in inputs.proteins
        ),
        "world_mortality_share_represented_pct": float(
            inputs.burden.set_index("burden_category")
            .loc[list(inputs.cohorts), "world_mortality_pct"]
            .sum()
        ),
    }
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
        "funnel": pd.DataFrame(funnel),
        "source_tables": source_tables,
        "design": design,
        "status": (
            "junction_review_required" if design and not design["clean_junctions"] else "designed"
        )
        if design
        else ("no_feasible_construct" if constraint_error else "no_ms_supported_construct"),
    }
    if output_dir is not None:
        from .vaccine_report import write_vaccine_report

        write_vaccine_report(result, output_dir)
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
