"""Canine evidence policy around Tsarina's shared native-segment constructor."""

from __future__ import annotations

import math
import sys
from dataclasses import asdict
from fnmatch import fnmatchcase
from pathlib import Path

import mhcgnomes
import pandas as pd

from .canine_inputs import canine_prevalence, canonical_hash, dla_alleles
from .vaccine_construct import (
    NoFeasibleConstruct,
    encode_construct,
    optimize_construct,
    reserve_target_segments,
)
from .vaccine_coverage import paired_target_coverage
from .vaccine_inputs import Protein
from .vaccine_sequences import (
    peptide_occurrences,
    shared_kmers,
    specific_intervals,
    supported_segments,
)
from .version import __version__


def _prediction_archive(bundle, affinity_fn):
    lookup = {
        (r["peptide"], dla_alleles([r["allele"]])[0]): r["affinity_nm"]
        for r in bundle.get("predictions", [])
    }
    used = {}

    def predict(peptides, alleles):
        wanted = {(p, a) for p in peptides for a in alleles}
        missing = wanted - used.keys()
        if missing:
            if affinity_fn is not None:
                frame = affinity_fn(sorted({p for p, _ in missing}), alleles)
                if frame.duplicated(["peptide", "allele"]).any():
                    raise ValueError("Duplicate DLA affinity predictions")
                supplied = {
                    (r.peptide, dla_alleles([r.allele])[0]): r.affinity_nm
                    for r in frame.itertuples(index=False)
                }
            else:
                supplied = lookup
            for key in sorted(missing):
                if key not in supplied:
                    raise ValueError(
                        f"Unassessed DLA affinity for {key[0]}/{key[1]}; supply a complete frozen archive or explicit affinity callback (mhctools#544, mhcflurry#490)"
                    )
                value = float(supplied[key])
                if not math.isfinite(value) or value <= 0:
                    raise ValueError(f"Invalid DLA affinity for {key}")
                used[key] = value
        return pd.DataFrame(
            [{"peptide": p, "allele": a, "affinity_nm": used[(p, a)]} for p, a in sorted(wanted)],
            columns=["peptide", "allele", "affinity_nm"],
        )

    return predict, used


def _ms_support(hits, alleles, peptides, config, predict):
    observed = set(hits.peptide) & set(peptides) if not hits.empty else set()
    hit_candidates = {}
    for hit in hits.to_dict("records"):
        candidates = set(dla_alleles(hit["sample_alleles"]))
        if hit["restriction_kind"] == "untyped":
            candidates = set(alleles) if config.allow_untyped_ms else set()
        hit_candidates[hit["observation_id"]] = candidates & set(alleles)
    scores = {}
    for allele in alleles:
        needed = sorted(
            {
                h["peptide"]
                for h in hits.to_dict("records")
                if h["peptide"] in observed and allele in hit_candidates[h["observation_id"]]
            }
        )
        scores.update(
            {
                (r.peptide, r.allele): r.affinity_nm
                for r in predict(needed, [allele]).itertuples(index=False)
            }
        )
    rows, decisions = [], []
    for hit in hits.to_dict("records"):
        peptide = hit["peptide"]
        reason = "outside_native_specific_intervals"
        matched = []
        if peptide in observed:
            candidates = hit_candidates[hit["observation_id"]]
            matched = [
                a
                for a in alleles
                if a in candidates and scores[(peptide, a)] < config.ms_affinity_nm
            ]
            reason = "qualified" if matched else "no_qualifying_sample_allele"
            for allele in matched:
                tier = {
                    "monoallelic": "monoallelic_ms",
                    "sample_genotype": "sample_allele_ms",
                    "untyped": "unrestricted_ms",
                }[hit["restriction_kind"]]
                rows.append(
                    {
                        "peptide": peptide,
                        "allele": allele,
                        "evidence_tier": tier,
                        "restriction_assignment": "observed_monoallelic"
                        if tier == "monoallelic_ms"
                        else "inferred_by_affinity",
                        "affinity_nm": float(scores[(peptide, allele)]),
                        "observation_id": hit["observation_id"],
                        "sample_id": hit["sample_id"],
                        "host_taxon": hit["host_taxon"],
                        "host_species": hit["host_species"],
                        "sample_kind": hit["sample_kind"],
                        "sample_health": hit["sample_health"],
                        "sample_tissue": hit["sample_tissue"],
                        "source_id": hit["source_id"],
                    }
                )
        decisions.append(
            {**hit, "support_decision": reason, "qualifying_alleles": ";".join(matched)}
        )
    return pd.DataFrame(
        rows,
        columns=[
            "peptide",
            "allele",
            "evidence_tier",
            "restriction_assignment",
            "affinity_nm",
            "observation_id",
            "sample_id",
            "host_taxon",
            "host_species",
            "sample_kind",
            "sample_health",
            "sample_tissue",
            "source_id",
        ],
    ), pd.DataFrame(decisions)


def design_canine_vaccine(
    config, inputs, *, output_dir=None, on_progress=None, affinity_fn=None, cleavage_fn=None
):
    bundle = inputs.canine_evidence
    if canonical_hash(bundle) != inputs.provenance.get("bundle_canonical_sha256"):
        raise ValueError("Canine bundle changed after import; reload the frozen input")
    if (
        inputs.provenance.get("taxon") != 9615
        or inputs.provenance.get("reference_key") != bundle["reference_key"]
    ):
        raise ValueError("Canine portable input reference/taxon provenance disagrees")
    expected = sorted(
        (r["gene_id"], r["protein_id"], r["transcript_id"], r["sequence"])
        for r in bundle["occurrences"]
    )
    actual = sorted((p.gene_id, p.protein_id, p.transcript_id, p.sequence) for p in inputs.proteins)
    if actual != expected:
        raise ValueError("Portable protein occurrences disagree with the frozen reference")
    if config.definition != bundle["restriction_policy"]["definition"]:
        raise ValueError("Canine restriction definition disagrees with the bundle policy")
    if bundle["genotype_population"]["cohort"] != config.canine_cohort:
        raise ValueError("Genotype sampling frame disagrees with the named RNA cohort")
    if config.selection_mode == "budget":
        raise ValueError(
            "Canine budget allocation requires a species-scoped marginal objective; use ranked/supported (tsarina#201)"
        )
    if config.codon_species == "h_sapiens":
        raise ValueError(
            "Specify a canine coding policy explicitly; generic reverse translation is available"
        )
    if output_dir is not None and Path(output_dir).exists() and any(Path(output_dir).iterdir()):
        raise ValueError("Output directory already contains files")
    ranking, donors, sample_audit = canine_prevalence(bundle, config.canine_cohort)
    alleles = dla_alleles(config.alleles if config.alleles is not None else bundle["panel"])
    capabilities = {dla_alleles([c["allele"]])[0]: c for c in bundle["capabilities"]}
    supported = []
    capability_rows = []
    for allele in alleles:
        cap = capabilities.get(
            allele, {"tier": "unsupported", "reason": "capability_missing", "lengths": []}
        )
        usable = cap["tier"] == "empirically_supported" or (
            config.allow_exploratory_dla and cap["tier"] == "sequence_extrapolated"
        )
        usable = usable and set(config.lengths) <= set(cap.get("lengths", []))
        if usable:
            supported.append(allele)
        capability_rows.append({**cap, "allele": allele, "selection_enabled": bool(usable)})
    if not config.allow_exploratory_dla and any(
        capabilities.get(a, {}).get("tier") == "sequence_extrapolated" for a in alleles
    ):
        # Preserve RNA output below; extrapolated alleles contribute no support.
        selection_limit = "sequence_extrapolated_DLA_requires_explicit_opt_in"
    else:
        selection_limit = "" if alleles else "no_DLA_panel_construct_unassessed"
    pairs = [
        {**p, "alleles": dla_alleles(p["alleles"])}
        for p in bundle["tumor_genotype_pairs"]
        if p["cohort"] == config.canine_cohort
    ]
    cohort_donors = set(donors.donor)
    if any(p["tumor_donor"] not in cohort_donors for p in pairs):
        raise ValueError(
            "Genotype pairs must refer to independent dogs admitted to the selected cohort"
        )
    missing_mass = bundle["genotype_population"]["missing_mass"]
    # All occurrence identities, not gene symbols or alternate-contig aliases,
    # govern sharing. A rejected/unknown group contributes every occurrence.
    admitted_keys = set(ranking.loc[ranking.restriction_status.eq("admitted"), "proteoform_key"])
    background = [
        Protein(r["occurrence_id"], r["gene_id"], r["protein_id"], r["sequence"])
        for r in bundle["occurrences"]
    ]
    admitted_occurrences = {
        r["occurrence_id"] for r in bundle["occurrences"] if r["sequence_id"] in admitted_keys
    }
    raw_hits = pd.DataFrame(bundle["ms_hits"])
    accepted, rejected = raw_hits.copy(), raw_hits.iloc[:0].copy()
    if not accepted.empty:
        # Consume the reviewed source decision with method/response precedence;
        # raw binding labels alone neither establish nor refute MS detection.
        positive = accepted.is_ms_observation.eq(True)
        unknown_result = accepted[~positive].assign(
            rejection_reason="positive_MS_source_decision_not_established"
        )
        rejected = pd.concat([rejected, unknown_result], ignore_index=True)
        accepted = accepted[positive].copy()
        elution = accepted.assay_context.eq("MHC_ligand_elution") & accepted.mhc_class.eq("I")
        rejected = pd.concat(
            [
                rejected,
                accepted[~elution].assign(rejection_reason="not_class_I_MHC_ligand_elution"),
            ],
            ignore_index=True,
        )
        accepted = accepted[elution].copy()
    allowed_tissues = set(bundle["restriction_policy"]["allowed_tissues"])
    normal_rows = [
        r
        for r in accepted.to_dict("records")
        if r["host_taxon"] == 9615
        and r["sample_kind"] == "primary_tissue"
        and r["sample_health"] == "healthy"
        and r.get("verified_normal") is True
        and r["sample_tissue"] not in allowed_tissues
    ]
    normal_kmers = {
        r["peptide"][i : i + 8] for r in normal_rows for i in range(len(r["peptide"]) - 7)
    }

    def background_check(sequences):
        return shared_kmers(background, admitted_occurrences, sequences) | normal_kmers

    exceptions = {s.upper() for s in config.allow_genes}
    patterns = [s.upper() for s in config.exclude_gene_patterns]
    ranking["eligible"] = ranking.apply(
        lambda r: (
            r.restriction_status == "admitted"
            and r.prevalence_lower > 0
            and not any(
                g.upper() not in exceptions and any(fnmatchcase(g.upper(), p) for p in patterns)
                for g in r["name"].split("/")
            )
        ),
        axis=1,
    )
    eligible = ranking[ranking.eligible].copy()
    candidates = eligible.head(config.top_k) if config.selection_mode == "ranked" else eligible
    forbidden = shared_kmers(background, admitted_occurrences, candidates.sequence)
    intervals_before_normal = {
        r.proteoform_key: specific_intervals(r.sequence, forbidden)
        for r in candidates.itertuples(index=False)
    }
    intervals = {
        r.proteoform_key: specific_intervals(r.sequence, forbidden | normal_kmers)
        for r in candidates.itertuples(index=False)
    }
    peptides = {
        r.sequence[i : i + k]
        for r in candidates.itertuples(index=False)
        for a, b in intervals[r.proteoform_key]
        for k in config.lengths
        for i in range(a, b - k + 1)
    }
    predict, used = _prediction_archive(bundle, affinity_fn)
    support, ms_decisions = _ms_support(accepted, supported, peptides, config, predict)
    segments, ligand_rows = supported_segments(
        candidates, intervals, support, score_column="prevalence_lower"
    )
    all_supported_segments = segments.copy()
    if config.selection_mode == "supported":
        reserved, _ = reserve_target_segments(segments, config)
        chosen = {p[0].proteoform_key for p in reserved}
        segments = [s for s in segments if s.proteoform_key in chosen]
    design, limit = None, selection_limit
    if segments:
        try:
            if on_progress:
                on_progress(
                    f"Assembling canine cohort {config.canine_cohort}: {len(segments)} MS-supported native pieces"
                )
            design = optimize_construct(
                segments,
                supported,
                config,
                predict,
                cleavage_fn,
                on_progress,
                background_fn=background_check,
            )
            design["coding_sequence"], design["nucleotide_sequence"] = encode_construct(
                design, config
            )
            design["length_aa"] = len(design["protein"])
            design["length_nt"] = len(design["nucleotide_sequence"])
            design["final_assessment"] = "unassessed_canine_presentation_and_processing"
            design["unassessed_alleles"] = sorted(set(alleles) - set(supported))
        except NoFeasibleConstruct as error:
            limit = str(error)
    placements = (
        {r["segment_id"]: r for r in design["layers"] if r["kind"] == "cta_segment"}
        if design
        else {}
    )
    for ligand in ligand_rows:
        layer = placements.get(ligand["segment_id"])
        retained = (
            layer is not None
            and layer["native_start"] <= ligand["start"]
            and ligand["end"] <= layer["native_end"]
        )
        ligand["assembled"] = retained
        ligand["construct_start"] = (
            layer["start_aa"] + ligand["start"] - layer["native_start"] if retained else None
        )
        ligand["construct_end"] = (
            layer["start_aa"] + ligand["end"] - layer["native_start"] if retained else None
        )
    retained_keys = {r["proteoform_key"] for r in placements.values()}
    ranking["selected"] = ranking.proteoform_key.isin(retained_keys)
    ranking["selection_reason"] = ranking.apply(
        lambda r: (
            "assembled"
            if r.selected
            else r.restriction_status
            if r.restriction_status != "admitted"
            else "not_eligible"
            if not r.eligible
            else "not_screened"
            if r.proteoform_key not in intervals
            else "no_specific_sequence"
            if not intervals[r.proteoform_key]
            else "not_retained_after_MS_and_length_gates"
        ),
        axis=1,
    )
    funnel = []
    for r in candidates.itertuples(index=False):
        native = [s for s in all_supported_segments if s.proteoform_key == r.proteoform_key]
        layers = [p for p in placements.values() if p["proteoform_key"] == r.proteoform_key]
        hits = [
            h for h in ligand_rows if h["proteoform_key"] == r.proteoform_key and h["assembled"]
        ]
        funnel.append(
            {
                "proteoform_key": r.proteoform_key,
                "name": r.name,
                "raw_aa": r.length_aa,
                "specific_aa": sum(b - a for a, b in intervals_before_normal[r.proteoform_key]),
                "specific_pieces": len(intervals_before_normal[r.proteoform_key]),
                "normal_ms_filtered_aa": sum(b - a for a, b in intervals[r.proteoform_key]),
                "normal_ms_filtered_pieces": len(intervals[r.proteoform_key]),
                "ms_supported_aa": sum(s.specific_end - s.specific_start for s in native),
                "ms_supported_pieces": len(native),
                "max_padded_aa": sum(
                    s.bounds(config.max_padding, config.max_padding)[1]
                    - s.bounds(config.max_padding, config.max_padding)[0]
                    for s in native
                ),
                "assembled_aa": sum(p["end_aa"] - p["start_aa"] for p in layers),
                "assembled_pieces": len(layers),
                "MS_observed_peptides": len({h["peptide"] for h in hits}),
                "MS_observations": len({h["observation_id"] for h in hits}),
                "peptide_DLA_pairs": len({(h["peptide"], h["allele"]) for h in hits}),
            }
        )
    ligands = pd.DataFrame(
        ligand_rows,
        columns=[
            "proteoform_key",
            "name",
            "segment_id",
            "specific_piece",
            "start",
            "end",
            *support.columns,
            "assembled",
            "construct_start",
            "construct_end",
        ],
    )
    expr = {
        sid: {
            r.donor: (r.expression_lower, r.expression_upper)
            for r in donors[donors.proteoform_key.eq(sid)].itertuples(index=False)
        }
        for sid in ranking.proteoform_key
    }
    curves, final_coverage = [], None
    selected_targets = {}
    aa = 0
    for step, layer in enumerate(placements.values(), 1):
        key = layer["proteoform_key"]
        selected_targets.setdefault(key, {"expression": expr[key], "alleles": set()})
        hit_rows = [
            r for r in ligand_rows if r["assembled"] and r["segment_id"] == layer["segment_id"]
        ]
        selected_targets[key]["alleles"].update(r["allele"] for r in hit_rows)
        aa = layer["end_aa"]
        coverage = paired_target_coverage(selected_targets, pairs, supported, missing_mass)
        prefix_hits = [r for r in ligand_rows if r["assembled"] and r["construct_end"] <= aa]
        curves.append(
            {
                "step": step,
                "name": layer["name"],
                "length_aa": aa,
                "n_proteins": len(selected_targets),
                "MS_observed_peptides": len({r["peptide"] for r in prefix_hits}),
                "observed_monoallelic_peptides": len(
                    {r["peptide"] for r in prefix_hits if r["evidence_tier"] == "monoallelic_ms"}
                ),
                "inferred_only_peptides": len(
                    {r["peptide"] for r in prefix_hits}
                    - {r["peptide"] for r in prefix_hits if r["evidence_tier"] == "monoallelic_ms"}
                ),
                **{k: v for k, v in coverage.items() if k not in {"pairs", "pairing_modes"}},
            }
        )
        final_coverage = coverage
    if final_coverage is None:
        final_coverage = paired_target_coverage({}, pairs, supported, missing_mass)
    status = "exploratory_unassessed" if design else "no_MS_supported_construct"
    if design and not design["background_safe"]:
        status = "rejected_background_overlap"
    elif design and not design["clean_junctions"]:
        status = "exploratory_junction_review_required"
    if config.selection_mode == "supported" and len(retained_keys) < config.top_k:
        status = "insufficient_supported_targets"
        limit = f"Only {len(retained_keys)} proteins assembled; requested {config.top_k}"
    protein_ms_map = []
    normal_ids = {r["observation_id"] for r in normal_rows}
    for protein in candidates.itertuples(index=False):
        for hit in accepted.to_dict("records"):
            matches = [
                (a, b, "exact_observed_peptide")
                for a, b in peptide_occurrences(protein.sequence, hit["peptide"])
            ]
            if hit["observation_id"] in normal_ids:
                kmers = {hit["peptide"][i : i + 8] for i in range(len(hit["peptide"]) - 7)}
                matches.extend(
                    (i, i + 8, "verified_healthy_MS_8mer_overlap")
                    for i in range(len(protein.sequence) - 7)
                    if protein.sequence[i : i + 8] in kmers
                )
            protein_ms_map.extend(
                {
                    "proteoform_key": protein.proteoform_key,
                    "name": protein.name,
                    "start": a,
                    "end": b,
                    "mapping": kind,
                    **hit,
                }
                for a, b, kind in matches
            )
    result = {
        "species": "canine",
        "taxon": 9615,
        "status": status,
        "limit": limit,
        "config": asdict(config),
        "alleles": alleles,
        "provenance": inputs.provenance,
        "reference": bundle["reference"],
        "sources": bundle["sources"],
        "versions": {
            "tsarina_source": __version__,
            "mhcgnomes_source": getattr(mhcgnomes, "__version__", "unknown"),
            "pandas_source": pd.__version__,
            "python": sys.version.split()[0],
        },
        "restriction_policy": bundle["restriction_policy"],
        "expression_policy": bundle["expression_policy"],
        "cohort": {"name": config.canine_cohort, **bundle["cohorts"][config.canine_cohort]},
        "genotype_population": bundle["genotype_population"],
        "objective": "independent-dog absolute-TPM prevalence lower bound in one named cohort",
        "prediction_provider": "injected_affinity_callback"
        if affinity_fn
        else "frozen_affinity_archive",
        "prediction_callbacks": {
            "affinity": affinity_fn is not None,
            "cleavage": cleavage_fn is not None,
        },
        "coordinate_system": "zero-based, half-open",
        "ranking": ranking,
        "ligands": ligands,
        "funnel": pd.DataFrame(funnel),
        "design": design,
        "coverage": final_coverage,
        "source_tables": {
            "occurrences": pd.DataFrame(bundle["occurrences"]),
            "samples": pd.DataFrame(bundle["samples"]),
            "rna_bounds": pd.DataFrame(bundle["rna_bounds"]),
            "donor_prevalence": donors,
            "cohort_sample_audit": sample_audit,
            "ms_raw": raw_hits,
            "ms_modality_rejected": rejected,
            "ms_support_decisions": ms_decisions,
            "normal_ms_exclusions": pd.DataFrame(normal_rows),
            "dla_capabilities": pd.DataFrame(capability_rows),
            "capabilities_raw": pd.DataFrame(bundle["capabilities"]),
            "tumor_genotype_pairs": pd.DataFrame(pairs),
            "paired_coverage": pd.DataFrame(final_coverage["pairs"]),
            "cumulative_coverage": pd.DataFrame(curves),
            "affinity_predictions_used": pd.DataFrame(
                [
                    {"peptide": p, "allele": a, "affinity_nm": v}
                    for (p, a), v in sorted(used.items())
                ]
            ),
            "native_specific_intervals": pd.DataFrame(
                [
                    {
                        "proteoform_key": r.proteoform_key,
                        "name": r.name,
                        "start": a,
                        "end": b,
                        "sequence": r.sequence[a:b],
                    }
                    for r in candidates.itertuples(index=False)
                    for a, b in intervals[r.proteoform_key]
                ]
            ),
            "protein_ms_map": pd.DataFrame(protein_ms_map),
        },
    }
    from .canine_report import canine_counts

    result["counts"] = canine_counts(result, bundle)
    if output_dir is not None:
        from .canine_report import write_canine_report

        write_canine_report(result, bundle, output_dir)
    return result
