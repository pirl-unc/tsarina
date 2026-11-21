"""Tests for vaccine search stop/continuation (#191), clean-junction backfill (#192),
target-specific HLA allocation (#201), and evidence tradeoff justification (#202)."""

from __future__ import annotations

import json
from dataclasses import replace

import pandas as pd
import pytest

from tsarina.vaccine import design_vaccine
from tsarina.vaccine_budget import (
    allocation_justification_table,
    compare_allocation_tradeoffs,
    select_budget_segments,
    target_hla_support_table,
)
from tsarina.vaccine_construct import VaccineConfig, optimize_construct
from tsarina.vaccine_inputs import Protein, VaccineInputs
from tsarina.vaccine_sequences import Segment

ALLELE = "HLA-A*02:01"


def _seg(
    key, seq, *, proteoform_key=None, name=None, rank=1, score=0.1, alleles=(ALLELE,), peptides=None
):
    pkey = proteoform_key or key
    return Segment(
        segment_id=key,
        proteoform_key=pkey,
        name=name or pkey,
        rank=rank,
        score=score,
        protein_sequence=seq,
        specific_start=0,
        specific_end=len(seq),
        ligand_start=0,
        ligand_end=len(seq),
        alleles=tuple(alleles),
        peptides=tuple(peptides) if peptides is not None else (seq[:8],),
    )


def _cleavage(sequences):
    return {s: [0.5] * len(s) for s in sequences}


@pytest.fixture
def codons(monkeypatch):
    from Bio.Data.CodonTable import standard_dna_table

    table = {aa: c for c, aa in standard_dna_table.forward_table.items()}
    monkeypatch.setattr(
        "tsarina.vaccine_elements.codon_optimize",
        lambda protein, species: "".join(table[aa] for aa in protein),
    )


def test_search_stop_reason_and_verified_checkpoint_continuation_191(tmp_path):
    """tsarina#191: report budget_exhausted_improving when improving at cap,
    produce bit-for-bit identical results when resumed from checkpoint while
    reusing cached predictions, and reject incompatible checkpoints."""
    segments = [
        _seg("s1", "ACDEFGHIK", rank=1),
        _seg("s2", "LMNPQRSTV", rank=2),
        _seg("s3", "WYACDEFGH", rank=3),
        _seg("s4", "IKLMNPQRS", rank=4),
    ]
    eval_calls = []

    def affinity_fn(peptides, alleles):
        eval_calls.extend(peptides)
        # Penalize every junction window containing 'K' or 'W' except when separated by AAY
        rows = []
        for p in peptides:
            for a in alleles:
                aff = 50.0 if ("K" in p or "W" in p) and "AAY" not in p else 5000.0
                rows.append({"peptide": p, "allele": a, "affinity_nm": aff})
        return pd.DataFrame(rows)

    cfg_full = VaccineConfig(
        alleles=(ALLELE,),
        min_padding=0,
        max_padding=0,
        linkers=("", "AAY"),
        beam_width=2,
        optimization_rounds=3,
    )
    uninterrupted = optimize_construct(
        segments, [ALLELE], cfg_full, affinity_fn=affinity_fn, cleavage_fn=_cleavage
    )

    ckpt_path = tmp_path / "search_ckpt.json"
    eval_calls.clear()
    cfg_part1 = replace(cfg_full, optimization_rounds=1)
    part1 = optimize_construct(
        segments,
        [ALLELE],
        cfg_part1,
        affinity_fn=affinity_fn,
        cleavage_fn=_cleavage,
        checkpoint_out=ckpt_path,
    )
    assert ckpt_path.exists()
    assert part1["search_stop_reason"] in {
        "clean",
        "converged",
        "budget_exhausted_improving",
        "budget_exhausted_plateau",
    }
    assert part1["search_progress"]["rounds_completed"] == 1
    assert part1["search_progress"]["rounds_budget"] == 1

    calls_before_resume = len(eval_calls)
    resumed = optimize_construct(
        segments,
        [ALLELE],
        cfg_full,
        affinity_fn=affinity_fn,
        cleavage_fn=_cleavage,
        checkpoint_in=ckpt_path,
    )
    assert resumed["protein"] == uninterrupted["protein"]
    assert resumed["final_objective"] == uninterrupted["final_objective"]
    assert resumed["search_history"] == uninterrupted["search_history"]
    assert resumed["search_progress"]["resumed_from_round"] == 1
    assert (
        resumed["search_progress"]["rounds_completed"]
        == uninterrupted["search_progress"]["rounds_completed"]
    )
    assert len(eval_calls) >= calls_before_resume  # Only new peptides predicted on resume.

    # Incompatible config/allele change must be rejected loudly.
    bad_cfg = replace(cfg_full, linkers=("", "GGG"))
    with pytest.raises(ValueError, match="Incompatible vaccine search checkpoint"):
        optimize_construct(
            segments,
            [ALLELE],
            bad_cfg,
            affinity_fn=affinity_fn,
            cleavage_fn=_cleavage,
            checkpoint_in=ckpt_path,
        )


def test_search_reports_budget_exhausted_improving_on_progressive_objective_191():
    """tsarina#191: an objective that still improves in the final round reports
    budget_exhausted_improving rather than converged or infeasible."""
    segments = [
        _seg("s1", "ACDEFGHI", rank=1),
        _seg("s2", "LMNPQRST", rank=2),
        _seg("s3", "VWYFGHIK", rank=3),
    ]

    def progressive_affinity(peptides, alleles):
        # Every direct join without AAY binds; inserting AAY at each join fixes one join per round.
        rows = []
        for p in peptides:
            for a in alleles:
                aff = 50.0 if "AAY" not in p else 5000.0
                rows.append({"peptide": p, "allele": a, "affinity_nm": aff})
        return pd.DataFrame(rows)

    cfg = VaccineConfig(
        alleles=(ALLELE,),
        min_padding=0,
        max_padding=0,
        linkers=("", "AAY"),
        beam_width=1,
        optimization_rounds=1,
    )
    res = optimize_construct(
        segments, [ALLELE], cfg, affinity_fn=progressive_affinity, cleavage_fn=_cleavage
    )
    assert res["search_progress"]["improved_in_final_round"] is True
    assert res["search_stop_reason"] == "budget_exhausted_improving"


def test_supported_clean_junction_backfill_replaces_unclean_target_192(
    monkeypatch, codons, tmp_path
):
    """tsarina#192 exact reproduction: CTA1+CTA2 fails clean-junction gate due to
    'K' in CTA2, so supported mode with require_clean_junctions=True backfills
    to clean CTA1+CTA3 and records panel_search_audit."""
    sequences = ["ACDEFGHI", "KLMNPQRS", "TVWYACDE"]
    inputs = VaccineInputs(
        proteins=[Protein(f"g{i}", f"CTA{i}", f"p{i}", s) for i, s in enumerate(sequences, 1)],
        cta_gene_ids={"g1", "g2", "g3"},
        gene_keys={f"g{i}": f"g{i}" for i in range(1, 4)},
        prevalence=pd.DataFrame(
            [
                {
                    "proteoform_key": f"g{i}",
                    "cancer_code": "C",
                    "prevalence_p95": 1 - i / 10,
                    "n_samples": 10,
                }
                for i in range(1, 4)
            ]
        ),
        burden=pd.DataFrame(
            [
                {
                    "burden_category": "cancer",
                    "world_mortality_pct": 20.0,
                    "world_incidence_pct": 12.0,
                    "world_mortality_count": float("nan"),
                }
            ]
        ),
        cohorts={"cancer": ["C"]},
        ms_hits=pd.DataFrame(
            [
                {
                    "peptide": s,
                    "mhc_restriction": ALLELE,
                    "is_monoallelic": True,
                    "assay_method": "cellular MHC/mass spectrometry",
                    "qualitative_measurement": "Positive",
                    "pmid": "synthetic",
                }
                for s in sequences
            ]
        ),
        provenance={"synthetic": True},
    )

    def ms_scores(peptides, alleles, predictor):
        return pd.DataFrame(
            [
                {
                    "peptide": p,
                    "allele": a,
                    "affinity_nm": 50.0,
                    "presentation_score": 0.9,
                    "presentation_percentile": 0.2,
                }
                for p in peptides
                for a in alleles
            ]
        )

    def junction_scores(peptides, alleles):
        return pd.DataFrame(
            [
                {"peptide": p, "allele": a, "affinity_nm": 10.0 if "K" in p else 5000.0}
                for p in peptides
                for a in alleles
            ]
        )

    monkeypatch.setattr("tsarina.scoring.score_presentation", ms_scores)
    cfg = VaccineConfig(
        top_k=2,
        selection_mode="supported",
        alleles=(ALLELE,),
        min_padding=0,
        max_padding=0,
        optimization_rounds=2,
        beam_width=4,
        max_length_aa=30,
        ms_support_mode="sample_affinity",
        require_clean_junctions=True,
    )

    result = design_vaccine(
        cfg, inputs=inputs, output_dir=tmp_path, affinity_fn=junction_scores, cleavage_fn=_cleavage
    )
    assert result["status"] == "designed"
    assert result["design"]["clean_junctions"] is True
    assert result["ranking"].query("selected")["name"].tolist() == ["CTA1", "CTA3"]
    assert (
        result["ranking"].set_index("name").loc["CTA2", "selection_reason"]
        == "junction_backfill_excluded"
    )
    assert result["provenance"]["panel_search_status"] == "clean_panel_found"
    audit = result["source_tables"]["panel_search_audit"]
    assert len(audit) == 2
    assert audit.iloc[0]["status"] == "junction_review_required"
    assert audit.iloc[1]["status"] == "clean"
    assert (tmp_path / "panel_search_audit.csv").exists()

    # Genuinely impossible clean-panel case: every peptide triggers a binder.
    def all_bind(peptides, alleles):
        return pd.DataFrame(
            [{"peptide": p, "allele": a, "affinity_nm": 10.0} for p in peptides for a in alleles]
        )

    out_fail = tmp_path / "impossible"
    with pytest.raises(ValueError, match="unresolved junction binders"):
        design_vaccine(
            cfg, inputs=inputs, output_dir=out_fail, affinity_fn=all_bind, cleavage_fn=_cleavage
        )
    manifest = json.loads((out_fail / "manifest.json").read_text())
    assert manifest["provenance"]["panel_search_status"] == "all_combinations_exhausted"
    assert (out_fail / "panel_search_audit.csv").exists()


def test_target_specific_hla_objective_and_support_table_201():
    """tsarina#201: target_specific objective credits new (proteoform, allele) pairs
    even when another CTA already covers those alleles in the global construct union."""
    other_cta = _seg(
        "other_1",
        "ACDEFGHIK",
        proteoform_key="other",
        name="OTHER_CTA",
        rank=1,
        score=0.20,
        alleles=("HLA-A*01:01", "HLA-A*68:01"),
        peptides=("ACDEFGHI",),
    )
    prame_1 = _seg(
        "prame_1",
        "LMNPQRSTV",
        proteoform_key="prame",
        name="PRAME",
        rank=2,
        score=0.18,
        alleles=("HLA-A*02:01",),
        peptides=("LMNPQRST",),
    )
    # prame_2 adds HLA-A*01:01 and HLA-A*68:01 specifically to PRAME (0 new global alleles!).
    prame_2 = _seg(
        "prame_2",
        "WYACDEFGH",
        proteoform_key="prame",
        name="PRAME",
        rank=2,
        score=0.18,
        alleles=("HLA-A*01:01", "HLA-A*68:01"),
        peptides=("WYACDEFG",),
    )
    # filler has 1 new peptide on HLA-A*02:01 (already on PRAME) and same length.
    filler = _seg(
        "filler_1",
        "IKLMNPQRS",
        proteoform_key="other",
        name="OTHER_CTA",
        rank=1,
        score=0.20,
        alleles=("HLA-A*01:01",),
        peptides=("IKLMNPQR", "KLMNPQRS"),
    )
    cancer = pd.DataFrame(
        [
            {
                "proteoform_key": "other",
                "burden_category": "lung",
                "prevalence_p95": 0.8,
                "complete_measurement": True,
                "world_mortality_pct": 20.0,
                "world_incidence_pct": 10.0,
            },
            {
                "proteoform_key": "prame",
                "burden_category": "breast",
                "prevalence_p95": 0.5,
                "complete_measurement": True,
                "world_mortality_pct": 15.0,
                "world_incidence_pct": 10.0,
            },
        ]
    )
    ligand_rows = [
        {
            "proteoform_key": "other",
            "segment_id": "other_1",
            "peptide": "ACDEFGHI",
            "allele": "HLA-A*01:01",
            "evidence_tier": "monoallelic_ms",
        },
        {
            "proteoform_key": "other",
            "segment_id": "other_1",
            "peptide": "ACDEFGHI",
            "allele": "HLA-A*68:01",
            "evidence_tier": "sample_allele_ms",
        },
        {
            "proteoform_key": "prame",
            "segment_id": "prame_1",
            "peptide": "LMNPQRST",
            "allele": "HLA-A*02:01",
            "evidence_tier": "monoallelic_ms",
        },
        {
            "proteoform_key": "prame",
            "segment_id": "prame_2",
            "peptide": "WYACDEFG",
            "allele": "HLA-A*01:01",
            "evidence_tier": "monoallelic_ms",
        },
        {
            "proteoform_key": "prame",
            "segment_id": "prame_2",
            "peptide": "WYACDEFG",
            "allele": "HLA-A*68:01",
            "evidence_tier": "monoallelic_ms",
        },
        {
            "proteoform_key": "other",
            "segment_id": "filler_1",
            "peptide": "IKLMNPQR",
            "allele": "HLA-A*01:01",
            "evidence_tier": "unrestricted_ms",
        },
    ]

    # Budget fits 3 segments (1 + 27 = 28 aa).
    cfg_global = VaccineConfig(
        selection_mode="budget",
        max_length_aa=28,
        min_padding=0,
        allocation_objective="global_union",
    )
    chosen_global, _ = select_budget_segments(
        [other_cta, prame_1, prame_2, filler], cancer, cfg_global, ligand_rows=ligand_rows
    )
    assert "filler_1" in {
        s.segment_id for s in chosen_global
    }  # 2 peptides beats 1 under global_union

    cfg_target = replace(cfg_global, allocation_objective="target_specific")
    chosen_target, hist_target = select_budget_segments(
        [other_cta, prame_1, prame_2, filler], cancer, cfg_target, ligand_rows=ligand_rows
    )
    assert "prame_2" in {s.segment_id for s in chosen_target}
    prame2_step = next(h for h in hist_target if h["segment_id"] == "prame_2")
    assert prame2_step["new_alleles"] == 0
    assert prame2_step["new_target_alleles"] == 2
    assert prame2_step["new_monoallelic_target_alleles"] == 2

    support_tbl = target_hla_support_table(
        [other_cta, prame_1, prame_2], cancer, ligand_rows=ligand_rows
    ).set_index("proteoform_key")
    assert support_tbl.loc["prame", "supported_panel_alleles"] == 3
    assert support_tbl.loc["prame", "monoallelic_alleles"] == 3
    assert support_tbl.loc["other", "monoallelic_alleles"] == 1
    assert support_tbl.loc["other", "sample_allele_alleles"] == 1


def test_evidence_tradeoff_tolerance_and_justification_table_202():
    """tsarina#202: explicit min_mortality_gain_pp prevents tiny expression gains
    from automatically displacing MS-rich common-CTA regions, and leave-one-out
    removal/replacement effects are quantified."""
    common_1 = _seg(
        "common_1",
        "ACDEFGHIK",
        proteoform_key="common",
        name="COMMON",
        rank=1,
        score=0.50,
        alleles=("HLA-A*02:01", "HLA-A*01:01"),
        peptides=("ACDEFGHI", "CDEFGHIK", "DEFGHIKL", "EFGHIKLM", "FGHIKLMN"),
    )
    common_2 = _seg(
        "common_2",
        "LMNPQRSTV",
        proteoform_key="common",
        name="COMMON",
        rank=1,
        score=0.50,
        alleles=("HLA-B*07:02",),
        peptides=("LMNPQRST", "MNPQRSTV", "NPQRSTVW", "PQRSTVWY"),
    )
    rare = _seg(
        "rare_1",
        "WYACDEFGH",
        proteoform_key="rare",
        name="RARE",
        rank=2,
        score=0.001,
        alleles=("HLA-A*02:01",),
        peptides=("WYACDEFG",),
    )
    cancer = pd.DataFrame(
        [
            {
                "proteoform_key": "common",
                "burden_category": "lung",
                "prevalence_p95": 0.50,
                "complete_measurement": True,
                "world_mortality_pct": 20.0,
                "world_incidence_pct": 10.0,
            },
            {
                "proteoform_key": "rare",
                "burden_category": "thyroid",
                "prevalence_p95": 0.02,  # Only +0.01 percentage point mortality gain (0.5% * 0.02)!
                "complete_measurement": True,
                "world_mortality_pct": 0.5,
                "world_incidence_pct": 0.5,
            },
        ]
    )
    all_segs = [common_1, common_2, rare]
    # Budget fits 2 segments (1 + 18 = 19 aa).
    cfg_global = VaccineConfig(
        selection_mode="budget",
        max_length_aa=19,
        min_padding=0,
        allocation_objective="global_union",
    )
    sel_global, _ = select_budget_segments(all_segs, cancer, cfg_global)
    assert {s.segment_id for s in sel_global} == {"common_1", "rare_1"}

    cfg_tradeoff = replace(
        cfg_global, allocation_objective="evidence_tradeoff", min_mortality_gain_pp=0.05
    )
    sel_tradeoff, _ = select_budget_segments(all_segs, cancer, cfg_tradeoff)
    assert {s.segment_id for s in sel_tradeoff} == {"common_1", "common_2"}

    # Justification table exposes tiny removal delta of RARE and superior MS replacement by COMMON_2.
    just = allocation_justification_table(sel_global, all_segs, cancer, cfg_global).set_index(
        "name"
    )
    assert just.loc["RARE", "delta_mortality_lower_pp"] == pytest.approx(0.01, abs=1e-4)
    assert just.loc["RARE", "replacement_segments"] == "common_2"
    assert just.loc["RARE", "replacement_delta_ms_peptides"] == 3  # 4 - 1 = +3 distinct MS peptides

    # Tradeoff comparison table compares global_union, target_specific, and tolerances cleanly.
    comp = compare_allocation_tradeoffs(all_segs, cancer, cfg_global, tolerances_pp=(0.0, 0.05))
    tol_row = comp[
        (comp["policy"] == "evidence_tradeoff") & (comp["min_mortality_gain_pp"] == 0.05)
    ].iloc[0]
    global_row = comp[comp["policy"] == "global_union"].iloc[0]
    assert tol_row["distinct_ms_peptides"] > global_row["distinct_ms_peptides"]
    assert tol_row["mortality_loss_vs_max_pp"] == pytest.approx(0.01, abs=1e-4)
