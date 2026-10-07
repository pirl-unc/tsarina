"""Scientific invariants for mortality ranking, subtraction and final assembly."""

import json
import sys
from dataclasses import replace
from types import SimpleNamespace

import pandas as pd
import pytest

from tsarina.vaccine import design_vaccine
from tsarina.vaccine_construct import (
    PredictionAudit,
    VaccineConfig,
    encode_construct,
    optimize_construct,
)
from tsarina.vaccine_inputs import (
    Protein,
    VaccineInputs,
    filter_ms_modality,
    load_vaccine_inputs,
    panel_ms_support,
    rank_proteoforms,
    require_current_hitlist,
    resolve_background_cta_ids,
)
from tsarina.vaccine_sequences import (
    Segment,
    assemble_layers,
    junction_windows,
    shared_kmers,
    specific_intervals,
)

ALLELE = "HLA-A*02:01"


def test_existing_positional_vaccine_config_keeps_panel_argument():
    config = VaccineConfig(10, "global54_abc")
    config.validate()
    assert config.panel == "global54_abc" and config.selection_mode == "ranked"


def test_live_design_requires_current_imported_hitlist(monkeypatch):
    monkeypatch.setattr("hitlist.version.__version__", "1.64.6")
    with pytest.raises(ImportError, match=r"Hitlist >=1\.64\.7; imported 1\.64\.6"):
        require_current_hitlist()
    monkeypatch.setattr("hitlist.version.__version__", "1.64.7")
    require_current_hitlist()


def test_alternate_haplotypes_inherit_gene_identity_but_independent_loci_do_not():
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    proteins = [
        Protein("cta", "PRAME", "p", sequence, contig="22"),
        Protein("alt", "PRAME", "p_alt", sequence, contig="HSCHR22_1_CTG3"),
        Protein("other", "OTHER", "p_other", sequence[:8], contig="1"),
        Protein("unknown_alt", "OTHER", "p_unknown", sequence[-8:], contig="HSCHR1_1_CTG"),
    ]
    ids, aliases = resolve_background_cta_ids(proteins, {"cta"})
    assert ids == {"cta", "alt"}
    assert aliases.iloc[0].cta_reference_gene_ids == "cta"
    assert shared_kmers(proteins[:2], ids, [sequence]) == set()
    assert shared_kmers(proteins, ids, [sequence]) == {sequence[:8], sequence[-8:]}
    # A same-symbol primary locus is never admitted by the alias resolver.
    independent = Protein("normal", "PRAME", "normal_p", sequence, contig="1")
    all_ids, _ = resolve_background_cta_ids([*proteins, independent], {"cta"})
    assert "normal" not in all_ids
    assert specific_intervals(sequence, shared_kmers([independent], all_ids, [sequence])) == []


@pytest.mark.parametrize("percentile", [float("nan"), float("inf"), -1, 101])
def test_invalid_ms_prediction_cannot_become_an_evidence_dropout(inputs, monkeypatch, percentile):
    def scores(peptides, alleles, predictor):
        return pd.DataFrame(
            [
                {"peptide": p, "allele": a, "presentation_percentile": percentile}
                for p in peptides
                for a in alleles
            ]
        )

    monkeypatch.setattr("tsarina.scoring.score_presentation", scores)
    with pytest.raises(ValueError, match="Invalid MS presentation percentiles"):
        panel_ms_support({"ACDEFGHIK"}, [ALLELE], inputs.ms_hits)


@pytest.fixture
def inputs():
    seq = "MMMMMMMMACDEFGHIKLLLLLLLL"
    return VaccineInputs(
        proteins=[
            Protein("g1", "CTA1", "p1", seq),
            Protein("g2", "CTA2", "p2", seq),
            Protein("g3", "CTA3", "p3", "MMMMMMMMTVWYACDEFGLLLLLLLL"),
            Protein("g4", "OTHER", "p4", "MMMMMMMMLLLLLLLL", canonical=False),
        ],
        cta_gene_ids={"g1", "g2", "g3"},
        gene_keys={"g1": "group", "g2": "group", "g3": "g3"},
        prevalence=pd.DataFrame(
            [
                {"proteoform_key": key, "cancer_code": code, "prevalence_p95": val, "n_samples": n}
                for key, code, val, n in [
                    ("group", "L1", 0.9, 10),
                    ("group", "L2", 0.3, 30),
                    ("g3", "L1", 0.2, 10),
                    ("g3", "L2", 0.2, 30),
                ]
            ]
        ),
        burden=pd.DataFrame(
            [
                {
                    "burden_category": "lung",
                    "world_mortality_pct": 20.0,
                    "world_incidence_pct": 12.0,
                    "world_mortality_count": float("nan"),
                }
            ]
        ),
        cohorts={"lung": ["L1", "L2"]},
        ms_hits=pd.DataFrame(
            [
                {
                    "peptide": "ACDEFGHIK",
                    "mhc_restriction": ALLELE,
                    "is_monoallelic": True,
                    "pmid": "123",
                    "assay_modality": "mass_spectrometry",
                }
            ]
        ),
        provenance={"synthetic": True},
    )


def affinity(peptides, alleles):
    return pd.DataFrame(
        [{"peptide": p, "allele": a, "affinity_nm": 5000.0} for p in peptides for a in alleles]
    )


def cleavage(sequences):
    return {s: [0.5] * len(s) for s in sequences}


def segment(key, sequence, ligand_start=0, ligand_end=None, rank=1):
    return Segment(
        key,
        key,
        key,
        rank,
        1 / rank,
        sequence,
        0,
        len(sequence),
        ligand_start,
        len(sequence) if ligand_end is None else ligand_end,
        (ALLELE,),
        (sequence[ligand_start:ligand_end],),
    )


@pytest.fixture
def codons(monkeypatch):
    # Independent genetic-code backtranslation; production uses Vaxrank/DnaChisel.
    from Bio.Data.CodonTable import standard_dna_table

    table = {aa: c for c, aa in standard_dna_table.forward_table.items()}
    monkeypatch.setattr(
        "tsarina.vaccine_elements.codon_optimize",
        lambda protein, species: "".join(table[aa] for aa in protein),
    )


@pytest.fixture
def ms_scores(monkeypatch):
    def scores(peptides, alleles, predictor):
        return pd.DataFrame(
            [
                {
                    "peptide": p,
                    "allele": a,
                    "presentation_percentile": 0.2,
                    "presentation_score": 0.9,
                    "affinity_nm": 50.0,
                }
                for p in peptides
                for a in alleles
            ]
        )

    monkeypatch.setattr("tsarina.scoring.score_presentation", scores)


def test_rank_identical_proteins_once_with_weighted_cohort_denominators(inputs):
    inputs.proteins.append(Protein("g1", "CTA1", "p1_alt", inputs.proteins[0].sequence, "t_alt"))
    ranking, detail, summary = rank_proteoforms(inputs)
    assert ranking.name.tolist() == ["CTA1/CTA2", "CTA3"]
    assert ranking.iloc[0].mortality_weighted_score == pytest.approx(
        0.2 * (0.9 * 10 + 0.3 * 30) / 40
    )
    assert len(detail) == 4
    assert set(ranking.iloc[0].protein_ids.split(";")) == {"p1", "p1_alt", "p2"}
    assert "t_alt" in ranking.iloc[0].transcript_ids.split(";")
    assert summary.iloc[0].n_samples == 40
    assert pd.isna(summary.iloc[0].world_mortality_count)


def test_missing_prevalence_is_explicit_not_asserted_zero(inputs):
    inputs.prevalence.loc[1, "prevalence_p95"] = float("nan")
    rank, detail, summary = rank_proteoforms(inputs)
    group = rank.set_index("proteoform_key").loc["group"]
    assert not group.complete_measurement
    assert detail.query("measurement_status == 'missing'").shape[0] == 1
    assert summary.query("proteoform_key == 'group'").iloc[0].n_measured_samples == 10


def test_sequence_registry_disagreement_refuses_percentile_arithmetic(inputs):
    inputs.gene_keys["g2"] = "g2"
    with pytest.raises(ValueError, match="prevalences cannot be added"):
        rank_proteoforms(inputs)


@pytest.mark.parametrize(
    "change,error",
    [
        (lambda x: x.cohorts.update({"breast": ["L1"]}), "counted twice"),
        (lambda x: x.prevalence.__setitem__("prevalence_p95", 1.1), "fraction"),
        (lambda x: x.prevalence.__setitem__("n_samples", 0), "denominators"),
    ],
)
def test_invalid_scientific_inputs_fail(inputs, change, error):
    change(inputs)
    with pytest.raises(ValueError, match=error):
        rank_proteoforms(inputs)


def test_subtraction_masks_entire_overlapping_shared_eightmers():
    seq = "ACDEFGHIKLMNPQRSTVWY"
    forbidden = {seq[2:10], seq[7:15]}
    assert specific_intervals(seq, forbidden) == [(0, 2), (15, 20)]
    assert specific_intervals("ACDXEFG", set()) == [(0, 3), (4, 7)]


def test_all_noncta_isoforms_block_but_withincta_sharing_does_not(inputs):
    seq = inputs.proteins[0].sequence
    blocked = shared_kmers(inputs.proteins, inputs.cta_gene_ids, [seq])
    assert specific_intervals(seq, blocked) == [(8, 17)]
    inputs.proteins.append(Protein("g5", "OTHER2", "p5", seq, canonical=False))
    blocked = shared_kmers(inputs.proteins, inputs.cta_gene_ids, [seq])
    assert specific_intervals(seq, blocked) == []


def test_ms_assignment_separates_measured_and_inferred_restrictions(ms_scores):
    other = "HLA-A*24:02"
    hits = pd.DataFrame(
        [
            {"peptide": "ACDEFGHIK", "mhc_restriction": ALLELE, "is_monoallelic": True},
            {"peptide": "TVWYACDEF", "mhc_restriction": "", "is_monoallelic": False},
        ]
    )
    hits["assay_modality"] = "mass_spectrometry"
    result = panel_ms_support(set(hits.peptide), [ALLELE, other], hits)
    assert set(result.query("peptide == 'ACDEFGHIK'").allele) == {ALLELE}
    assert result.query("peptide == 'ACDEFGHIK'").iloc[0].evidence_tier == "monoallelic_ms"
    assert set(result.query("peptide == 'TVWYACDEF'").evidence_tier) == {"unrestricted_ms"}


def test_positive_ms_modality_rejects_nonbinding_fluorescence_and_structures():
    hits = pd.DataFrame(
        [
            {"peptide": "a", "assay_method": "cellular MHC/mass spectrometry"},
            {"peptide": "b", "assay_method": "", "source": "supplement"},
            {"peptide": "c", "assay_modality": "mass_spectrometry"},
            {
                "peptide": "d",
                "assay_method": "cellular MHC/direct/fluorescence",
                "assay_modality": "mass_spectrometry",
                "is_binding_assay": False,
            },
            {"peptide": "e", "assay_method": "x-ray crystallography", "is_binding_assay": False},
            {"peptide": "f", "assay_method": "purified MHC/direct/fluorescence"},
            {"peptide": "g"},
        ]
    )
    accepted, rejected = filter_ms_modality(hits)
    assert set(accepted.peptide) == {"a", "b", "c"}
    assert set(rejected.peptide) == {"d", "e", "f", "g"}
    assert rejected.set_index("peptide").loc["g", "rejection_reason"] == "unknown_assay_modality"


def test_sample_affinity_accepts_every_binding_sample_allele_without_percentile_gate(monkeypatch):
    b, outside = "HLA-A*24:02", "HLA-B*07:02"

    def score(peptides, alleles, predictor):
        return pd.DataFrame(
            [
                {
                    "peptide": p,
                    "allele": a,
                    "affinity_nm": {ALLELE: 100, b: 800, outside: 10}[a],
                    "presentation_percentile": 99.0,
                }
                for p in peptides
                for a in alleles
            ]
        )

    monkeypatch.setattr("tsarina.scoring.score_presentation", score)
    hits = pd.DataFrame(
        [
            {
                "peptide": "ACDEFGHIK",
                "assay_method": "mass spectrometry",
                "mhc_allele_provenance": "sample_allele_match",
                "mhc_allele_set": f"{ALLELE};{b}",
                "provenance_id": "sample:1",
            }
        ]
    )
    result = panel_ms_support(set(hits.peptide), [ALLELE, b, outside], hits, mode="sample_affinity")
    assert set(result.allele) == {ALLELE, b}
    assert set(result.evidence_tier) == {"sample_allele_ms"}
    assert len(result.attrs["ms_assignments"]) == 2
    assert set(result.attrs["ms_assignments"].provenance_id) == {"sample:1"}


def test_untyped_sample_affinity_is_explicit_and_study_pool_is_not_sample_typing(monkeypatch):
    def score(peptides, alleles, predictor):
        return pd.DataFrame(
            [{"peptide": p, "allele": a, "affinity_nm": 999.0} for p in peptides for a in alleles]
        )

    monkeypatch.setattr("tsarina.scoring.score_presentation", score)
    hits = pd.DataFrame(
        [
            {
                "peptide": "ACDEFGHIK",
                "assay_method": "mass spectrometry",
                "mhc_allele_provenance": "pmid_class_pool",
                "mhc_allele_set": ALLELE,
            }
        ]
    )
    assert panel_ms_support(set(hits.peptide), [ALLELE], hits, mode="sample_affinity").empty
    result = panel_ms_support(
        set(hits.peptide), [ALLELE], hits, mode="sample_affinity", allow_untyped=True
    )
    assert result.iloc[0].evidence_tier == "unrestricted_ms"
    assert result.attrs["ms_assignments"].iloc[0].sample_alleles == ""
    # A different nested peptide is not an observation of the requested peptide.
    assert panel_ms_support(
        {"CDEFGHIK"}, [ALLELE], hits, mode="sample_affinity", allow_untyped=True
    ).empty


@pytest.mark.parametrize("bad", [0, -1, float("nan"), float("inf")])
def test_sample_affinity_rejects_missing_or_invalid_predictions(inputs, monkeypatch, bad):
    monkeypatch.setattr(
        "tsarina.scoring.score_presentation",
        lambda p, a, predictor: pd.DataFrame(
            [{"peptide": x, "allele": y, "affinity_nm": bad} for x in p for y in a]
        ),
    )
    with pytest.raises(ValueError, match="Invalid MS affinity"):
        panel_ms_support({"ACDEFGHIK"}, [ALLELE], inputs.ms_hits, mode="sample_affinity")


def test_ms_support_validates_cutoff_even_without_observations():
    with pytest.raises(ValueError, match="cutoff must be finite"):
        panel_ms_support(
            set(), [ALLELE], pd.DataFrame(), mode="sample_affinity", affinity_nm=float("nan")
        )


def test_supported_selection_backfills_and_records_the_rejection(inputs, ms_scores, codons):
    inputs.ms_hits.loc[0, "peptide"] = "TVWYACDEF"
    result = design_vaccine(
        VaccineConfig(top_k=1, selection_mode="supported", alleles=(ALLELE,)),
        inputs,
        affinity_fn=affinity,
        cleavage_fn=cleavage,
    )
    assert list(result["funnel"].name) == ["CTA3"]
    screen = result["source_tables"]["selection_screen"].set_index("name")
    assert screen.loc["CTA1/CTA2", "selection_reason"] == "no_qualified_panel_ms_ligands"
    assert screen.loc["CTA3", "selected"]


def test_exclusion_vetoes_any_group_member_without_changing_cta_background(
    inputs, ms_scores, codons
):
    inputs.proteins[0] = replace(inputs.proteins[0], gene_name="MAGEA4")
    inputs.proteins[1] = replace(inputs.proteins[1], gene_name="MAGEA3")
    inputs.ms_hits.loc[0, "peptide"] = "TVWYACDEF"
    result = design_vaccine(
        VaccineConfig(
            top_k=1, alleles=(ALLELE,), exclude_gene_patterns=("MAGE*",), allow_genes=("MAGEA4",)
        ),
        inputs,
        affinity_fn=affinity,
        cleavage_fn=cleavage,
    )
    assert list(result["funnel"].name) == ["CTA3"]
    row = result["ranking"].query("name == 'MAGEA3/MAGEA4'").iloc[0]
    assert not row.eligible and row.selection_reason == "excluded_gene_pattern"
    assert result["provenance"]["background_translated_occurrences"] == 1


def test_supported_selection_continues_after_an_empty_first_batch(inputs, ms_scores, codons):
    inputs.proteins = [
        Protein(f"g{i}", f"CTA{i}", f"p{i}", "M" * i + "ACDEFGHIK") for i in range(1, 11)
    ]
    inputs.proteins.append(Protein("g11", "CTA11", "p11", "MTVWYACDEF"))
    inputs.cta_gene_ids = {p.gene_id for p in inputs.proteins}
    inputs.gene_keys = {g: g for g in inputs.cta_gene_ids}
    inputs.prevalence = pd.DataFrame(
        [
            {
                "proteoform_key": f"g{i}",
                "cancer_code": c,
                "prevalence_p95": 1 - i * 0.03,
                "n_samples": n,
            }
            for i in range(1, 12)
            for c, n in [("L1", 10), ("L2", 30)]
        ]
    )
    inputs.ms_hits.loc[0, "peptide"] = "TVWYACDEF"
    result = design_vaccine(
        VaccineConfig(top_k=1, selection_mode="supported", alleles=(ALLELE,)),
        inputs,
        affinity_fn=affinity,
        cleavage_fn=cleavage,
    )
    assert list(result["funnel"].name) == ["CTA11"]
    assert result["provenance"]["inspected_proteoforms"] == 11


def test_supported_shortfall_writes_audit_instead_of_claiming_target_count(
    inputs, ms_scores, codons, tmp_path
):
    with pytest.raises(ValueError, match="Only 1 supported proteoforms fit; requested 2"):
        design_vaccine(
            VaccineConfig(top_k=2, selection_mode="supported", alleles=(ALLELE,)),
            inputs,
            output_dir=tmp_path,
            affinity_fn=affinity,
            cleavage_fn=cleavage,
        )
    assert (
        json.loads((tmp_path / "manifest.json").read_text())["status"]
        == "insufficient_supported_proteoforms"
    )
    assert (tmp_path / "selection_screen.csv").exists()


def test_supported_zero_ms_writes_empty_selection_audit(inputs, tmp_path):
    inputs.ms_hits.loc[0, "assay_method"] = "cellular MHC/direct/fluorescence"
    with pytest.raises(ValueError, match="Only 0 supported proteoforms fit"):
        design_vaccine(
            VaccineConfig(top_k=1, selection_mode="supported", alleles=(ALLELE,)),
            inputs,
            output_dir=tmp_path,
            affinity_fn=affinity,
            cleavage_fn=cleavage,
        )
    assert (
        json.loads((tmp_path / "manifest.json").read_text())["status"]
        == "insufficient_supported_proteoforms"
    )
    assert (
        pd.read_csv(tmp_path / "rejected_ms_observations.csv").iloc[0].assay_method
        == "cellular MHC/direct/fluorescence"
    )


def test_supported_construct_reserves_one_piece_per_protein_at_length_cap():
    a1 = segment("a1", "CCCCCCCCCCC")
    a2 = replace(segment("a2", "DDDDDDDDDDD"), proteoform_key=a1.proteoform_key)
    b = segment("b", "EEEEEEEEEEEE", rank=2)
    result = optimize_construct(
        [a1, a2, b],
        [ALLELE],
        VaccineConfig(
            top_k=2,
            selection_mode="supported",
            max_length_aa=24,
            max_padding=0,
            optimization_rounds=0,
        ),
        affinity,
        cleavage,
    )
    keys = {r["proteoform_key"] for r in result["layers"] if r["kind"] == "cta_segment"}
    assert keys == {"a1", "b"}
    assert len(result["protein"]) == 24
    assert len(result["excluded_segments"]) == 1


def test_junction_audit_includes_multiboundary_and_synthetic_start():
    a, b = segment("a", "ACDEFGHI"), segment("b", "KLMNPQRS")
    seq, layers = assemble_layers([(a, 0, 0, ""), (b, 0, 0, "AAY")])
    windows = junction_windows(seq, layers)
    assert seq.startswith("MACDE")
    assert any(w["start"] == 0 for w in windows)
    assert any(";" in w["boundaries"] for w in windows)
    assert all(len(w["peptide"]) in {8, 9, 10, 11} for w in windows)


def test_order_search_removes_binder_and_preserves_ligands():
    a, b = segment("a", "MCCCCCCCC"), segment("b", "DDDDDDDDD", rank=2)

    def scorer(peptides, alleles):
        frame = affinity(peptides, alleles)
        frame.loc[frame.peptide.str.contains("CD"), "affinity_nm"] = 100.0
        return frame

    result = optimize_construct(
        [a, b], [ALLELE], VaccineConfig(max_padding=0, optimization_rounds=1), scorer, cleavage
    )
    assert result["clean_junctions"]
    assert result["layers"][1]["segment_id"] == "b"
    assert a.protein_sequence in result["protein"] and b.protein_sequence in result["protein"]


def test_clamped_padding_keeps_beam_slots_for_distinct_constructs():
    a, b = segment("a", "MCCCCCCCC"), segment("b", "DDDDDDDDD", rank=2)
    result = optimize_construct(
        [a, b],
        [ALLELE],
        VaccineConfig(padding_step=2, optimization_rounds=2),
        affinity,
        cleavage,
    )
    # No terminal context exists: there are only two orders and two joins.
    # Requested padding 0..10 must not create additional construct states.
    assert all(row["candidates"] <= 4 for row in result["search_history"])
    assert result["clean_junctions"]


def test_padding_changes_only_unsupported_edges():
    a, b = segment("a", "MCCCCCCCCQ", ligand_end=9), segment("b", "DDDDDDDD", rank=2)

    def scorer(peptides, alleles):
        frame = affinity(peptides, alleles)
        frame.loc[frame.peptide.str.contains("Q"), "affinity_nm"] = 100.0
        return frame

    result = optimize_construct(
        [a, b],
        [ALLELE],
        VaccineConfig(max_padding=1, padding_step=1, optimization_rounds=1),
        scorer,
        cleavage,
    )
    assert result["clean_junctions"]
    assert "MCCCCCCCC" in result["protein"]
    assert "Q" not in result["protein"]


def test_aay_rescue_versus_direct_joins():
    a, b = segment("a", "MCCCCCCCC"), segment("b", "DDDDDDDDD", rank=2)

    def scorer(peptides, alleles):
        frame = affinity(peptides, alleles)
        frame.loc[frame.peptide.str.contains("CD|DC|DM"), "affinity_nm"] = 100.0
        return frame

    result = optimize_construct(
        [a, b], [ALLELE], VaccineConfig(max_padding=0, optimization_rounds=1), scorer, cleavage
    )
    assert result["clean_junctions"]
    assert any(layer["sequence"] == "AAY" for layer in result["layers"])


def test_linker_can_trim_padding_in_same_step_at_length_cap():
    a = segment("a", "MCCCCCCCCQQQ", ligand_end=9)
    b = segment("b", "DDDDDDDDD", rank=2)

    def scorer(peptides, alleles):
        frame = affinity(peptides, alleles)
        frame.loc[frame.peptide.str.contains("QD"), "affinity_nm"] = 100.0
        frame.loc[frame.peptide.str.contains("CD|DM"), "affinity_nm"] = 1.0
        frame.loc[frame.peptide.str.contains("AAY"), "affinity_nm"] = 10000.0
        return frame

    cfg = VaccineConfig(
        max_padding=3,
        padding_step=3,
        beam_width=1,
        optimization_rounds=1,
        max_length_aa=21,
    )
    result = optimize_construct([a, b], [ALLELE], cfg, scorer, cleavage)
    assert result["clean_junctions"]
    assert result["protein"] == "MCCCCCCCCAAYDDDDDDDDD"


def test_limits_include_start_stop_and_polya(codons):
    a, b = segment("a", "ACDEFGHIK"), segment("b", "TVWYACDEF", rank=2)
    cfg = VaccineConfig(max_padding=0, max_length_aa=10, max_length_nt=43, poly_a_length=10)
    result = optimize_construct([a, b], [ALLELE], cfg, affinity, cleavage)
    assert len(result["protein"]) == 10
    assert result["excluded_segments"] == [{"segment_id": "b", "reason": "construct_length_limit"}]
    cds, full = encode_construct(result, cfg)
    assert len(cds) == 33 and len(full) == 43
    assert "T" not in full and full.endswith("A" * 10)
    dna, _ = encode_construct(result, replace(cfg, vaccine_type="dna"))
    assert "U" not in dna and dna.endswith("TAA")


def test_predictor_missing_pairs_cannot_look_clean():
    audit = PredictionAudit(
        [ALLELE, "HLA-A*24:02"], VaccineConfig(), lambda p, a: affinity(p, a[:1]), cleavage
    )
    a = segment("a", "ACDEFGHIK")
    with pytest.raises(ValueError, match="Missing junction affinity"):
        audit.prepare([assemble_layers([(a, 0, 0, "")])])


def test_full_pipeline_artifacts_and_dropout(inputs, ms_scores, codons, tmp_path):
    cfg = VaccineConfig(top_k=2, alleles=(ALLELE,), max_padding=0)
    result = design_vaccine(
        cfg, inputs, output_dir=tmp_path, affinity_fn=affinity, cleavage_fn=cleavage
    )
    funnel = result["funnel"].set_index("name")
    assert funnel.loc["CTA1/CTA2", "assembled_aa"] == 9
    assert funnel.loc["CTA3", "status"] == "no_qualified_panel_ms_ligands"
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    from tsarina.vaccine_website import render_saved_reports

    site = json.loads((tmp_path / "website" / "data.json").read_text())
    view = site["designs"][cfg.definition]
    assert view["panel_size"] == 1
    assert view["protein_aa"] == len(result["design"]["protein"])
    assert view["sequences"]["full"] == result["design"]["full_nt"]
    index = render_saved_reports(
        {"custom": tmp_path}, tmp_path / "rendered", analysis_date="2026-10-07"
    )
    replay = json.loads(index.with_name("data.json").read_text())
    assert replay["analysis_date"] == "2026-10-07"
    assert replay["designs"]["custom"]["coverage"] == view["coverage"]
    assert manifest["provenance"]["synthetic"]
    assert manifest["provenance"]["ms_input_kind"] == "supplied_observations"
    assert "MS evidence source: **supplied observations**" in (tmp_path / "report.md").read_text()
    assert manifest["prediction_callbacks"] == {"affinity": True, "cleavage": True}
    assert manifest["cancer_summary"][0]["world_mortality_count"] is None
    assert {
        "report.md",
        "protein.fasta",
        "cds.fasta",
        "full.fasta",
        "funnel.csv",
        "hla_panel.csv",
        "sequence-funnel.svg",
        "cancer-priorities.svg",
        "construct-map.svg",
    } <= set(manifest["artifact_sha256"])
    with pytest.raises(ValueError, match="already contains"):
        design_vaccine(cfg, inputs, output_dir=tmp_path, affinity_fn=affinity, cleavage_fn=cleavage)


def test_no_ms_support_writes_explanatory_audit(inputs, tmp_path):
    inputs.ms_hits = inputs.ms_hits.iloc[:0]
    with pytest.raises(ValueError, match="No selected CTA retains"):
        design_vaccine(
            VaccineConfig(alleles=(ALLELE,)),
            inputs,
            output_dir=tmp_path,
            affinity_fn=affinity,
            cleavage_fn=cleavage,
        )
    assert (
        json.loads((tmp_path / "manifest.json").read_text())["status"]
        == "no_ms_supported_construct"
    )


def test_impossible_length_writes_funnel(inputs, ms_scores, tmp_path):
    with pytest.raises(ValueError, match="length constraints"):
        design_vaccine(
            VaccineConfig(alleles=(ALLELE,), max_length_nt=5),
            inputs,
            output_dir=tmp_path,
            affinity_fn=affinity,
            cleavage_fn=cleavage,
        )
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["status"] == "no_feasible_construct"
    assert manifest["funnel"][0]["status"] == "construct_length_limit"
    assert not (tmp_path / "full.fasta").exists()


def test_empty_ms_input_and_duplicate_alleles(inputs, ms_scores, codons):
    assert panel_ms_support({"ACDEFGHIK"}, [ALLELE], pd.DataFrame()).empty
    result = design_vaccine(
        VaccineConfig(alleles=(ALLELE, ALLELE)), inputs, affinity_fn=affinity, cleavage_fn=cleavage
    )
    assert result["alleles"] == [ALLELE]
    ligand = result["ligands"].iloc[0]
    assert (
        result["design"]["protein"][ligand.construct_start : ligand.construct_end] == ligand.peptide
    )


def test_strict_unresolved_junctions_fail_with_audit(inputs, ms_scores, codons, tmp_path):
    def strong(peptides, alleles):
        frame = affinity(peptides, alleles)
        frame["affinity_nm"] = 20.0
        return frame

    with pytest.raises(ValueError, match="unresolved junction binders"):
        design_vaccine(
            VaccineConfig(alleles=(ALLELE,), require_clean_junctions=True),
            inputs,
            output_dir=tmp_path,
            affinity_fn=strong,
            cleavage_fn=cleavage,
        )
    assert not json.loads((tmp_path / "manifest.json").read_text())["design"]["clean_junctions"]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_padding": -1},
        {"min_padding": 11},
        {"max_length_nt": 0},
        {"junction_affinity_nm": float("nan")},
        {"lengths": (7,)},
        {"shared_k": 9},
        {"linkers": ("X",)},
    ],
)
def test_invalid_design_constraints(kwargs):
    with pytest.raises(ValueError):
        VaccineConfig(**kwargs).validate()


def test_vaccine_cli_both_definitions(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "tsarina.vaccine_website.render_saved_reports", lambda reports, output: output
    )
    from tsarina.cli import main

    calls = []

    def design(config, **kwargs):
        calls.append((config, kwargs))
        return {"design": {"length_aa": 20, "length_nt": 63}}

    monkeypatch.setattr("tsarina.vaccine.design_vaccine", design)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "tsarina",
            "vaccine",
            "--cta-definition",
            "both",
            "--hla",
            ALLELE,
            "--vaccine-type",
            "mrna",
            "--max-length-nt",
            "100",
            "--selection-mode",
            "supported",
            "--exclude-gene-pattern",
            "MAGE*",
            "--allow-gene",
            "MAGEA4",
            "--ms-support-mode",
            "sample-affinity",
            "--allow-untyped-ms",
            "-o",
            str(tmp_path),
        ],
    )
    main()
    assert [c.definition for c, _ in calls] == ["strict", "loose"]
    assert all(c.vaccine_type == "rna" and c.max_length_nt == 100 for c, _ in calls)
    assert calls[1][1]["output_dir"] == tmp_path / "loose"
    assert all(
        c.selection_mode == "supported"
        and c.exclude_gene_patterns == ("MAGE*",)
        and c.allow_genes == ("MAGEA4",)
        and c.ms_support_mode == "sample_affinity"
        and c.allow_untyped_ms
        for c, _ in calls
    )


def test_loader_collapses_genome_before_percentiles_and_keeps_repeated_sources(inputs, monkeypatch):
    from oncoref import cta, expression, proteoforms

    calls = []
    monkeypatch.setattr("tsarina.vaccine_inputs.require_current_hitlist", lambda: None)

    def prevalence(code, **kwargs):
        calls.append(kwargs)
        return inputs.prevalence.query("cancer_code == @code").rename(
            columns={"prevalence_p95": "frac_samples_top5pct"}
        )

    monkeypatch.setattr(cta, "cta_gene_ids", lambda: inputs.cta_gene_ids)
    monkeypatch.setattr(
        proteoforms,
        "proteoform_key",
        lambda gene, scope: inputs.gene_keys[gene] if scope == "genome" else gene,
    )
    monkeypatch.setattr(expression, "proteoform_within_sample_top_fraction", prevalence)
    genes = []
    for p in inputs.proteins:
        transcript = SimpleNamespace(
            id=p.transcript_id,
            protein_id="shared_protein_id",
            protein_sequence=p.sequence,
            biotype="protein_coding",
        )
        genes.append(
            SimpleNamespace(
                id=p.gene_id, name=p.gene_name, biotype="protein_coding", transcripts=[transcript]
            )
        )
    monkeypatch.setattr(
        "pyensembl.EnsemblRelease", lambda release: SimpleNamespace(genes=lambda: genes)
    )
    result = load_vaccine_inputs(cohorts=inputs.cohorts)
    assert all(c["scope"] == "genome" and c["threshold"] == 0.95 for c in calls)
    assert len(result.proteins) == len(inputs.proteins)
    assert {p.gene_id for p in result.proteins if p.protein_id == "shared_protein_id"} == {
        p.gene_id for p in inputs.proteins
    }


def test_both_continues_after_one_definition_fails(monkeypatch, tmp_path):
    from tsarina.cli import main

    definitions = []

    def design(config, **kwargs):
        definitions.append(config.definition)
        if config.definition == "strict":
            raise ValueError("No supported strict sequence")
        return {"design": {"length_aa": 20, "length_nt": 63}}

    monkeypatch.setattr("tsarina.vaccine.design_vaccine", design)
    monkeypatch.setattr(
        sys, "argv", ["tsarina", "vaccine", "--cta-definition", "both", "-o", str(tmp_path)]
    )
    with pytest.raises(SystemExit) as error:
        main()
    assert error.value.code == 1
    assert definitions == ["strict", "loose"]
