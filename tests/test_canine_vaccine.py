"""Offline scientific acceptance for canine evidence migration (#208)."""

import copy
import json
from dataclasses import replace
from hashlib import sha256
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq

from tsarina.canine_inputs import (
    canine_prevalence,
    dla_alleles,
    load_canine_vaccine_inputs,
    validate_canine_bundle,
)
from tsarina.vaccine import design_vaccine
from tsarina.vaccine_construct import VaccineConfig
from tsarina.vaccine_coverage import paired_target_coverage

FIXTURE = Path(__file__).parent / "fixtures" / "canine-vaccine.json"
A = "DLA-88*001:01"
B = "DLA-88*501:01"
U = "DLA-88*003:02"


def config(**kwargs):
    return replace(
        VaccineConfig.for_canine(
            "osteosarcoma",
            allow_exploratory_dla=True,
            top_k=2,
            lengths=(8, 9),
            min_padding=0,
            max_padding=0,
            optimization_rounds=1,
            max_length_nt=90,
            poly_a_length=6,
        ),
        **kwargs,
    )


def modified_inputs(tmp_path, modify):
    bundle = json.loads(FIXTURE.read_text())
    modify(bundle)
    path = tmp_path / "bundle.json"
    path.write_text(json.dumps(bundle))
    return load_canine_vaccine_inputs(path)


def test_independent_dog_bounds_and_identical_heart_source():
    inputs = load_canine_vaccine_inputs(FIXTURE)
    rank, donors, audit = canine_prevalence(inputs.canine_evidence, "osteosarcoma")
    one = rank[rank.name.eq("DOG_CTA1")].iloc[0]
    assert one.total_identified_dogs == 3
    assert one.positive_dogs == 1 and one.possible_positive_dogs == 2
    assert one.prevalence_lower == pytest.approx(1 / 3)
    assert one.prevalence_upper == pytest.approx(2 / 3)
    assert (
        donors[(donors.proteoform_key.eq(one.proteoform_key)) & donors.donor.eq("dog1")]
        .iloc[0]
        .n_biopsies
        == 2
    )
    assert audit[audit.sample_id.eq("unknown-donor")].iloc[0].inclusion_reason == "unresolved_donor"
    heart = rank[rank.name.eq("DOG_CTA3")].iloc[0]
    assert heart.restriction_status == "rejected"
    assert "heart_alias" in heart.occurrence_ids
    assert "healthy_primary_RNA:heart" in heart.restriction_reason
    assert heart.gene_ids == "g3;g4"
    assert "p95" not in rank.columns and "mortality_weighted_score" not in rank.columns


def test_offline_design_funnel_exact_MS_and_whole_product(tmp_path):
    result = design_vaccine(
        config(), load_canine_vaccine_inputs(FIXTURE), output_dir=tmp_path / "report"
    )
    assert result["status"] == "exploratory_unassessed"
    assert set(result["ranking"].query("selected").name) == {"DOG_CTA1", "DOG_CTA2"}
    assert result["design"]["final_assessment"].startswith("unassessed")
    assert result["design"]["unassessed_alleles"] == [U]
    assert result["design"]["cleavage_model"]["name"] == "unassessed"
    assert all(r["cleavage_probability"] is None for r in result["design"]["cleavage"])
    assert result["design"]["background_safe"]
    dna = result["design"]["coding_sequence"].replace("U", "T")
    assert str(Seq(dna).translate()) == result["design"]["protein"] + "*"
    assert result["design"]["length_nt"] <= 90
    first = result["funnel"].set_index("name").loc["DOG_CTA1"]
    assert first.raw_aa == 16 and first.normal_ms_filtered_aa == 8
    assert first.assembled_pieces == 1 and first.assembled_aa == 8
    rejected = result["source_tables"]["ms_modality_rejected"].set_index("observation_id")
    assert set(rejected.index) == {"binding", "negative", "unknown-result"}
    observed = result["ligands"].query("assembled")
    assert "CDEFGHI" not in set(observed.peptide)
    assert set(observed.query("evidence_tier == 'monoallelic_ms'").allele) == {A}
    assert set(observed.query("evidence_tier == 'sample_allele_ms'").allele) == {B}
    assert set(observed.host_species) == {"Homo sapiens", "Canis lupus familiaris"}
    assert len(set(result["source_tables"]["tumor_genotype_pairs"].iloc[0].alleles)) == 3
    assert result["coverage"]["joint_1_lower"] == pytest.approx(0.8 * 2 / 3)
    assert result["coverage"]["joint_1_upper"] == pytest.approx(1)
    assert result["coverage"]["unsupported_genotype_mass"] == pytest.approx(0.8 * 2 / 3)
    manifest = json.loads((tmp_path / "report" / "manifest.json").read_text())
    assert manifest["counts"]["assembled_proteins"] == 2
    assert manifest["counts"]["construct_length_aa"] == result["design"]["length_aa"]
    assert manifest["counts"]["construct_length_nt"] == result["design"]["length_nt"]
    for file, expected in manifest["artifact_sha256"].items():
        assert sha256((tmp_path / "report" / file).read_bytes()).hexdigest() == expected
    assert "HWE" in (tmp_path / "report" / "report.md").read_text()


def test_no_nested_MS_or_silent_untyped_assignment(tmp_path):
    def modify(b):
        b["ms_hits"] = [r for r in b["ms_hits"] if r["observation_id"] == "long1"]

    result = design_vaccine(config(lengths=(8,)), modified_inputs(tmp_path, modify))
    assert result["design"] is None
    assert result["ligands"].empty
    assert (
        result["source_tables"]["ms_support_decisions"].iloc[0].support_decision
        == "outside_native_specific_intervals"
    )

    def untype(b):
        b["ms_hits"] = [r for r in b["ms_hits"] if r["observation_id"] == "typed1"]
        b["ms_hits"][0].update(sample_alleles=[], restriction_kind="untyped")

    inputs = modified_inputs(tmp_path, untype)
    assert design_vaccine(config(), inputs)["design"] is None
    enabled = design_vaccine(config(allow_untyped_ms=True), inputs)
    assert set(enabled["ligands"].evidence_tier) == {"unrestricted_ms"}
    assert set(enabled["ligands"].restriction_assignment) == {"inferred_by_affinity"}


def test_exploratory_opt_in_preserves_RNA_report(tmp_path):
    result = design_vaccine(
        config(allow_exploratory_dla=False),
        load_canine_vaccine_inputs(FIXTURE),
        output_dir=tmp_path / "rna",
    )
    assert result["design"] is None and not result["ranking"].empty
    assert result["limit"] == "sequence_extrapolated_DLA_requires_explicit_opt_in"
    assert (tmp_path / "rna" / "donor_prevalence.csv").exists()


@pytest.mark.parametrize(
    "change,match",
    [
        (lambda b: b.update(taxon=9606), "taxon 9615"),
        (lambda b: b["occurrences"][0].update(reference_key="other"), "Mixed reference"),
        (lambda b: b["samples"][0].update(unit="counts"), "unit mismatch"),
        (
            lambda b: b["rna_bounds"].append(copy.deepcopy(b["rna_bounds"][0])),
            "Duplicate/unresolved",
        ),
        (lambda b: b["occurrences"][0].update(normal_assessment="unknown"), "Admission requires"),
        (lambda b: b["occurrences"][0].update(sequence_id="bad"), "sequence identity"),
        (lambda b: b["ms_hits"][0].update(mhc_taxon=9606), "independent peptide-source"),
        (lambda b: b["ms_hits"][0].update(sample_alleles=[A, B]), "Monoallelic"),
        (
            lambda b: b["tumor_genotype_pairs"].append(copy.deepcopy(b["tumor_genotype_pairs"][0])),
            "Duplicate pair_id",
        ),
    ],
)
def test_reference_and_evidence_rejections(change, match):
    bundle = json.loads(FIXTURE.read_text())
    change(bundle)
    with pytest.raises(ValueError, match=match):
        validate_canine_bundle(bundle)


def test_species_boundary_and_frozen_identity():
    inputs = load_canine_vaccine_inputs(FIXTURE)
    with pytest.raises(ValueError, match="species='canine'"):
        design_vaccine(VaccineConfig(), inputs)
    with pytest.raises(ValueError, match="frozen canine"):
        design_vaccine(config())
    with pytest.raises(ValueError, match="Human Atlas"):
        config(normal_ms_atlas_dir="human").validate()
    with pytest.raises(ValueError, match="explicit sequences"):
        config(include_utrs=True).validate()
    with pytest.raises(ValueError, match="certification unavailable"):
        config(require_clean_junctions=True).validate()
    with pytest.raises(ValueError, match="coding policy explicitly"):
        design_vaccine(config(codon_species="h_sapiens"), inputs)
    inputs.canine_evidence["sources"]["synthetic"]["version"] = "changed"
    with pytest.raises(ValueError, match="changed after import"):
        design_vaccine(config(), inputs)
    with pytest.raises(ValueError, match="non-class-I"):
        dla_alleles(["HLA-A*02:01"])
    assert dla_alleles([A, B, U]) == [A, B, U]


def test_missing_predictions_fail_closed(tmp_path):
    inputs = modified_inputs(tmp_path, lambda b: b.update(predictions=[]))
    with pytest.raises(ValueError, match="Unassessed DLA affinity"):
        design_vaccine(config(), inputs)


def test_shared_coverage_credits_target_specific_alleles():
    pairs = [
        {
            "pair_id": "x",
            "tumor_donor": "dog",
            "genotype_donor": "dog",
            "alleles": [B],
            "weight": 1,
            "pairing": "observed",
        }
    ]
    targets = {
        "common": {"expression": {"dog": (0, 0)}, "alleles": {B}},
        "expressed": {"expression": {"dog": (1, 1)}, "alleles": {A}},
    }
    assert paired_target_coverage(targets, pairs, {A, B})["joint_1_lower"] == 0
    targets["expressed"]["alleles"].add(B)
    assert paired_target_coverage(targets, pairs, {A, B})["joint_1_lower"] == 1
    assert paired_target_coverage({}, [], set(), missing_mass=1)["joint_1_upper"] == 0


def test_entire_product_background_audit_includes_added_methionine(tmp_path):
    def modify(b):
        source = copy.deepcopy(b["occurrences"][0])
        source.update(
            occurrence_id="initiation-background",
            gene_id="g5",
            transcript_id="t5",
            protein_id="p5",
            gene_symbol="NORMAL",
            sequence="MACDEFGH",
            sequence_id=sha256(b"MACDEFGH").hexdigest(),
            restriction_status="rejected",
            restriction_reason="normal synthetic initiation window",
        )
        b["occurrences"].append(source)

    result = design_vaccine(
        config(top_k=1, optimization_rounds=0), modified_inputs(tmp_path, modify)
    )
    assert result["status"] == "rejected_background_overlap"
    assert result["design"]["background_overlaps"] == [
        {"start": 0, "end": 8, "peptide": "MACDEFGH"}
    ]


def test_cli_defaults_route_frozen_bundle(tmp_path, capsys, monkeypatch):
    import sys

    from tsarina.cli import main

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "tsarina",
            "vaccine",
            "--species",
            "canine",
            "--input-bundle",
            str(FIXTURE),
            "--canine-cohort",
            "osteosarcoma",
            "--allow-exploratory-dla",
            "--top-k",
            "2",
            "--lengths",
            "8,9",
            "--max-padding",
            "0",
            "--optimization-rounds",
            "1",
            "-o",
            str(tmp_path / "cli"),
        ],
    )
    main()
    assert "exploratory_unassessed" in capsys.readouterr().out
    manifest = json.loads((tmp_path / "cli" / "manifest.json").read_text())
    assert manifest["config"]["codon_species"] == "generic"
    assert manifest["config"]["panel"] == "bundle"


def test_verified_saved_canine_report_and_tamper_rejection(tmp_path):
    from tsarina.vaccine_website import render_saved_reports

    source = tmp_path / "source"
    design_vaccine(config(), load_canine_vaccine_inputs(FIXTURE), output_dir=source)
    target = render_saved_reports({"canine": source}, tmp_path / "copy")
    assert target.read_bytes() == (source / "index.html").read_bytes()
    collection = render_saved_reports(
        {"strict": source, "loose demonstration": source}, tmp_path / "collection"
    )
    assert "loose demonstration" in collection.read_text()
    (source / "ranking.csv").write_text("tampered")
    with pytest.raises(ValueError, match="failed verification"):
        render_saved_reports({"canine": source}, tmp_path / "tampered")


def test_empty_DLA_panel_produces_RNA_evidence(tmp_path):
    inputs = modified_inputs(
        tmp_path,
        lambda b: b.update(
            panel=[],
            tumor_genotype_pairs=[],
            genotype_population={
                "name": "unassessed",
                "description": "No genotypes available",
                "cohort": "osteosarcoma",
                "missing_mass": 1,
            },
        ),
    )
    result = design_vaccine(config(), inputs, output_dir=tmp_path / "RNA-only")
    assert result["design"] is None
    assert result["counts"]["selection_enabled_alleles"] == 0
    assert result["coverage"]["missing_mass"] == 1
    assert (tmp_path / "RNA-only" / "donor_prevalence.csv").exists()


def test_cohort_quantification_and_prediction_identity_rejected():
    bundle = json.loads(FIXTURE.read_text())
    bundle["samples"][0]["quantification"] = "gene-level counts"
    with pytest.raises(ValueError, match="quantification mismatch"):
        validate_canine_bundle(bundle)
    bundle = json.loads(FIXTURE.read_text())
    bundle["capabilities"][0]["version"] = "different release"
    with pytest.raises(ValueError, match="capability fingerprint"):
        validate_canine_bundle(bundle)


@pytest.mark.parametrize("kind", ["bulk_proteomics", "II"])
def test_no_bulk_proteomics_or_class_II_support(tmp_path, kind):
    def change(b):
        b["ms_hits"] = [r for r in b["ms_hits"] if r["host_species"] == "Canis lupus familiaris"]
        for row in b["ms_hits"]:
            if kind == "II":
                row.update(mhc_class="II", restriction_kind="untyped", sample_alleles=[])
            else:
                row["assay_context"] = kind

    result = design_vaccine(config(), modified_inputs(tmp_path, change))
    assert result["design"] is None and result["ligands"].empty
    assert "not_class_I_MHC_ligand_elution" in set(
        result["source_tables"]["ms_modality_rejected"].rejection_reason
    )


def test_identical_admitted_occurrences_are_one_target(tmp_path):
    def change(b):
        other = copy.deepcopy(b["occurrences"][0])
        other.update(
            occurrence_id="cta-alias",
            gene_id="alias",
            transcript_id="alias-t",
            protein_id="alias-p",
            gene_symbol="ALIAS",
        )
        b["occurrences"].append(other)

    inputs = modified_inputs(tmp_path, change)
    rank, _, _ = canine_prevalence(inputs.canine_evidence, "osteosarcoma")
    assert len(rank) == 3
    group = rank[rank.name.eq("ALIAS/DOG_CTA1")].iloc[0]
    assert group.prevalence_lower == pytest.approx(1 / 3)
    assert group.gene_ids == "alias;g1"


def test_complete_nucleotide_budget_and_supported_shortfall():
    inputs = load_canine_vaccine_inputs(FIXTURE)
    result = design_vaccine(config(include_utrs=True, utr_5p="ACGT", utr_3p="AACCGG"), inputs)
    design = result["design"]
    assert design["length_nt"] == design["length_aa"] * 3 + 3 + 4 + 6 + 6
    short = design_vaccine(config(selection_mode="supported", max_length_nt=45), inputs)
    assert short["status"] == "insufficient_supported_targets"
    assert short["counts"]["assembled_proteins"] == 1
    assert short["design"]["length_nt"] <= 45


def test_exact_affinity_cutoff_and_runtime_callback_provenance():
    inputs = load_canine_vaccine_inputs(FIXTURE)
    assert design_vaccine(config(ms_affinity_nm=100), inputs)["design"] is None

    def runtime(peptides, alleles):
        return pd.DataFrame(
            [
                {
                    "peptide": p,
                    "allele": a,
                    "affinity_nm": 100 if p in {"ACDEFGHI", "TVWYHGFE", "TVWYHGFED"} else 10000,
                }
                for p in peptides
                for a in alleles
            ],
            columns=["peptide", "allele", "affinity_nm"],
        )

    result = design_vaccine(config(), inputs, affinity_fn=runtime)
    assert result["prediction_provider"] == "injected_affinity_callback"
    assert result["design"]["final_assessment"].startswith("unassessed")


def test_positive_MS_source_decision_precedes_raw_binding_label(tmp_path):
    def change(bundle):
        hit = next(r for r in bundle["ms_hits"] if r["observation_id"] == "typed1")
        hit["qualitative_measurement"] = "Negative binding label"
        hit["assay_response"] = "explicit positive MHC ligand elution detection"
        hit["ms_admission_reason"] = (
            "reviewed positive MS response; raw binding label is independent"
        )
        assert hit["is_ms_observation"] is True

    result = design_vaccine(config(), modified_inputs(tmp_path, change))
    assert "typed1" in set(result["ligands"].observation_id)
    assert "negative" in set(result["source_tables"]["ms_modality_rejected"].observation_id)


def test_portable_reference_cannot_drift():
    inputs = load_canine_vaccine_inputs(FIXTURE)
    inputs.proteins.pop()
    with pytest.raises(ValueError, match="occurrences disagree"):
        design_vaccine(config(), inputs)
    inputs = load_canine_vaccine_inputs(FIXTURE)
    inputs.provenance["reference_key"] = "different"
    with pytest.raises(ValueError, match="provenance disagrees"):
        design_vaccine(config(), inputs)


def test_only_requested_panel_capabilities_set_selection_limit(tmp_path):
    from tsarina.canine_inputs import canonical_hash

    def change(bundle):
        cap = next(c for c in bundle["capabilities"] if c["allele"] == A)
        cap.update(tier="empirically_supported", validation_taxon=9615, validation_sha256="a" * 64)
        bundle["prediction_capability_sha256"] = canonical_hash(bundle["capabilities"])

    result = design_vaccine(
        config(alleles=(A,), allow_exploratory_dla=False, top_k=1),
        modified_inputs(tmp_path, change),
    )
    assert result["design"] is not None
    assert result["limit"] == ""


def test_canine_API_rejects_named_human_panel():
    with pytest.raises(ValueError, match="bundle DLA panel"):
        config(panel="global54_abc").validate()
