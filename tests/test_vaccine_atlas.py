"""Analytical coverage and full-protein source auditing invariants."""

import json
from dataclasses import replace

import pandas as pd
import pytest

from tsarina.regions import allele_frequency_audit
from tsarina.vaccine_budget import select_budget_segments
from tsarina.vaccine_construct import VaccineConfig
from tsarina.vaccine_coverage import carrier_reach, ciwd_frequencies, coverage_tables, union_bounds
from tsarina.vaccine_inputs import filter_ms_modality
from tsarina.vaccine_sequences import Segment
from tsarina.vaccine_tissues import sample_group, tissue_map, verified_normal_ms
from tsarina.vaccine_website import load_saved_report


def test_carrier_probability_sums_same_locus_and_never_normalizes_panel():
    frequencies = {"HLA-A*01:01": 0.1, "HLA-A*02:01": 0.2, "HLA-B*07:02": 0.1}
    result = carrier_reach([*frequencies, "HLA-A*01:01"], frequencies)
    assert result["loci"]["A"] == pytest.approx(0.51)
    assert result["reach"] == pytest.approx(1 - 0.7**2 * 0.9**2)
    assert carrier_reach(["HLA-A*01:01"], frequencies)["reach"] == pytest.approx(0.19)


def test_missing_allele_is_explicit_known_subset_estimate():
    result = carrier_reach(["HLA-A*01:01", "HLA-C*04:03"], {"HLA-A*01:01": 0.1})
    assert result["missing_alleles"] == ["HLA-C*04:03"]
    assert result["reach"] == pytest.approx(0.19)


@pytest.mark.parametrize(
    "frequencies",
    [
        {"HLA-A*01:01": -1},
        {"HLA-A*01:01": float("inf")},
        {"HLA-A*01:01": 0.7, "HLA-A*02:01": 0.5},
    ],
)
def test_invalid_frequency_cannot_be_silently_clipped(frequencies):
    with pytest.raises(ValueError):
        carrier_reach(frequencies, frequencies)


def test_ciwd_primary_table_values_and_coherent_reference():
    values = ciwd_frequencies(
        allele_frequency_audit(["HLA-A*02:01", "HLA-A*01:01", "HLA-C*14:03", "HLA-C*04:03"])
    )
    assert values == {"HLA-A*02:01": 0.24065, "HLA-A*01:01": 0.13343, "HLA-C*14:03": 0.00067}


def test_union_bounds_do_not_assume_patient_independence():
    assert union_bounds([0.3, 0.3]) == pytest.approx((0.3, 0.6, 0))
    assert union_bounds([0.7, 0.8]) == pytest.approx((0.8, 1, 0))
    assert union_bounds([0.3, None]) == pytest.approx((0.3, 1, 1))
    assert union_bounds([]) == (0, 0, 0)
    with pytest.raises(ValueError):
        union_bounds([1.1])


def test_cumulative_evidence_deduplicates_and_upgrades_strongest_tier():
    proteins = pd.DataFrame(
        [
            {"rank": i, "name": key, "proteoform_key": key, "assembled_aa": 10}
            for i, key in enumerate(["p1", "p2"], 1)
        ]
    )
    ligands = pd.DataFrame(
        [
            {
                "proteoform_key": "p1",
                "segment_id": "s1",
                "peptide": "AAAAAAAA",
                "allele": "HLA-A*01:01",
                "evidence_tier": "unrestricted_ms",
                "assembled": True,
            },
            {
                "proteoform_key": "p2",
                "segment_id": "s2",
                "peptide": "AAAAAAAA",
                "allele": "HLA-A*01:01",
                "evidence_tier": "monoallelic_ms",
                "assembled": True,
            },
        ]
    )
    cancer = pd.DataFrame(
        [
            {
                "proteoform_key": p,
                "burden_category": "lung",
                "prevalence_p95": 0.3,
                "complete_measurement": True,
                "world_mortality_pct": 20,
                "world_incidence_pct": 10,
            }
            for p in ["p1", "p2"]
        ]
    )
    layers = [
        {"kind": "cta_segment", "name": p, "proteoform_key": p, "segment_id": s, "end_aa": end}
        for p, s, end in [("p1", "s1", 11), ("p2", "s2", 21)]
    ]
    totals, details = coverage_tables(proteins, ligands, cancer, layers, {"HLA-A*01:01": 0.1})
    final = totals.query("axis == 'protein' and step == 2").iloc[0]
    assert final.ms_peptides == final.pmhc == final.alleles == 1
    assert final.peptides_monoallelic_ms == 1 and final.peptides_unrestricted_ms == 0
    assert final.mortality_lower == pytest.approx(0.06)
    assert final.mortality_upper == pytest.approx(0.12)
    assert totals.query("axis == 'segment'").length_aa.tolist() == [0, 11, 21]
    cancer.loc[0, "complete_measurement"] = False
    totals, details = coverage_tables(proteins, ligands, cancer, layers, {})
    assert details.query("axis == 'protein' and step == 2").missing_proteins.iloc[0] == 1
    assert totals.query("axis == 'protein' and step == 2").mortality_upper.iloc[0] == pytest.approx(
        0.2
    )


def test_explicit_ms_method_overrides_legacy_binding_flag():
    hits = pd.DataFrame(
        [
            {
                "peptide": "ACDEFGHIK",
                "assay_method": "cellular MHC/mass spectrometry",
                "is_binding_assay": True,
            },
            {"peptide": "AAAAAAAAA", "assay_method": "fluorescence", "is_binding_assay": False},
            {
                "peptide": "KRFSLDFNL",
                "assay_method": "cellular MHC/mass spectrometry",
                "qualitative_measurement": "Negative",
                "is_binding_assay": False,
            },
        ]
    )
    accepted, rejected = filter_ms_modality(hits)
    assert accepted.peptide.tolist() == ["ACDEFGHIK"]
    assert rejected.peptide.tolist() == ["AAAAAAAAA", "KRFSLDFNL"]
    assert rejected.rejection_reason.tolist() == ["explicit_non_ms_assay", "negative_assay_result"]


def test_budget_has_no_protein_cap_and_recomputes_marginal_gain():
    segments = [
        Segment(
            f"s{i}",
            f"p{i}",
            f"CTA{i}",
            i,
            0.1,
            "ACDEFGHIK",
            0,
            9,
            0,
            9,
            ("HLA-A*01:01",),
            ("ACDEFGHIK",),
        )
        for i in range(1, 4)
    ]
    cancer = pd.DataFrame(
        [
            {
                "proteoform_key": f"p{i}",
                "burden_category": "lung",
                "prevalence_p95": p,
                "complete_measurement": True,
                "world_mortality_pct": 20,
                "world_incidence_pct": 10,
            }
            for i, p in [(1, 0.8), (2, 0.8), (3, 0.9)]
        ]
    )
    config = VaccineConfig(top_k=1, selection_mode="budget", max_length_aa=30)
    config.validate()
    chosen, history = select_budget_segments(segments, cancer, config)
    assert [s.name for s in chosen] == ["CTA3"]  # Others add no expression or evidence gain.
    segments[1] = replace(segments[1], peptides=("CDEFGHIKL",))
    chosen, history = select_budget_segments(segments, cancer, config)
    assert {s.name for s in chosen} == {"CTA3", "CTA2"}
    assert history[1]["mortality_gain"] == 0 and history[1]["new_ms_peptides"] == 1
    with pytest.raises(ValueError, match="requires"):
        replace(config, max_length_aa=None).validate()
    short, _ = select_budget_segments(segments, cancer, replace(config, max_length_nt=30))
    assert len(short) == 0  # Added M + stop consume the actual nt budget.
    short, _ = select_budget_segments(segments, cancer, replace(config, max_length_nt=33))
    assert len(short) == 1


def test_full_protein_map_preserves_removed_healthy_regions_and_cancer_context():
    proteins = pd.DataFrame(
        [{"proteoform_key": "p", "name": "CTA", "sequence": "ACDEFGHIKLMNPQRSTVWY"}]
    )
    hits = pd.DataFrame(
        [
            {"peptide": "ACDEFGHIK", "src_healthy_tissue": True, "source_tissue": "Heart"},
            {
                "peptide": "MNPQRSTVW",
                "src_cancer": True,
                "source_tissue": "Lung",
                "disease": "lung cancer",
            },
        ]
    )
    intervals = pd.DataFrame([{"proteoform_key": "p", "start": 10, "end": 20}])
    layers = [{"proteoform_key": "p", "native_start": 10, "native_end": 20}]
    mapped = tissue_map(proteins, hits, intervals, layers)
    assert mapped.source_group.tolist() == ["healthy_nonreproductive", "cancer"]
    assert mapped.retained.tolist() == [False, True]
    assert mapped.specific.tolist() == [False, True]
    assert mapped.cancer_only_in_corpus.tolist() == [False, True]
    assert sample_group({"source_tissue": "Lung"}) == "unknown_or_other"
    hits.loc[0, "peptide"] = "MNPQRSTVW"
    mapped = tissue_map(proteins, hits, intervals, layers)
    assert not mapped.cancer_only_in_corpus.any()
    assert mapped.retained.all()  # Normal overlap must not disappear from the map.


def test_saved_report_checks_artifact_hash_and_path_escape(tmp_path):
    manifest = {"artifact_sha256": {"../outside": "x"}}
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="escapes"):
        load_saved_report(tmp_path)
    manifest["artifact_sha256"] = {"protein.fasta": "bad"}
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="verification"):
        load_saved_report(tmp_path)


def test_verified_normal_gate_requires_primary_resolved_donors_and_preserves_typing(monkeypatch):
    rows = pd.DataFrame(
        [
            {
                "peptide": "AETSYVKV",
                "donor_id": donor,
                "donor_status": status,
                "source_tissue": "Heart",
                "tissue_status": tissue,
                "assay_modality": "mass_spectrometry",
                "is_cell_line": line,
                "source_record_id": f"atlas:{i}",
                "mhc_class": "I",
                "sample_alleles": "A*01:01;B*49:01",
                "pmid": 33858848,
            }
            for i, (donor, status, tissue, line) in enumerate(
                [
                    ("AUT01-DN11", "resolved", "nonmalignant", False),
                    ("AUT01-DN11", "resolved", "nonmalignant", False),
                    ("cell-line", "resolved", "nonmalignant", True),
                    ("tumor-adjacent", "resolved", "adjacent", False),
                    ("", "unknown", "nonmalignant", False),
                ]
            )
        ]
    )
    monkeypatch.setattr(
        "hitlist.tissue_blacklist.load_atlas_tissue_evidence", lambda path: (rows, {})
    )
    excluded, observations, summary, audit, provenance = verified_normal_ms(
        "snapshot", {"AETSYVKV"}, 1
    )
    assert summary.n_donors.iloc[0] == 1
    assert provenance["min_donors"] == 1
    assert len(excluded) == 2
    assert not observations.is_monoallelic.any()
    assert observations.sample_alleles.eq("A*01:01;B*49:01").all()
    assert set(audit.exclusion_reason) == {
        "",
        "cell_line",
        "not_nonmalignant_tissue",
        "unresolved_donor",
    }
    excluded, *_ = verified_normal_ms("snapshot", {"AETSYVKV"}, 2)
    assert excluded.empty  # Repeated records from one person are not two donors.
    rows["peptide"] = "SSSSSAETSYVKVVVVVV"  # Longer ligand, absent intact from the CTA.
    excluded, *_ = verified_normal_ms("snapshot", {"AETSYVKV"}, 1)
    assert len(excluded) == 2
    assert excluded.peptide.eq("SSSSSAETSYVKVVVVVV").all()
    with pytest.raises(ValueError, match="verified Atlas"):
        VaccineConfig(normal_ms_policy="exclude").validate()
