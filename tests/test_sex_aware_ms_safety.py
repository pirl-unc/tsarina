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

"""Tests for sex-aware MS safety filtering using Hitlist female/male reproductive flags."""

from __future__ import annotations

import argparse
from unittest.mock import patch

import pandas as pd

import tsarina
from tsarina import cli_hits
from tsarina import personalize as pmod
from tsarina.ms_attribution import attribute_ms_evidence, summarize_ms_attribution
from tsarina.ms_evidence import aggregate_ms_hits_by_peptide


def _mock_observations() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "peptide": "OVARYPEP1",
                "mhc_restriction": "HLA-A*02:01",
                "mhc_species": "Homo sapiens",
                "is_ms_observation": True,
                "source_tissue": "Ovary",
                "src_cancer": True,
                "src_healthy_tissue": False,
                "src_healthy_reproductive": False,
                "src_healthy_reproductive_female": True,
                "src_healthy_reproductive_male": False,
                "pmid": "111",
            },
            {
                "peptide": "PROSTPEP1",
                "mhc_restriction": "HLA-A*02:01",
                "mhc_species": "Homo sapiens",
                "is_ms_observation": True,
                "source_tissue": "Prostate",
                "src_cancer": True,
                "src_healthy_tissue": False,
                "src_healthy_reproductive": False,
                "src_healthy_reproductive_female": False,
                "src_healthy_reproductive_male": True,
                "pmid": "222",
            },
            {
                "peptide": "TESTISPEP",
                "mhc_restriction": "HLA-A*02:01",
                "mhc_species": "Homo sapiens",
                "is_ms_observation": True,
                "source_tissue": "Testis",
                "src_cancer": True,
                "src_healthy_tissue": False,
                "src_healthy_reproductive": True,
                "src_healthy_reproductive_female": False,
                "src_healthy_reproductive_male": False,
                "pmid": "333",
            },
            {
                "peptide": "HEARTPEP1",
                "mhc_restriction": "HLA-A*02:01",
                "mhc_species": "Homo sapiens",
                "is_ms_observation": True,
                "source_tissue": "Heart",
                "src_cancer": True,
                "src_healthy_tissue": True,
                "src_healthy_reproductive": False,
                "src_healthy_reproductive_female": False,
                "src_healthy_reproductive_male": False,
                "pmid": "444",
            },
        ]
    )


def test_ms_evidence_and_attribution_propagate_sex_reproductive_flags():
    obs = _mock_observations()
    agg = aggregate_ms_hits_by_peptide(obs).set_index("peptide")
    assert bool(agg.loc["OVARYPEP1", "ms_healthy_reproductive_female"]) is True
    assert bool(agg.loc["OVARYPEP1", "ms_healthy_reproductive_male"]) is False
    assert bool(agg.loc["PROSTPEP1", "ms_healthy_reproductive_male"]) is True
    assert bool(agg.loc["PROSTPEP1", "ms_healthy_reproductive_female"]) is False

    cands = pd.DataFrame({"peptide": ["OVARYPEP1", "PROSTPEP1", "TESTISPEP", "HEARTPEP1"]})
    attributed = attribute_ms_evidence(
        cands, obs, requested_alleles=["HLA-A*02:01"], evidence_loaded=True
    )
    summary = summarize_ms_attribution(attributed).set_index("peptide")
    assert bool(summary.loc["OVARYPEP1", "ms_in_healthy_reproductive_female"]) is True
    assert bool(summary.loc["OVARYPEP1", "ms_in_healthy_reproductive_male"]) is False
    assert bool(summary.loc["PROSTPEP1", "ms_in_healthy_reproductive_male"]) is True
    assert bool(summary.loc["TESTISPEP", "ms_in_healthy_reproductive"]) is True
    assert bool(summary.loc["TESTISPEP", "ms_in_healthy_reproductive_female"]) is False
    assert bool(summary.loc["TESTISPEP", "ms_in_healthy_reproductive_male"]) is False


def test_personalized_targets_enforces_same_sex_reproductive_ms_under_all_adult(monkeypatch):
    def fake_score(peptides, alleles, predictor="mhcflurry"):
        return pd.DataFrame(
            [
                {
                    "peptide": pep,
                    "allele": allele,
                    "presentation_percentile": 0.1,
                    "presentation_score": 0.9,
                    "affinity_nm": 20.0,
                }
                for pep in peptides
                for allele in alleles
            ]
        )

    monkeypatch.setattr("tsarina.scoring.score_presentation", fake_score)
    custom = pd.DataFrame(
        {
            "peptide": ["OVARYPEP1", "PROSTPEP1", "TESTISPEP", "HEARTPEP1"],
            "category": ["cta", "cta", "cta", "cta"],
            "source": ["G1", "G2", "G3", "G4"],
        }
    )
    obs = _mock_observations()

    # Female all_adult: drops HEARTPEP1 (somatic) and OVARYPEP1 (same-sex female repro),
    # keeps PROSTPEP1 (opposite-sex male repro) and TESTISPEP (core CTA repro).
    fem_all = tsarina.personalized_targets(
        hla_alleles=["HLA-A*02:01"],
        custom_targets=custom,
        ms_observations=obs,
        patient_sex="female",
        sex_eligibility="all_adult",
        show_progress=False,
    )
    assert set(fem_all["peptide"]) == {"PROSTPEP1", "TESTISPEP"}

    # Male all_adult: drops HEARTPEP1 and PROSTPEP1 (same-sex male repro),
    # keeps OVARYPEP1 (opposite-sex female repro) and TESTISPEP.
    male_all = tsarina.personalized_targets(
        hla_alleles=["HLA-A*02:01"],
        custom_targets=custom,
        ms_observations=obs,
        patient_sex="male",
        sex_eligibility="all_adult",
        show_progress=False,
    )
    assert set(male_all["peptide"]) == {"OVARYPEP1", "TESTISPEP"}

    # Female sex_organ_only: retains OVARYPEP1 (flagged ms_in_healthy_reproductive_female=True),
    # still drops HEARTPEP1 (somatic).
    fem_sex_org = tsarina.personalized_targets(
        hla_alleles=["HLA-A*02:01"],
        custom_targets=custom,
        ms_observations=obs,
        patient_sex="female",
        sex_eligibility="sex_organ_only",
        show_progress=False,
    )
    assert set(fem_sex_org["peptide"]) == {"OVARYPEP1", "PROSTPEP1", "TESTISPEP"}
    ov_row = fem_sex_org[fem_sex_org["peptide"] == "OVARYPEP1"].iloc[0]
    assert bool(ov_row["ms_in_healthy_reproductive_female"]) is True

    table = pmod.format_table(fem_sex_org)
    assert "cancer+repro_f" in table
    assert "cancer+repro_m" in table


def test_cta_healthy_tissue_ms_hits_respects_patient_sex(monkeypatch):
    import tsarina.indexing

    rows = _mock_observations().assign(
        mhc_allele_set="HLA-A*02:01",
        mhc_allele_provenance="exact",
        cell_name="",
    )
    monkeypatch.setattr(tsarina.indexing, "load_ms_evidence", lambda **kw: rows, raising=True)

    # Default (no patient_sex): somatic only
    base = tsarina.cta_healthy_tissue_ms_hits("PRAME")
    assert list(base["peptide"]) == ["HEARTPEP1"]

    # Female all_adult: somatic + female reproductive tract
    fem = tsarina.cta_healthy_tissue_ms_hits("PRAME", patient_sex="female")
    assert set(fem["peptide"]) == {"HEARTPEP1", "OVARYPEP1"}

    # Male all_adult: somatic + male accessory reproductive tract
    male = tsarina.cta_healthy_tissue_ms_hits("PRAME", patient_sex="male")
    assert set(male["peptide"]) == {"HEARTPEP1", "PROSTPEP1"}

    # sex_organ_only: somatic only (same-sex reproductive permitted)
    fem_org = tsarina.cta_healthy_tissue_ms_hits(
        "PRAME", patient_sex="female", sex_eligibility="sex_organ_only"
    )
    assert list(fem_org["peptide"]) == ["HEARTPEP1"]


def test_cli_hits_healthy_tissue_forwards_sex_flags(capsys):
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd")
    p_hits = cli_hits.build_parser(sub)
    args = p_hits.parse_args(
        ["--gene", "PRAME", "--healthy-tissue", "--sex", "female", "--sex-eligibility", "all_adult"]
    )
    empty = pd.DataFrame(columns=tsarina.ms_evidence.CTA_HEALTHY_TISSUE_MS_COLUMNS)
    with patch("tsarina.ms_evidence.cta_healthy_tissue_ms_hits", return_value=empty) as mock_fn:
        cli_hits.handle(args)
    mock_fn.assert_called_once_with(
        "PRAME",
        mhc_class=None,
        mhc_species="Homo sapiens",
        patient_sex="female",
        sex_eligibility="all_adult",
    )
