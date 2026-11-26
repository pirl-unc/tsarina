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

"""Tests for OncoRef sex-restricted CTA eligibility in API, personalize, and catalog CLI."""

from __future__ import annotations

import json
import os
import subprocess
import sys

import pandas as pd
import pytest

import tsarina
from tsarina import (
    CTA_sex_eligible_gene_ids,
    CTA_sex_eligible_gene_names,
    CTA_sex_restricted_antigens,
    personalized_targets,
)


def test_sex_restricted_antigens_catalog_and_gene_sets():
    all_rows = CTA_sex_restricted_antigens()
    assert len(all_rows) >= 400
    assert {"entity_id", "label", "gene_ids", "category", "sex_restriction_class"} <= set(
        all_rows.columns
    )

    female_all = CTA_sex_restricted_antigens(sex="female", eligibility="all_adult")
    male_all = CTA_sex_restricted_antigens(sex="male", eligibility="all_adult")
    assert female_all["female_all_adult_eligible"].all()
    assert male_all["male_all_adult_eligible"].all()
    assert len(female_all) > len(male_all) > 0

    female_cta_only = CTA_sex_restricted_antigens(
        sex="female", eligibility="all_adult", include_reproductive_only=False
    )
    assert female_cta_only["is_cta"].all()
    assert len(female_cta_only) < len(female_all)

    female_names = CTA_sex_eligible_gene_names("female", eligibility="all_adult")
    male_names = CTA_sex_eligible_gene_names("male", eligibility="all_adult")
    assert "PRM2" in female_names and "PRM2" not in male_names
    assert "OR1N1" in male_names and "OR1N1" not in female_names
    assert "DUX4" in female_names and "DUX4" in male_names
    assert len(CTA_sex_eligible_gene_ids("female")) >= len(female_all)

    with pytest.raises(ValueError, match="Unsupported patient sex"):
        CTA_sex_restricted_antigens(sex="unknown")
    with pytest.raises(ValueError, match="Unsupported sex eligibility"):
        CTA_sex_restricted_antigens(sex="female", eligibility="invalid")

    for name in (
        "CTA_sex_restricted_antigens",
        "CTA_sex_eligible_gene_ids",
        "CTA_sex_eligible_gene_names",
    ):
        assert name in tsarina.__all__


def test_personalized_targets_admits_sex_safe_reproductive_and_filters_ineligible(monkeypatch):
    # Pick one male-restricted Reproductive Antigen only (PRM2 -> female all_adult eligible)
    # and one shared male & female CTA (ineligible under strict all_adult).
    sex_df = CTA_sex_restricted_antigens()
    shared_cta = sex_df[sex_df["is_cta"] & (~sex_df["female_all_adult_eligible"])].iloc[0]["label"]

    # Stub peptide enumeration & background screening so the test is fast and offline.
    monkeypatch.setattr(
        "tsarina.peptides.gene_peptides",
        lambda gene_ids, **kw: pd.DataFrame(
            [
                {
                    "peptide": "PRM2PEPT1",
                    "length": 9,
                    "gene_name": "PRM2",
                    "gene_id": next(iter(gene_ids)),
                },
                {
                    "peptide": "BGOVERLAP",
                    "length": 9,
                    "gene_name": "PRM2",
                    "gene_id": next(iter(gene_ids)),
                },
            ]
        ),
    )
    monkeypatch.setattr(
        "tsarina.partition.CTA_partition_gene_ids",
        lambda **kw: type("P", (), {"non_cta": {"BG1", "BG2"}})(),
    )
    monkeypatch.setattr(
        "tsarina.peptides._reference_overlapping_peptides",
        lambda candidate_peptides, **kw: candidate_peptides & {"BGOVERLAP"},
    )

    with pytest.warns(UserWarning, match=f"{shared_cta}.*not female all_adult eligible"):
        res_f = personalized_targets(
            hla_alleles=["HLA-A*02:01"],
            cta_expression={"PRM2": 50.0, shared_cta: 50.0},
            patient_sex="female",
            sex_eligibility="all_adult",
            score_presentation=False,
            skip_ms_evidence=True,
            drop_weak_tier=False,
            show_progress=False,
        )
    assert list(res_f["peptide"]) == ["PRM2PEPT1"]  # BGOVERLAP excluded by self-background!
    assert list(res_f["source"]) == ["PRM2"]
    assert res_f.attrs["patient_sex"] == "female"
    assert res_f.attrs["sex_eligibility"] == "all_adult"

    # In a male patient under all_adult, PRM2 (Male-restricted) is NOT admitted.
    with pytest.warns(UserWarning, match="PRM2"):
        res_m = personalized_targets(
            hla_alleles=["HLA-A*02:01"],
            cta_expression={"PRM2": 50.0},
            patient_sex="male",
            sex_eligibility="all_adult",
            score_presentation=False,
            skip_ms_evidence=True,
            drop_weak_tier=False,
            show_progress=False,
        )
    assert res_m.empty
    assert res_m.attrs["patient_sex"] == "male"


def test_catalog_ctas_sex_cli_filter(tmp_path):
    env = {**os.environ, "HITLIST_DATA_DIR": str(tmp_path)}
    r = subprocess.run(
        [
            sys.executable,
            "-m",
            "tsarina.cli",
            "catalog",
            "ctas",
            "--sex",
            "female",
            "--top",
            "5",
            "--json",
        ],
        capture_output=True,
        text=True,
        env=env,
    )
    assert r.returncode == 0, r.stderr
    data = json.loads(r.stdout)
    assert data["sex"] == "female"
    assert data["sex_eligibility"] == "all_adult"
    assert data["shown"] == 5
    assert all("sex_restriction_class" in row for row in data["ctas"])
    assert all("sex_restriction_code" in row for row in data["ctas"])

    dog_ok = subprocess.run(
        [
            sys.executable,
            "-m",
            "tsarina.cli",
            "catalog",
            "ctas",
            "--species",
            "dog",
            "--sex",
            "female",
            "--top",
            "5",
            "--json",
        ],
        capture_output=True,
        text=True,
        env=env,
    )
    assert dog_ok.returncode == 0, dog_ok.stderr
    dog_data = json.loads(dog_ok.stdout)
    assert dog_data["species"] == "dog"
    assert dog_data["sex"] == "female"
    assert dog_data["total"] >= 1
    assert dog_data["shown"] == 5
    assert all("sex_restriction_code" in row for row in dog_data["ctas"])
