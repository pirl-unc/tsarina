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

"""Tests for native mouse CTA & mouse sex-restricted workflows (oncoref>=1.49.0)."""

from __future__ import annotations

import json
import os
import subprocess
import sys

import pandas as pd
import pytest

import tsarina


def _run_cli(*args, cache_dir):
    env = {**os.environ, "HITLIST_DATA_DIR": str(cache_dir)}
    return subprocess.run(
        [sys.executable, "-m", "tsarina.cli", *args],
        capture_output=True,
        text=True,
        env=env,
    )


def test_mouse_sex_restricted_antigens_and_gene_helpers():
    all_mouse = tsarina.CTA_sex_restricted_antigens(species="mouse")
    assert len(all_mouse) >= 100
    assert {"entity_id", "label", "gene_ids", "category", "sex_restriction_class"} <= set(
        all_mouse.columns
    )

    fem_ids = tsarina.CTA_sex_eligible_gene_ids("female", species="mouse")
    male_ids = tsarina.CTA_sex_eligible_gene_ids("male", species="mouse")
    fem_names = tsarina.CTA_sex_eligible_gene_names("female", species="mouse")
    male_names = tsarina.CTA_sex_eligible_gene_names("male", species="mouse")

    assert fem_ids and all(g.startswith("ENSMUSG") for g in fem_ids)
    assert male_ids and all(g.startswith("ENSMUSG") for g in male_ids)
    assert fem_names and male_names
    assert fem_names != male_names

    # Unevaluated species (dog) raises a clear ValueError
    with pytest.raises(ValueError, match="not yet available"):
        tsarina.CTA_sex_restricted_antigens(sex="female", species="dog")


def test_cli_catalog_ctas_mouse_sex_and_dog_rejection(tmp_path):
    r_mouse = _run_cli(
        "catalog",
        "ctas",
        "--species",
        "mouse",
        "--sex",
        "female",
        "--top",
        "10",
        "--json",
        cache_dir=tmp_path,
    )
    assert r_mouse.returncode == 0, r_mouse.stderr
    payload = json.loads(r_mouse.stdout)
    assert payload["species"] == "mouse"
    assert payload["sex"] == "female"
    assert payload["total"] >= 1
    assert all("sex_restriction_class" in row for row in payload["ctas"])

    r_dog = _run_cli("catalog", "ctas", "--species", "dog", "--sex", "female", cache_dir=tmp_path)
    assert r_dog.returncode == 1
    assert "not yet available" in r_dog.stderr.lower()


def test_personalized_targets_native_mouse_cta_and_sex_eligibility(monkeypatch):
    def fake_gene_peptides(gene_ids, ensembl_release=112, lengths=(8, 9, 10, 11), **kwargs):
        assert ensembl_release == 115
        assert kwargs.get("species") == "mus_musculus"
        rows = []
        for gid in gene_ids:
            rows.append(
                {
                    "gene_name": "Plac1" if gid == "ENSMUSG00000061082" else gid,
                    "gene_id": gid,
                    "transcript_id": f"T_{gid}",
                    "peptide": "SIINFEKL" if gid == "ENSMUSG00000061082" else "SELFPEPT",
                    "length": 8,
                    "start": 1,
                    "end": 8,
                    "n_flank": "",
                    "c_flank": "",
                }
            )
            rows.append(
                {
                    "gene_name": "Plac1" if gid == "ENSMUSG00000061082" else gid,
                    "gene_id": gid,
                    "transcript_id": f"T_{gid}",
                    "peptide": "SELFPEPT",
                    "length": 8,
                    "start": 9,
                    "end": 16,
                    "n_flank": "",
                    "c_flank": "",
                }
            )
        return pd.DataFrame(rows)

    def fake_score(peptides, alleles, predictor="mhcflurry"):
        return pd.DataFrame(
            [
                {
                    "peptide": pep,
                    "allele": allele,
                    "presentation_percentile": 0.2,
                    "presentation_score": 0.9,
                    "affinity_nm": 18.0,
                }
                for pep in peptides
                for allele in alleles
            ]
        )

    monkeypatch.setattr("tsarina.peptides.gene_peptides", fake_gene_peptides)
    monkeypatch.setattr("tsarina.scoring.score_presentation", fake_score)

    # Plac1 is an evaluated core mouse CTA (ENSMUSG00000061082).
    # SELFPEPT matches host_proteins and must be screened out!
    out = tsarina.personalized_targets(
        hla_alleles=["H-2-Kb"],
        cta_expression={"Plac1": 45.0},
        host_species="mouse",
        lengths=(8,),
        host_proteins={"MOUSE_SELF": "AAASELFPEPTAAA"},
        skip_ms_evidence=True,
        show_progress=False,
    )
    assert len(out) == 1
    assert out.iloc[0]["peptide"] == "SIINFEKL"
    assert out.iloc[0]["source"] == "Plac1"
    assert out.iloc[0]["category"] == "cta"
    assert out.attrs["ensembl_release"] == 115
    assert out.attrs["species_context"]["host_taxon"] == 10090
