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

import tsarina


def _run_cli(*args, cache_dir):
    env = {**os.environ, "HITLIST_DATA_DIR": str(cache_dir)}
    return subprocess.run(
        [sys.executable, "-m", "tsarina.cli", *args],
        capture_output=True,
        text=True,
        env=env,
    )


def test_mouse_and_dog_sex_restricted_antigens_and_gene_helpers():
    for sp, prefix in (("mouse", "ENSMUSG"), ("dog", "ENSCAFG")):
        all_rows = tsarina.CTA_sex_restricted_antigens(species=sp)
        assert len(all_rows) >= 100
        assert {
            "entity_id",
            "label",
            "gene_ids",
            "category",
            "sex_restriction_class",
            "sex_restriction_code",
        } <= set(all_rows.columns)

        fem_ids = tsarina.CTA_sex_eligible_gene_ids("female", species=sp)
        male_ids = tsarina.CTA_sex_eligible_gene_ids("male", species=sp)
        fem_names = tsarina.CTA_sex_eligible_gene_names("female", species=sp)
        male_names = tsarina.CTA_sex_eligible_gene_names("male", species=sp)

        assert fem_ids and all(g.startswith(prefix) for g in fem_ids)
        assert male_ids and all(g.startswith(prefix) for g in male_ids)
        assert fem_names and male_names
        assert fem_names != male_names


def test_cli_catalog_ctas_mouse_and_dog_sex(tmp_path):
    for sp in ("mouse", "dog"):
        r = _run_cli(
            "catalog",
            "ctas",
            "--species",
            sp,
            "--sex",
            "female",
            "--top",
            "10",
            "--json",
            cache_dir=tmp_path,
        )
        assert r.returncode == 0, r.stderr
        payload = json.loads(r.stdout)
        assert payload["species"] == sp
        assert payload["sex"] == "female"
        assert payload["total"] >= 1
        assert all("sex_restriction_class" in row for row in payload["ctas"])
        assert all("sex_restriction_code" in row for row in payload["ctas"])


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
