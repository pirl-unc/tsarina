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

"""Tests for curated cancer fusions, orfdata viral subsets, and Cys MS-detectability caveats."""

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


def test_cancer_fusions_filtering_matches_oncoref():
    all_df = tsarina.cancer_fusions()
    assert len(all_df) >= 100
    assert {"cancer_code", "gene_5prime", "gene_3prime", "is_defining"} <= set(all_df.columns)

    ewing = tsarina.cancer_fusions(cancer_code="sarc_ews", defining_only=True)
    assert not ewing.empty
    assert (ewing["cancer_code"] == "SARC_EWS").all()
    assert ewing["is_defining"].all()
    assert "FLI1" in set(ewing["gene_3prime"])

    alk = tsarina.cancer_fusions(gene="alk", rnaseq_detectable_only=True)
    assert not alk.empty
    assert ((alk["gene_5prime"] == "ALK") | (alk["gene_3prime"] == "ALK")).all()
    assert (alk["rnaseq_detectable"] == "yes").all()

    ets = tsarina.cancer_fusions(fusion_family="ewsr1-ets")
    assert not ets.empty
    assert (ets["fusion_family"] == "EWSR1-ETS").all()


def test_fusion_junction_peptides_spans_breakpoint_and_screens_host():
    # Split (n_terminal_seq, c_terminal_seq) form
    split_df = tsarina.fusion_junction_peptides(
        {
            "gene_5prime": "EWSR1",
            "gene_3prime": "FLI1",
            "n_terminal_seq": "SYGQQSSY",
            "c_terminal_seq": "NPSYDSVR",
        },
        lengths=(9,),
        host_proteins={"SELF1": "MMMYGQQSSYNPMMM"},
    )
    # 8 + 8 = 16 aa, junction at index 8 -> 8 spanning 9-mers before self-screening.
    # "YGQQSSYNP" matches SELF1 and must be screened out -> 7 remain.
    assert len(split_df) == 7
    assert "YGQQSSYNP" not in set(split_df["peptide"])
    assert (split_df["label"] == "EWSR1-FLI1").all()
    assert (split_df["n_terminal_residues"] >= 1).all()
    assert (split_df["c_terminal_residues"] >= 1).all()
    assert (~split_df["in_reference_proteome"]).all()

    # Joined (sequence, junction_pos) form with require_host_exclusive=False keeps & flags overlap
    raw_df = tsarina.fusion_junction_peptides(
        [
            {
                "label": "EWSR1-FLI1",
                "sequence": "SYGQQSSYNPSYDSVR",
                "junction_pos": 8,
            }
        ],
        lengths=(9,),
        require_host_exclusive=False,
        host_proteins={"SELF1": "MMMYGQQSSYNPMMM"},
    )
    assert len(raw_df) == 8
    flagged = raw_df.loc[raw_df["peptide"] == "YGQQSSYNP", "in_reference_proteome"]
    assert flagged.tolist() == [True]

    with pytest.raises(ValueError, match="junction_pos"):
        tsarina.fusion_junction_peptides({"sequence": "ACDEFGHIK", "junction_pos": 0})


def test_personalized_targets_accepts_fusions_and_annotates_cys_detectability(monkeypatch):
    from tsarina import personalize as pmod

    def fake_score(peptides, alleles, predictor="mhcflurry"):
        rows = []
        for pep in peptides:
            for allele in alleles:
                rows.append(
                    {
                        "peptide": pep,
                        "allele": allele,
                        "presentation_percentile": 0.15,
                        "presentation_score": 0.9,
                        "affinity_nm": 25.0,
                    }
                )
        return pd.DataFrame(rows)

    monkeypatch.setattr("tsarina.scoring.score_presentation", fake_score)

    out = tsarina.personalized_targets(
        hla_alleles=["HLA-A*02:01"],
        fusions=[
            {
                "label": "CUSTOM-FUS",
                "gene_5prime": "GENE5",
                "gene_3prime": "GENE3",
                "n_terminal_seq": "SLLMW",
                "c_terminal_seq": "ITQC",
            }
        ],
        lengths=(9,),
        host_proteins={"HOST": "MMMMMMMMMMMM"},
        skip_ms_evidence=True,
        show_progress=False,
    )
    assert len(out) == 1
    row = out.iloc[0]
    assert row["peptide"] == "SLLMWITQC"
    assert row["category"] == "fusion"
    assert row["source"] == "CUSTOM-FUS"
    assert row["tier"] == 1
    assert "C-terminal Cys" in row["ms_detectability_note"]

    rendered = pmod.format_table(out)
    assert "contain Cys" in rendered


def test_annotate_ms_detectability_distinguishes_cterm_internal_and_observed():
    df = pd.DataFrame(
        {
            "peptide": ["SLLMWITQC", "ACDEFGHIK", "SLLQHLIGL", "ACDEFGHIK"],
            "ms_hit_count": [0, 0, 0, 5],
        }
    )
    ann = tsarina.annotate_ms_detectability(df)
    assert ann["has_cysteine"].tolist() == [True, True, False, True]
    assert ann["c_terminal_cysteine"].tolist() == [True, False, False, False]
    assert "C-terminal Cys" in ann.loc[0, "ms_detectability_note"]
    assert "Contains Cys" in ann.loc[1, "ms_detectability_note"]
    assert ann.loc[2, "ms_detectability_note"] == ""
    # Already observed by MS -> no zero-MS detectability caveat
    assert ann.loc[3, "ms_detectability_note"] == ""


def test_orfdata_viral_subsets_and_bundled_refseq_orfs():
    pytest.importorskip("orfdata")

    subsets = tsarina.orfdata_viral_subsets()
    names = {s["name"] for s in subsets}
    assert {"cancer-associated", "cancer-causing", "cancer-integrated-driver"} <= names
    by_name = {s["name"]: s for s in subsets}
    assert by_name["cancer-causing"]["virus_count"] >= 15

    # Bundled RefSeq HPV16 ORFs work completely offline
    hpv16 = tsarina.orfdata_viral_peptides(
        reference_id="NC_001526.4",
        lengths=(9,),
        installed=False,
        host_exclusive=True,
        host_proteins={"HUMAN_STUB": "MMMMMMMMMMMM"},
    )
    assert not hpv16.empty
    assert (hpv16["virus"] == "NC_001526.4").all()
    assert {"E6", "E7"} <= set(hpv16["gene_name"])


def test_cli_catalog_fusions_and_orfdata_viruses(tmp_path):
    ov = _run_cli("catalog", "--json", cache_dir=tmp_path)
    assert ov.returncode == 0, ov.stderr
    ov_data = json.loads(ov.stdout)
    assert ov_data["fusions_count"] >= 100
    assert ov_data["defining_fusions_count"] >= 50

    fus = _run_cli(
        "catalog",
        "fusions",
        "--cancer",
        "SARC_EWS",
        "--defining-only",
        "--json",
        cache_dir=tmp_path,
    )
    assert fus.returncode == 0, fus.stderr
    fus_data = json.loads(fus.stdout)
    assert fus_data["total"] >= 1
    assert any(r["fusion"] == "EWSR1-FLI1" for r in fus_data["fusions"])

    fus_text = _run_cli("catalog", "fusions", "--gene", "ALK", "--top", "5", cache_dir=tmp_path)
    assert fus_text.returncode == 0, fus_text.stderr
    assert "ALK" in fus_text.stdout

    # Default viruses --json stays backward-compatible list
    v_default = _run_cli("catalog", "viruses", "--json", cache_dir=tmp_path)
    assert v_default.returncode == 0, v_default.stderr
    assert isinstance(json.loads(v_default.stdout), list)

    pytest.importorskip("orfdata")
    v_orf = _run_cli(
        "catalog",
        "viruses",
        "--source",
        "orfdata",
        "--subset",
        "cancer-causing",
        "--json",
        cache_dir=tmp_path,
    )
    assert v_orf.returncode == 0, v_orf.stderr
    orf_payload = json.loads(v_orf.stdout)
    assert len(orf_payload) == 1
    assert orf_payload[0]["name"] == "cancer-causing"

    v_all = _run_cli("catalog", "viruses", "--source", "all", "--json", cache_dir=tmp_path)
    assert v_all.returncode == 0, v_all.stderr
    all_payload = json.loads(v_all.stdout)
    assert "builtin" in all_payload and "orfdata_subsets" in all_payload


def test_top_level_exports_for_1_47_0():
    for name in (
        "cancer_fusions",
        "fusion_junction_peptides",
        "orfdata_viral_subsets",
        "orfdata_viral_peptides",
        "annotate_ms_detectability",
    ):
        assert name in tsarina.__all__
        assert hasattr(tsarina, name)
