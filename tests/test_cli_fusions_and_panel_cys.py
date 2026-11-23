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

"""Tests for personalize --fusion/--fusions CLI and panel Cys MS-detectability caveats."""

from __future__ import annotations

import argparse
import json
from unittest.mock import patch

import pandas as pd
import pytest

from tsarina import cli_personalize, cli_spanning, spanning


def test_parse_inline_fusion_specs_and_file_loaders(tmp_path):
    inline = cli_personalize._parse_fusion_specs(
        ["EWSR1-FLI1:sygqqssy|npsydsvr", "AAAABBBB|CCCCDDDD"]
    )
    assert inline == [
        {
            "n_terminal_seq": "SYGQQSSY",
            "c_terminal_seq": "NPSYDSVR",
            "label": "EWSR1-FLI1",
            "gene_5prime": "EWSR1",
            "gene_3prime": "FLI1",
        },
        {
            "n_terminal_seq": "AAAABBBB",
            "c_terminal_seq": "CCCCDDDD",
        },
    ]

    with pytest.raises(argparse.ArgumentTypeError, match="must contain '\\|'"):
        cli_personalize._parse_fusion_specs(["BADSPEC"])

    csv_path = tmp_path / "fusions.csv"
    csv_path.write_text(
        "label,gene_5prime,gene_3prime,n_terminal_seq,c_terminal_seq\n"
        "SS18-SSX1,SS18,SSX1,QQGPPPQQ,PQMTFGGR\n"
    )
    loaded_csv = cli_personalize._load_fusions_file(str(csv_path))
    assert len(loaded_csv) == 1
    assert loaded_csv[0]["label"] == "SS18-SSX1"

    json_path = tmp_path / "fusions.json"
    json_path.write_text(
        json.dumps([{"label": "BCR-ABL1", "sequence": "KQSSKALQRPVASDF", "junction_pos": 7}])
    )
    loaded_json = cli_personalize._load_fusions_file(str(json_path))
    assert loaded_json == [{"label": "BCR-ABL1", "sequence": "KQSSKALQRPVASDF", "junction_pos": 7}]


def test_cli_personalize_forwards_inline_and_file_fusions(tmp_path, capsys):
    json_path = tmp_path / "f.json"
    json_path.write_text(
        json.dumps({"label": "F2", "n_terminal_seq": "SLLMW", "c_terminal_seq": "ITQC"})
    )

    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd")
    p_pers = cli_personalize.build_parser(sub)
    args = p_pers.parse_args(
        [
            "--hla",
            "HLA-A*02:01",
            "--fusion",
            "F1:SYGQQ|NPSYD",
            "--fusions",
            str(json_path),
            "--no-require-host-exclusive-fusions",
            "--skip-ms-evidence",
            "--no-score",
        ]
    )

    empty = pd.DataFrame(columns=["peptide", "category", "source"])
    with patch("tsarina.personalize.personalized_targets", return_value=empty) as mock_fn:
        cli_personalize.handle(args)

    mock_fn.assert_called_once()
    kwargs = mock_fn.call_args.kwargs
    assert len(kwargs["fusions"]) == 2
    assert kwargs["fusions"][0]["label"] == "F1"
    assert kwargs["fusions"][1]["label"] == "F2"
    assert kwargs["require_host_exclusive_fusions"] is False


def test_panel_long_output_and_table_surface_cys_ms_detectability_caveat():
    selected = pd.DataFrame(
        [
            {
                "cta": "NY-ESO-1",
                "cta_members": "CTAG1A;CTAG1B",
                "allele": "HLA-A*02:01",
                "peptide": "SLLMWITQC",
                "length": 9,
                "evidence_tier": "predicted_only",
                "ms_hit_count": 0,
                "ms_source_count": 0,
                "ms_alleles": "",
                "ms_pmids": "",
                "ms_samples": "",
                "ms_in_cancer": False,
                "ms_in_healthy_tissue": False,
                "ms_detectability_note": "",
                "presentation_percentile": 0.15,
                "presentation_score": 0.91,
                "affinity_nm": 28.0,
                "affinity_percentile": 0.20,
                "netmhcpan_affinity_nm": pd.NA,
                "netmhcpan_affinity_percentile": pd.NA,
                "peptide_rank_in_cell": 1,
            },
            {
                "cta": "PRAME",
                "cta_members": "PRAME",
                "allele": "HLA-A*02:01",
                "peptide": "SLLQHLIGL",
                "length": 9,
                "evidence_tier": "monoallelic_ms",
                "ms_hit_count": 32,
                "ms_source_count": 10,
                "ms_alleles": "HLA-A*02:01",
                "ms_pmids": "111;222",
                "ms_samples": "s1",
                "ms_in_cancer": True,
                "ms_in_healthy_tissue": False,
                "ms_detectability_note": "",
                "presentation_percentile": 0.01,
                "presentation_score": 0.99,
                "affinity_nm": 15.9,
                "affinity_percentile": 0.05,
                "netmhcpan_affinity_nm": pd.NA,
                "netmhcpan_affinity_percentile": pd.NA,
                "peptide_rank_in_cell": 1,
            },
        ]
    )
    long_df = spanning._to_long(selected, ["NY-ESO-1", "PRAME"], ["HLA-A*02:01"])
    assert "ms_detectability_note" in long_df.columns
    ny_note = long_df.loc[long_df["peptide"] == "SLLMWITQC", "ms_detectability_note"].iloc[0]
    prame_note = long_df.loc[long_df["peptide"] == "SLLQHLIGL", "ms_detectability_note"].iloc[0]
    assert "C-terminal Cys" in ny_note
    assert prame_note == ""

    long_df.attrs["cta_order"] = ["NY-ESO-1", "PRAME"]
    long_df.attrs["allele_order"] = ["HLA-A*02:01"]
    rendered = cli_spanning.format_panel_table(long_df)
    assert "1 zero-MS panel peptide(s) contain Cys" in rendered
