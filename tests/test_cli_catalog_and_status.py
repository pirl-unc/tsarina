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

"""Tests for `tsarina status`, `tsarina catalog`, CLI flag symmetry, and top-level exports."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

import tsarina


def _run(*args, cache_dir):
    env = {**os.environ, "HITLIST_DATA_DIR": str(cache_dir)}
    return subprocess.run(
        [sys.executable, "-m", "tsarina.cli", *args],
        capture_output=True,
        text=True,
        env=env,
    )


def test_status_text_and_json_are_read_only(tmp_path):
    r_text = _run("status", cache_dir=tmp_path)
    assert r_text.returncode == 0, r_text.stderr
    assert f"tsarina {tsarina.__version__}" in r_text.stdout
    assert "OncoRef:" in r_text.stdout
    assert "Hitlist:" in r_text.stdout
    assert "MHCflurry:" in r_text.stdout
    assert list(tmp_path.iterdir()) == []

    r_json = _run("status", "--json", cache_dir=tmp_path)
    assert r_json.returncode == 0, r_json.stderr
    payload = json.loads(r_json.stdout)
    assert payload["tsarina_version"] == tsarina.__version__
    assert payload["oncoref"]["human_cta_count"] > 0
    assert payload["hitlist"]["observations_built"] is False
    assert list(tmp_path.iterdir()) == []


def test_catalog_overview_and_subcommands(tmp_path):
    ov = _run("catalog", "--json", cache_dir=tmp_path)
    assert ov.returncode == 0, ov.stderr
    data = json.loads(ov.stdout)
    assert data["ctas"]["human_cta_count"] > 0
    assert "global54_abc" in data["panels"]
    assert data["mutations_count"] == len(tsarina.HOTSPOT_MUTATIONS)
    assert "hpv16" in data["viruses"]

    ctas = _run("catalog", "ctas", "--tier", "core", "--top", "5", "--json", cache_dir=tmp_path)
    assert ctas.returncode == 0, ctas.stderr
    cta_data = json.loads(ctas.stdout)
    assert cta_data["shown"] == 5
    assert all(row["tier"] == "core" for row in cta_data["ctas"])

    mouse = _run("catalog", "ctas", "--species", "mouse", "--json", cache_dir=tmp_path)
    assert mouse.returncode == 0, mouse.stderr
    mouse_data = json.loads(mouse.stdout)
    assert mouse_data["total"] >= 1
    assert mouse_data["coverage"]["selection_status"] == "pending_reconciled_panels"

    panels = _run("catalog", "panels", cache_dir=tmp_path)
    assert panels.returncode == 0, panels.stderr
    assert "global54_abc" in panels.stdout

    panel_detail = _run("catalog", "panels", "--panel", "iedb27_ab", "--json", cache_dir=tmp_path)
    assert panel_detail.returncode == 0, panel_detail.stderr
    pd_json = json.loads(panel_detail.stdout)
    assert pd_json["panel"] == "iedb27_ab"
    assert len(pd_json["alleles"]) == 27

    bad_panel = _run("catalog", "panels", "--panel", "not_a_panel", cache_dir=tmp_path)
    assert bad_panel.returncode == 1
    assert "unknown panel" in bad_panel.stderr.lower()

    muts = _run("catalog", "mutations", "--json", cache_dir=tmp_path)
    assert muts.returncode == 0, muts.stderr
    assert any(m["label"] == "KRAS G12D" for m in json.loads(muts.stdout))

    virs = _run("catalog", "viruses", "--json", cache_dir=tmp_path)
    assert virs.returncode == 0, virs.stderr
    assert any(v["key"] == "hpv16" for v in json.loads(virs.stdout))


def test_cli_allele_flag_symmetry_and_grouped_panel_help():
    from tsarina import cli_hits, cli_personalize, cli_spanning

    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd")
    p_pers = cli_personalize.build_parser(sub)
    p_hits = cli_hits.build_parser(sub)
    p_panel = cli_spanning.build_parser(sub)

    assert p_pers.parse_args(["--alleles", "HLA-A02:01"]).hla == ["HLA-A02:01"]
    assert p_hits.parse_args(["--gene", "PRAME", "--hla", "HLA-A*02:01", "--top", "5"]).allele == [
        "HLA-A*02:01"
    ]
    assert p_hits.parse_args(["--gene", "PRAME", "--alleles", "HLA-A*02:01"]).allele == [
        "HLA-A*02:01"
    ]
    assert p_panel.parse_args(["--hla", "HLA-A*02:01"]).alleles == ["HLA-A*02:01"]

    panel_help = p_panel.format_help()
    assert "CTA selection and ranking:" in panel_help
    assert "HLA panel and peptide enumeration:" in panel_help
    assert "MS evidence tiers and presentation scoring:" in panel_help


def test_top_level_headline_exports_are_available_and_in_all():
    for name in (
        "spanning_pmhc_set",
        "design_vaccine",
        "VaccineConfig",
        "CTA_selection_coverage",
    ):
        assert name in tsarina.__all__
        assert hasattr(tsarina, name)
