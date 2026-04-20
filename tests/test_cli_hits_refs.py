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

"""Tests for the ``--format refs`` aggregation in ``tsarina hits``."""

from __future__ import annotations

import pandas as pd

from tsarina.cli_hits import _aggregate_refs


def _scan_fixture() -> pd.DataFrame:
    """Two peptides x one allele, with two references for peptide A and
    mixed cancer / healthy-tissue provenance."""
    return pd.DataFrame(
        {
            "peptide": [
                "QYIAQFTSQF",
                "QYIAQFTSQF",
                "LYVDSLFFL",
            ],
            "mhc_restriction": [
                "HLA-A*24:02",
                "HLA-A*24:02",
                "HLA-A*24:02",
            ],
            "pmid": [27869121, 38000001, 33858848],
            "source_tissue": ["Lymph Node", "Skin", "Thymus"],
            "disease": ["skin melanoma", "skin melanoma", ""],
            "cell_line_name": ["", "", ""],
            "src_cancer": [True, True, False],
            "src_healthy_tissue": [False, False, False],
            "src_healthy_thymus": [False, False, True],
            "is_monoallelic": [False, True, False],
        }
    )


def test_refs_aggregates_per_pmhc_with_counts_and_lists():
    out = _aggregate_refs(_scan_fixture())
    assert set(out.columns) >= {
        "peptide",
        "length",
        "mhc_restriction",
        "hit_count",
        "ref_count",
        "pmids",
        "tissues",
        "diseases",
        "in_cancer",
        "in_healthy_tissue",
        "mono_allelic_hit_count",
    }
    assert len(out) == 2  # two distinct (peptide, allele) rows

    q_row = out[out["peptide"] == "QYIAQFTSQF"].iloc[0]
    assert int(q_row["length"]) == 10
    assert int(q_row["hit_count"]) == 2
    assert int(q_row["ref_count"]) == 2
    # pmids are sorted strings; both PMIDs appear
    assert q_row["pmids"] == "27869121;38000001"
    assert q_row["tissues"] == "Lymph Node;Skin"
    assert q_row["diseases"] == "skin melanoma"
    assert bool(q_row["in_cancer"]) is True
    assert bool(q_row["in_healthy_tissue"]) is False
    assert int(q_row["mono_allelic_hit_count"]) == 1

    l_row = out[out["peptide"] == "LYVDSLFFL"].iloc[0]
    assert int(l_row["length"]) == 9
    assert int(l_row["ref_count"]) == 1
    assert l_row["pmids"] == "33858848"
    assert l_row["tissues"] == "Thymus"
    assert l_row["diseases"] == ""
    assert bool(l_row["in_cancer"]) is False


def test_refs_empty_input_returns_canonical_columns():
    out = _aggregate_refs(pd.DataFrame(columns=["peptide", "mhc_restriction"]))
    assert out.empty
    assert list(out.columns) == [
        "peptide",
        "length",
        "mhc_restriction",
        "hit_count",
        "ref_count",
        "pmids",
        "tissues",
        "diseases",
        "cell_lines",
        "in_cancer",
        "in_healthy_tissue",
        "mono_allelic_hit_count",
    ]


def test_refs_missing_optional_columns_skips_them():
    """Raw scan vs cached observations produce slightly different column
    sets. The aggregator must silently skip any optional columns absent
    from the input rather than crashing."""
    skinny = pd.DataFrame(
        {
            "peptide": ["ABCDEFGHI", "ABCDEFGHI"],
            "mhc_restriction": ["HLA-A*02:01", "HLA-A*02:01"],
            # no pmid / source_tissue / src_cancer / is_monoallelic columns
        }
    )
    out = _aggregate_refs(skinny)
    assert len(out) == 1
    assert int(out.iloc[0]["hit_count"]) == 2
    # Optional columns that weren't supplied are absent from the output,
    # not None / NaN artifacts.
    assert "pmids" not in out.columns
    assert "tissues" not in out.columns
    assert "in_cancer" not in out.columns


def test_refs_is_in_cli_format_choices():
    from tsarina.cli_hits import _SUPPORTED_FORMATS

    assert "refs" in _SUPPORTED_FORMATS
