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

"""Tsarina-side tests for the ``--format refs`` path.

The aggregation itself lives in hitlist 1.12.0+ as
``hitlist.aggregate.aggregate_per_pmhc_with_refs`` and is tested there
(``tests/test_aggregate.py``). These tests only cover the tsarina-side
wiring: CLI format registration + end-to-end dispatch to the hitlist
function with ``ms_pmhc_*`` columns surfaced through the handler.
"""

from __future__ import annotations

import pandas as pd
from hitlist.aggregate import aggregate_per_pmhc_with_refs

from tsarina.cli_hits import _SUPPORTED_FORMATS


def test_refs_is_in_cli_format_choices():
    assert "refs" in _SUPPORTED_FORMATS


def test_hitlist_aggregator_produces_ms_pmhc_prefixed_columns():
    """Contract check against hitlist.aggregate.aggregate_per_pmhc_with_refs:
    the output columns tsarina's --format refs path surfaces should carry
    the ``ms_pmhc_*`` prefix.  If hitlist ever renames these, this test
    catches it before the CLI output silently drifts."""
    hits = pd.DataFrame(
        {
            "peptide": ["QYIAQFTSQF", "QYIAQFTSQF"],
            "mhc_restriction": ["HLA-A*24:02", "HLA-A*24:02"],
            "pmid": [27869121, 38000001],
            "source_tissue": ["Lymph Node", "Skin"],
            "src_cancer": [True, True],
            "src_healthy_tissue": [False, False],
            "is_monoallelic": [False, True],
        }
    )
    out = aggregate_per_pmhc_with_refs(hits)
    assert len(out) == 1
    row = out.iloc[0]
    # Core hitlist-contract columns the tsarina CLI depends on.
    for col in (
        "ms_pmhc_hit_count",
        "ms_pmhc_ref_count",
        "ms_pmhc_pmids",
        "ms_pmhc_tissues",
        "ms_pmhc_in_cancer",
        "ms_pmhc_in_healthy_tissue",
        "ms_pmhc_mono_hit_count",
    ):
        assert col in out.columns, f"hitlist aggregator changed — missing {col}"
    assert row["ms_pmhc_pmids"] == "27869121;38000001"
    assert bool(row["ms_pmhc_in_cancer"]) is True
    assert int(row["ms_pmhc_mono_hit_count"]) == 1
