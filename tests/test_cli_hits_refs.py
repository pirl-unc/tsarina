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


# ── Handler merge-path coverage ────────────────────────────────────────


def test_refs_handler_merges_gene_ident_cols_from_cached_path(tmp_path):
    """When the cached-observations fast path returns rows with
    ``gene_names`` / ``gene_ids`` / ``protein_ids`` columns, the refs
    branch of cli_hits.handle() must join those identity columns onto
    the aggregator output so output shape stays consistent with the
    raw-CSV override path."""
    import argparse
    from unittest.mock import patch

    import pandas as pd

    from tsarina import cli_hits

    hits_frame = pd.DataFrame(
        {
            "peptide": ["QYIAQFTSQF", "LYVDSLFFL"],
            "mhc_restriction": ["HLA-A*24:02", "HLA-A*24:02"],
            "gene_names": ["PRAME", "PRAME"],
            "gene_ids": ["ENSG00000185686", "ENSG00000185686"],
            "protein_ids": ["ENST00000405655", "ENST00000405655"],
            "pmid": [27869121, 33858848],
            "source_tissue": ["Lymph Node", "Thymus"],
            "disease": ["skin melanoma", ""],
            "cell_line_name": ["", ""],
            "src_cancer": [True, False],
            "src_healthy_tissue": [False, False],
            "is_monoallelic": [False, False],
        }
    )

    output_csv = tmp_path / "refs.csv"
    args = argparse.Namespace(
        gene="PRAME",
        uniprot=None,
        allele=[],
        serotype=[],
        species="Homo sapiens",
        mhc_class="I",
        min_resolution=None,
        lengths=(8, 9, 10, 11),
        ensembl_release=112,
        include_binding_assays=False,
        mono_allelic_only=False,
        format="refs",
        predict=False,
        predictor="mhcflurry",
        iedb_path=None,
        cedar_path=None,
        skip_ms_evidence=False,
        output=str(output_csv),
    )
    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_observations", return_value=hits_frame),
    ):
        cli_hits.handle(args)

    out = pd.read_csv(output_csv)
    # Gene-ident columns from the cached path must be merged onto the refs output.
    for col in ("gene_names", "gene_ids", "protein_ids"):
        assert col in out.columns, f"gene_ident_cols merge dropped {col}"
    # Aggregator output shape should also be intact alongside the idents.
    assert "ms_pmhc_hit_count" in out.columns
    assert set(out["peptide"]) == {"QYIAQFTSQF", "LYVDSLFFL"}
    # Each peptide's gene identity should survive the merge without duplication.
    assert all(out["gene_names"] == "PRAME")


def test_refs_handler_passes_through_when_no_gene_ident_cols(tmp_path):
    """When the scan result lacks gene_names / gene_ids / protein_ids
    (e.g. raw-scan override path without the enumeration frame), the
    refs branch should still produce the aggregator output rather than
    erroring or dropping rows."""
    import argparse
    from unittest.mock import patch

    import pandas as pd

    from tsarina import cli_hits

    hits_frame = pd.DataFrame(
        {
            "peptide": ["QYIAQFTSQF"],
            "mhc_restriction": ["HLA-A*24:02"],
            "pmid": [27869121],
            "source_tissue": ["Lymph Node"],
            "disease": ["skin melanoma"],
            "cell_line_name": [""],
            "src_cancer": [True],
            "src_healthy_tissue": [False],
            "is_monoallelic": [False],
        }
    )

    output_csv = tmp_path / "refs.csv"
    args = argparse.Namespace(
        gene="PRAME",
        uniprot=None,
        allele=[],
        serotype=[],
        species="Homo sapiens",
        mhc_class="I",
        min_resolution=None,
        lengths=(8, 9, 10, 11),
        ensembl_release=112,
        include_binding_assays=False,
        mono_allelic_only=False,
        format="refs",
        predict=False,
        predictor="mhcflurry",
        iedb_path=None,
        cedar_path=None,
        skip_ms_evidence=False,
        output=str(output_csv),
    )
    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_observations", return_value=hits_frame),
    ):
        cli_hits.handle(args)

    out = pd.read_csv(output_csv)
    # Output has the aggregator columns but no gene_ident_cols since
    # the input didn't provide them.
    assert "ms_pmhc_hit_count" in out.columns
    assert "gene_names" not in out.columns
    assert list(out["peptide"]) == ["QYIAQFTSQF"]
