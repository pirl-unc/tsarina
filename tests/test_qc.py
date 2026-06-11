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

"""Fragment-model QC (tsarina#108).

Flagging logic is exercised with an injected ``length_fn`` so it runs in CI,
which does not download pyensembl reference data. A single integration test
hits real pyensembl and skips when the data is unavailable.
"""

import pandas as pd
import pytest

from tsarina.qc import (
    MIN_CTA_PROTEIN_AA,
    find_fragment_gene_models,
    gene_max_protein_length,
)


def _df(rows):
    return pd.DataFrame(rows, columns=["Symbol", "Ensembl_Gene_ID"])


def _lengths(mapping):
    """Build a length_fn from a {gene_id: length|None} dict."""
    return lambda gene_id: mapping.get(gene_id)


def test_floor_separates_real_cta_from_fragment():
    # The floor must sit above the canonical fragment (GAGE12B, 7 aa) and below
    # the smallest real CTA in the bundled table (PRM1, 51 aa).
    assert 7 < MIN_CTA_PROTEIN_AA < 51


def test_flags_fragment_not_real_genes():
    df = _df(
        [
            ("GAGE12B", "ENSG00000236737"),  # 7-aa stub
            ("GAGE10", "ENSG00000215274"),  # 116-aa real GAGE
            ("MAGEA1", "ENSG00000198681"),  # full-length
        ]
    )
    length_fn = _lengths(
        {
            "ENSG00000236737": 7,
            "ENSG00000215274": 116,
            "ENSG00000198681": 309,
        }
    )
    flagged = find_fragment_gene_models(df, length_fn=length_fn)
    assert len(flagged) == 1
    assert flagged[0]["symbol"] == "GAGE12B"
    assert flagged[0]["gene_id"] == "ENSG00000236737"
    assert flagged[0]["protein_length"] == 7
    assert flagged[0]["reason"] == "fragment"


def test_small_but_real_cta_not_flagged():
    # PRM1 (protamine 1) is a legitimate 51-aa CTA; it must not be flagged.
    df = _df([("PRM1", "ENSG00000175646")])
    flagged = find_fragment_gene_models(df, length_fn=_lengths({"ENSG00000175646": 51}))
    assert flagged == []


def test_flags_gene_with_no_protein():
    df = _df([("FOO", "ENSG09999999999")])
    flagged = find_fragment_gene_models(df, length_fn=_lengths({"ENSG09999999999": None}))
    assert len(flagged) == 1
    assert flagged[0]["reason"] == "no_protein"
    assert flagged[0]["protein_length"] is None


def test_handles_versioned_and_multivalue_ids():
    # Versioned IDs are stripped; ;-separated cells are split and each checked.
    df = _df([("X", "ENSG00000236737.3"), ("Y", "ENSG00000215274;ENSG00000236737")])
    length_fn = _lengths({"ENSG00000236737": 7, "ENSG00000215274": 116})
    flagged = find_fragment_gene_models(df, length_fn=length_fn)
    # Both X and Y's second ID resolve to the 7-aa fragment.
    assert {f["symbol"] for f in flagged} == {"X", "Y"}
    assert all(f["gene_id"] == "ENSG00000236737" for f in flagged)


def test_clean_table_returns_empty():
    df = _df([("A", "ENSG1"), ("B", "ENSG2")])
    flagged = find_fragment_gene_models(df, length_fn=_lengths({"ENSG1": 120, "ENSG2": 300}))
    assert flagged == []


def test_custom_floor_is_respected():
    df = _df([("A", "ENSG1")])
    # 80 aa passes the default floor but trips a stricter 100-aa floor.
    assert find_fragment_gene_models(df, length_fn=_lengths({"ENSG1": 80})) == []
    flagged = find_fragment_gene_models(df, length_fn=_lengths({"ENSG1": 80}), min_protein_aa=100)
    assert len(flagged) == 1


# ── Integration: real pyensembl (skips when reference data is absent) ──────────


def _ensembl_or_skip():
    try:
        from pyensembl import EnsemblRelease

        ensembl = EnsemblRelease(112)
        ensembl.gene_by_id("ENSG00000198681")  # probe: triggers data load
        return ensembl
    except Exception as exc:  # data not downloaded in this environment (e.g. CI)
        pytest.skip(f"pyensembl release 112 data unavailable: {exc}")


def test_real_pyensembl_lengths_for_known_genes():
    ensembl = _ensembl_or_skip()
    # GAGE12B is the 7-aa fragment; GAGE10 is the 116-aa real GAGE.
    assert gene_max_protein_length("ENSG00000236737", ensembl=ensembl) == 7
    assert gene_max_protein_length("ENSG00000215274", ensembl=ensembl) == 116


def test_shipped_table_has_no_fragment_models():
    _ensembl_or_skip()
    # The bundled CTA universe must contain no fragment / no-protein gene models.
    flagged = find_fragment_gene_models()
    assert flagged == [], f"fragment gene models in shipped table: {flagged}"
