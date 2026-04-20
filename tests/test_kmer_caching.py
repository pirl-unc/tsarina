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

"""Tests for k-mer enumeration caching (tsarina#12).

``_non_cta_proteome_kmers`` and ``_human_proteome_kmers`` walk every
protein-coding transcript in an Ensembl release and produce a frozenset
of all k-mers at the requested lengths. The walk is ~10-60s depending on
release + length mix, so both helpers are ``@lru_cache``'d keyed on
``(ensembl_release, lengths)``.

These tests exercise the caching behavior without actually loading
Ensembl data: pyensembl is monkeypatched to a stub that tracks how many
times the proteome is traversed, so we can directly assert that two
calls with the same arguments produce a single traversal and that a
different argument triggers a second.
"""

from __future__ import annotations

import pytest

from tsarina import peptides, viral


class _FakeTranscript:
    def __init__(self, seq: str, biotype: str = "protein_coding"):
        self.protein_sequence = seq
        self.biotype = biotype


class _FakeGene:
    def __init__(self, gene_id: str, seq: str, biotype: str = "protein_coding"):
        self.gene_id = gene_id
        self.biotype = biotype
        self.transcripts = [_FakeTranscript(seq, biotype)]


class _FakeEnsembl:
    """Records how many times .genes() and .gene_by_id() are invoked."""

    call_counts: dict[str, int] = {"genes": 0, "gene_by_id": 0}

    def __init__(self, release: int):
        self.release = release
        # Fake universe: a handful of "protein" sequences long enough to
        # yield several k-mers at length 8.
        self._genes_by_id = {
            "ENSG_A": _FakeGene("ENSG_A", "MAAAAAAAAA"),
            "ENSG_B": _FakeGene("ENSG_B", "MCCCCCCCCC"),
            "ENSG_C": _FakeGene("ENSG_C", "MDDDDDDDDD"),
        }

    def genes(self):
        type(self).call_counts["genes"] += 1
        return list(self._genes_by_id.values())

    def gene_by_id(self, gene_id: str):
        type(self).call_counts["gene_by_id"] += 1
        if gene_id not in self._genes_by_id:
            raise ValueError(gene_id)
        return self._genes_by_id[gene_id]


class _FakePartition:
    cta: set[str] = {"ENSG_A"}
    non_cta: set[str] = {"ENSG_B", "ENSG_C"}


@pytest.fixture(autouse=True)
def _reset_caches_and_counts(monkeypatch):
    peptides._non_cta_proteome_kmers.cache_clear()
    viral._human_proteome_kmers.cache_clear()
    _FakeEnsembl.call_counts = {"genes": 0, "gene_by_id": 0}
    monkeypatch.setattr("pyensembl.EnsemblRelease", _FakeEnsembl)
    monkeypatch.setattr("tsarina.partition.CTA_partition_gene_ids", lambda r: _FakePartition())
    yield


# ── _non_cta_proteome_kmers (peptides.py) ──────────────────────────────


def test_non_cta_helper_is_lru_cached():
    assert hasattr(peptides._non_cta_proteome_kmers, "cache_info")
    assert hasattr(peptides._non_cta_proteome_kmers, "cache_clear")


def test_non_cta_returns_frozenset_of_kmers_from_non_cta_genes():
    result = peptides._non_cta_proteome_kmers(112, (8,))
    assert isinstance(result, frozenset)
    # ENSG_B = "MCCCCCCCCC" (10 chars) → three 8-mers: MCCCCCCC, CCCCCCCC, CCCCCCCC
    # But CCCCCCCC appears as both second and third (same k-mer, deduped).
    # ENSG_C = "MDDDDDDDDD" → analogous.
    # ENSG_A is CTA → excluded from this helper.
    assert "MCCCCCCC" in result
    assert "MDDDDDDD" in result
    assert "MAAAAAAA" not in result  # ENSG_A is CTA, excluded


def test_non_cta_caches_across_calls_with_same_args():
    r1 = peptides._non_cta_proteome_kmers(112, (8,))
    gene_by_id_calls_after_first = _FakeEnsembl.call_counts["gene_by_id"]
    assert gene_by_id_calls_after_first > 0  # first call did real work

    r2 = peptides._non_cta_proteome_kmers(112, (8,))
    # Second call hit cache — no additional traversal
    assert _FakeEnsembl.call_counts["gene_by_id"] == gene_by_id_calls_after_first
    # Identical frozenset object (lru_cache returns cached value)
    assert r1 is r2


def test_non_cta_different_args_compute_separately():
    peptides._non_cta_proteome_kmers(112, (8,))
    peptides._non_cta_proteome_kmers(112, (9,))  # different lengths → new cache key
    info = peptides._non_cta_proteome_kmers.cache_info()
    assert info.misses == 2
    assert info.hits == 0


# ── _human_proteome_kmers (viral.py) ───────────────────────────────────


def test_human_helper_is_lru_cached():
    assert hasattr(viral._human_proteome_kmers, "cache_info")
    assert hasattr(viral._human_proteome_kmers, "cache_clear")


def test_human_returns_frozenset_of_kmers_from_all_protein_coding_genes():
    result = viral._human_proteome_kmers(112, (8,))
    assert isinstance(result, frozenset)
    # ENSG_A (CTA) + ENSG_B + ENSG_C are all protein_coding, all included here
    # (the human helper does NOT exclude CTAs — that's the whole proteome).
    assert "MAAAAAAA" in result
    assert "MCCCCCCC" in result
    assert "MDDDDDDD" in result


def test_human_caches_across_calls_with_same_args():
    r1 = viral._human_proteome_kmers(112, (8,))
    genes_calls_after_first = _FakeEnsembl.call_counts["genes"]
    assert genes_calls_after_first == 1

    r2 = viral._human_proteome_kmers(112, (8,))
    # Second call hit cache — .genes() not called a second time
    assert _FakeEnsembl.call_counts["genes"] == 1
    assert r1 is r2


def test_human_different_release_triggers_recompute():
    viral._human_proteome_kmers(112, (8,))
    viral._human_proteome_kmers(113, (8,))  # different release → new cache key
    info = viral._human_proteome_kmers.cache_info()
    assert info.misses == 2
    assert info.hits == 0
