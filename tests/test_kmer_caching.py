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

"""Tests for the unified k-mer enumeration primitive (tsarina#12 + #19).

``peptides._proteome_kmers(release, lengths, gene_ids)`` walks protein-
coding transcripts in an Ensembl release and produces a ``frozenset`` of
every k-mer at the requested lengths within the gene subset.  The walk
is ~10-60s, so the helper is ``@lru_cache``'d on
``(release, lengths, gene_ids)``.

These tests exercise:

1. The unified ``_proteome_kmers`` primitive directly — cache behavior,
   gene-subset respect, full-proteome path (``gene_ids=None``).
2. The thin frozenset wrappers (``_non_cta_gene_ids`` /
   ``_cta_gene_ids``) that memoize the partition-derived sets.
3. End-to-end caching through the public callers
   (``cta_exclusive_peptides``, ``human_exclusive_viral_peptides``,
   ``cancer_specific_viral_peptides``).

A stubbed ``EnsemblRelease`` tracks ``.genes()`` / ``.gene_by_id()``
call counts so we can assert cache hits vs traversal without pulling
real Ensembl data.
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
    peptides._proteome_kmers.cache_clear()
    peptides._non_cta_gene_ids.cache_clear()
    peptides._cta_gene_ids.cache_clear()
    _FakeEnsembl.call_counts = {"genes": 0, "gene_by_id": 0}
    monkeypatch.setattr("pyensembl.EnsemblRelease", _FakeEnsembl)
    monkeypatch.setattr("tsarina.partition.CTA_partition_gene_ids", lambda r: _FakePartition())
    yield


# ── _proteome_kmers primitive ──────────────────────────────────────────


def test_proteome_kmers_is_lru_cached():
    assert hasattr(peptides._proteome_kmers, "cache_info")
    assert hasattr(peptides._proteome_kmers, "cache_clear")


def test_proteome_kmers_full_proteome_uses_genes_iteration():
    """gene_ids=None walks every gene via ensembl.genes() rather than
    looking up IDs one-by-one."""
    result = peptides._proteome_kmers(112, (8,), None)
    assert isinstance(result, frozenset)
    assert _FakeEnsembl.call_counts["genes"] == 1
    assert _FakeEnsembl.call_counts["gene_by_id"] == 0
    # All three fake genes are protein_coding → all appear.
    assert "MAAAAAAA" in result
    assert "MCCCCCCC" in result
    assert "MDDDDDDD" in result


def test_proteome_kmers_gene_subset_uses_gene_by_id():
    """A supplied gene_ids frozenset goes through gene_by_id()."""
    result = peptides._proteome_kmers(112, (8,), frozenset({"ENSG_B", "ENSG_C"}))
    assert "MCCCCCCC" in result
    assert "MDDDDDDD" in result
    assert "MAAAAAAA" not in result  # ENSG_A excluded by filter
    assert _FakeEnsembl.call_counts["genes"] == 0
    assert _FakeEnsembl.call_counts["gene_by_id"] == 2


def test_proteome_kmers_caches_by_gene_ids():
    """Same (release, lengths) with different gene_ids keys must produce
    separate cache entries; repeat lookups are cache hits."""
    peptides._proteome_kmers(112, (8,), frozenset({"ENSG_B"}))
    peptides._proteome_kmers(112, (8,), frozenset({"ENSG_C"}))
    info_after_two_distinct = peptides._proteome_kmers.cache_info()
    assert info_after_two_distinct.misses == 2
    assert info_after_two_distinct.hits == 0

    # Repeat — both should hit the cache
    peptides._proteome_kmers(112, (8,), frozenset({"ENSG_B"}))
    peptides._proteome_kmers(112, (8,), frozenset({"ENSG_C"}))
    info_after_repeat = peptides._proteome_kmers.cache_info()
    assert info_after_repeat.misses == 2
    assert info_after_repeat.hits == 2


def test_proteome_kmers_caches_by_lengths():
    peptides._proteome_kmers(112, (8,), None)
    peptides._proteome_kmers(112, (9,), None)
    info = peptides._proteome_kmers.cache_info()
    assert info.misses == 2
    assert info.hits == 0


# ── Gene-ID frozenset wrappers ─────────────────────────────────────────


def test_non_cta_gene_ids_returns_frozenset_of_non_cta_partition():
    result = peptides._non_cta_gene_ids(112)
    assert isinstance(result, frozenset)
    assert result == frozenset({"ENSG_B", "ENSG_C"})


def test_cta_gene_ids_returns_frozenset_of_cta_partition():
    result = peptides._cta_gene_ids(112)
    assert isinstance(result, frozenset)
    assert result == frozenset({"ENSG_A"})


def test_gene_id_wrappers_cache_per_release():
    """Calling the wrapper twice with the same release returns the same
    frozenset instance — so downstream _proteome_kmers lookups reuse
    the cached hash rather than rebuilding."""
    r1 = peptides._non_cta_gene_ids(112)
    r2 = peptides._non_cta_gene_ids(112)
    assert r1 is r2


# ── End-to-end caching through public callers ──────────────────────────


def test_cta_exclusive_peptides_hits_cache_on_repeat(monkeypatch):
    """The public wrapper must dispatch through _proteome_kmers — second
    call with the same args does NOT re-walk the proteome."""
    import pandas as pd

    monkeypatch.setattr(
        "tsarina.peptides.cta_peptides",
        lambda **kw: pd.DataFrame(
            {
                "peptide": ["MAAAAAAA", "NEVERSEEN"],
                "length": [8, 8],
                "gene_name": ["G1", "G1"],
                "gene_id": ["X", "X"],
            }
        ),
        raising=True,
    )
    peptides.cta_exclusive_peptides(ensembl_release=112, lengths=(8,))
    gene_by_id_after_first = _FakeEnsembl.call_counts["gene_by_id"]
    peptides.cta_exclusive_peptides(ensembl_release=112, lengths=(8,))
    assert _FakeEnsembl.call_counts["gene_by_id"] == gene_by_id_after_first


def test_human_exclusive_viral_peptides_hits_cache_on_repeat(monkeypatch):
    import pandas as pd

    monkeypatch.setattr(
        "tsarina.viral.viral_peptides",
        lambda **kw: pd.DataFrame(
            {
                "peptide": ["VIRALMER", "MAAAAAAA"],
                "length": [8, 8],
                "virus": ["test", "test"],
                "protein_id": ["P1", "P1"],
            }
        ),
        raising=True,
    )
    viral.human_exclusive_viral_peptides(virus="test", lengths=(8,), ensembl_release=112)
    genes_after_first = _FakeEnsembl.call_counts["genes"]
    viral.human_exclusive_viral_peptides(virus="test", lengths=(8,), ensembl_release=112)
    assert _FakeEnsembl.call_counts["genes"] == genes_after_first


def test_cancer_specific_viral_peptides_hits_cache_on_repeat(monkeypatch):
    """PR D headline: cancer_specific_viral_peptides now shares the
    cached _proteome_kmers primitive instead of walking the proteome
    on every call."""
    import pandas as pd

    monkeypatch.setattr(
        "tsarina.viral.viral_peptides",
        lambda **kw: pd.DataFrame(
            {
                "peptide": ["VIRALMER", "MAAAAAAA", "MCCCCCCC"],
                "length": [8, 8, 8],
                "virus": ["test", "test", "test"],
                "protein_id": ["P1", "P1", "P1"],
            }
        ),
        raising=True,
    )
    result_1 = viral.cancer_specific_viral_peptides(virus="test", lengths=(8,), ensembl_release=112)
    gene_by_id_after_first = _FakeEnsembl.call_counts["gene_by_id"]
    assert gene_by_id_after_first > 0  # first call did real work

    result_2 = viral.cancer_specific_viral_peptides(virus="test", lengths=(8,), ensembl_release=112)
    # Second call hits the cache — no additional gene_by_id lookups
    assert _FakeEnsembl.call_counts["gene_by_id"] == gene_by_id_after_first

    # Functional equivalence between the two calls
    assert list(result_1["peptide"]) == list(result_2["peptide"])
    # in_cta_protein flag set correctly: MAAAAAAA comes from ENSG_A (CTA),
    # VIRALMER is viral-only, MCCCCCCC comes from ENSG_B (non-CTA → dropped
    # by the noncta subtraction).
    kept = set(result_1["peptide"])
    assert "VIRALMER" in kept
    assert "MAAAAAAA" in kept  # in CTA proteins, not in non-CTA — kept
    assert "MCCCCCCC" not in kept  # in non-CTA proteins — filtered
    maaaaaaa_row = result_1[result_1["peptide"] == "MAAAAAAA"].iloc[0]
    assert bool(maaaaaaa_row["in_cta_protein"]) is True
    viralmer_row = result_1[result_1["peptide"] == "VIRALMER"].iloc[0]
    assert bool(viralmer_row["in_cta_protein"]) is False


def test_cta_and_non_cta_share_proteome_kmers_cache(monkeypatch):
    """cta_exclusive_peptides and cancer_specific_viral_peptides both
    query _proteome_kmers(non_cta). The second function should hit the
    cache populated by the first — no re-walk across functions."""
    import pandas as pd

    monkeypatch.setattr(
        "tsarina.peptides.cta_peptides",
        lambda **kw: pd.DataFrame(
            {
                "peptide": ["MAAAAAAA"],
                "length": [8],
                "gene_name": ["G1"],
                "gene_id": ["X"],
            }
        ),
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.viral.viral_peptides",
        lambda **kw: pd.DataFrame(
            {
                "peptide": ["VIRALMER"],
                "length": [8],
                "virus": ["test"],
                "protein_id": ["P1"],
            }
        ),
        raising=True,
    )
    peptides.cta_exclusive_peptides(ensembl_release=112, lengths=(8,))
    calls_after_cta = _FakeEnsembl.call_counts["gene_by_id"]
    viral.cancer_specific_viral_peptides(virus="test", lengths=(8,), ensembl_release=112)
    # cancer_specific_viral_peptides additionally needs the CTA k-mer set,
    # so it will do ONE new walk (for the CTA partition) — but NOT re-walk
    # the non-CTA partition that cta_exclusive_peptides already cached.
    # CTA partition has 1 gene (ENSG_A) → 1 additional gene_by_id call.
    assert _FakeEnsembl.call_counts["gene_by_id"] == calls_after_cta + 1
