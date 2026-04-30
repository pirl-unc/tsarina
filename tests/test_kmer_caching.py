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

"""Tests for tsarina's proteome k-mer filtering glue.

Tsarina owns:

1. The tsarina-specific CTA / non-CTA gene-ID partition frozensets
   (``_non_cta_gene_ids`` / ``_cta_gene_ids``) — memoized wrappers
   around :func:`tsarina.partition.CTA_partition_gene_ids`.
2. Candidate-driven CTA exclusivity filtering that streams non-CTA proteins
   against the much smaller CTA peptide set.
3. Thin viral call-site glue in ``human_exclusive_viral_peptides`` /
   ``cancer_specific_viral_peptides`` that invokes hitlist's primitive with
   the right gene filter.

These tests exercise partition memoization and pin down those dispatch
contracts.  The full hitlist proteome index is tested in hitlist's own suite.
"""

from __future__ import annotations

import pandas as pd
import pytest

from tsarina import peptides, viral


class _FakePartition:
    cta: set[str] = {"ENSG_A"}
    non_cta: set[str] = {"ENSG_B", "ENSG_C"}


@pytest.fixture(autouse=True)
def _reset_partition_caches(monkeypatch):
    peptides._non_cta_gene_ids.cache_clear()
    peptides._cta_gene_ids.cache_clear()
    monkeypatch.setattr("tsarina.partition.CTA_partition_gene_ids", lambda r: _FakePartition())
    yield


# ── Gene-ID partition frozenset wrappers ───────────────────────────────


def test_non_cta_gene_ids_returns_frozenset_of_non_cta_partition():
    result = peptides._non_cta_gene_ids(112)
    assert isinstance(result, frozenset)
    assert result == frozenset({"ENSG_B", "ENSG_C"})


def test_cta_gene_ids_returns_frozenset_of_cta_partition():
    result = peptides._cta_gene_ids(112)
    assert isinstance(result, frozenset)
    assert result == frozenset({"ENSG_A"})


def test_gene_id_wrappers_memoize_per_release():
    """Same release returns the same frozenset instance — avoids
    rebuilding the partition frozenset on every call."""
    assert peptides._non_cta_gene_ids(112) is peptides._non_cta_gene_ids(112)
    assert peptides._cta_gene_ids(112) is peptides._cta_gene_ids(112)


def test_gene_id_wrappers_separate_entries_per_release():
    """Different release → recompute (separate cache key)."""
    r1 = peptides._non_cta_gene_ids(112)
    r2 = peptides._non_cta_gene_ids(113)
    # Both frozensets should exist and match the stub (same content, but
    # distinct cache entries keyed on release).
    info = peptides._non_cta_gene_ids.cache_info()
    assert info.misses == 2
    assert r1 == r2  # stub ignores release; values are equal


def test_cta_gene_ids_for_names_limits_to_requested_symbols(monkeypatch):
    monkeypatch.setattr(
        "tsarina.gene_sets.CTA_gene_ids",
        lambda: {"ENSG_A", "ENSG_B", "ENSG_C"},
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.loader.cta_dataframe",
        lambda: pd.DataFrame(
            {
                "Symbol": ["MAGEA4", "PRAME;PRAME_ALT", "NONCTA"],
                "Ensembl_Gene_ID": ["ENSG_A", "ENSG_B;ENSG_C", "ENSG_D"],
            }
        ),
        raising=True,
    )

    assert peptides._cta_gene_ids_for_names(["PRAME"]) == ["ENSG_B", "ENSG_C"]


# ── Dispatch to hitlist.proteome.proteome_kmer_set ────────────────────


def test_cta_exclusive_peptides_streams_candidate_overlap(monkeypatch):
    """cta_exclusive_peptides should filter by streaming non-CTA overlaps,
    avoiding hitlist's full proteome_kmer_set posting-index build."""
    calls: list[dict] = []
    cta_calls: list[dict] = []

    def _hitlist_should_not_be_called(*args, **kwargs):
        raise AssertionError("cta_exclusive_peptides should not build full hitlist k-mer set")

    def _overlap_spy(candidate_peptides, **kw):
        calls.append(
            {
                "candidate_peptides": set(candidate_peptides),
                "ensembl_release": kw["ensembl_release"],
                "lengths": kw["lengths"],
            }
        )
        return {"UNWANTEDK"}

    monkeypatch.setattr("hitlist.proteome.proteome_kmer_set", _hitlist_should_not_be_called)
    monkeypatch.setattr(
        "tsarina.peptides._non_cta_overlapping_peptides",
        _overlap_spy,
        raising=True,
    )

    def _cta_peptides_spy(**kw):
        cta_calls.append(kw)
        return pd.DataFrame(
            {
                "peptide": ["GOODPEPTI", "UNWANTEDK"],
                "length": [9, 9],
                "gene_name": ["MAGEA4", "MAGEA4"],
                "gene_id": ["E1", "E1"],
            }
        )

    monkeypatch.setattr(
        "tsarina.peptides.cta_peptides",
        _cta_peptides_spy,
        raising=True,
    )

    out = peptides.cta_exclusive_peptides(
        ensembl_release=112,
        lengths=(9,),
        gene_names=["MAGEA4"],
    )

    assert cta_calls[0]["gene_names"] == ["MAGEA4"]
    assert len(calls) == 1
    assert calls[0]["ensembl_release"] == 112
    assert calls[0]["lengths"] == (9,)
    assert calls[0]["candidate_peptides"] == {"GOODPEPTI", "UNWANTEDK"}
    # The peptide that appeared in the streamed non-CTA overlap set must be filtered out.
    assert list(out["peptide"]) == ["GOODPEPTI"]


def test_non_cta_overlap_streaming_finds_candidate_in_transcript(monkeypatch):
    class _Transcript:
        biotype = "protein_coding"

        def __init__(self, protein_sequence: str):
            self.protein_sequence = protein_sequence

    class _Gene:
        def __init__(self, transcripts):
            self.transcripts = transcripts

    class _FakeEnsembl:
        def __init__(self, release):
            self.release = release

        def gene_by_id(self, gene_id):
            if gene_id == "ENSG_B":
                return _Gene([_Transcript("AAAAUNWANTEDKBBBB")])
            if gene_id == "ENSG_C":
                return _Gene([_Transcript("CCCCCCCCCCCC")])
            raise ValueError(gene_id)

    messages: list[str] = []
    monkeypatch.setattr("pyensembl.EnsemblRelease", _FakeEnsembl)

    seen = peptides._non_cta_overlapping_peptides(
        {"GOODPEPTI", "UNWANTEDK"},
        ensembl_release=112,
        lengths=(9,),
        on_progress=messages.append,
    )

    assert seen == {"UNWANTEDK"}
    assert any("Scanning 2 non-CTA genes" in msg for msg in messages)
    assert any("Found 1 CTA candidate peptides" in msg for msg in messages)


def test_human_exclusive_viral_peptides_delegates_with_none_gene_ids(monkeypatch):
    """The full-proteome call passes gene_ids=None so hitlist iterates
    every protein-coding gene rather than a subset."""
    calls: list[dict] = []

    def _spy(release, lengths, gene_ids=None, **kw):
        calls.append({"release": release, "lengths": lengths, "gene_ids": gene_ids})
        return frozenset({"HUMANSELF"})

    monkeypatch.setattr("hitlist.proteome.proteome_kmer_set", _spy, raising=True)
    monkeypatch.setattr(
        "tsarina.viral.viral_peptides",
        lambda **kw: pd.DataFrame(
            {
                "peptide": ["VIRALONLY", "HUMANSELF"],
                "length": [9, 9],
                "virus": ["test", "test"],
                "protein_id": ["P1", "P1"],
            }
        ),
        raising=True,
    )

    out = viral.human_exclusive_viral_peptides(virus="test", lengths=(9,), ensembl_release=112)

    assert len(calls) == 1
    assert calls[0]["gene_ids"] is None
    assert calls[0]["release"] == 112
    assert calls[0]["lengths"] == (9,)
    assert list(out["peptide"]) == ["VIRALONLY"]


def test_cancer_specific_viral_peptides_delegates_twice_non_cta_and_cta(monkeypatch):
    """cancer_specific_viral_peptides subtracts non-CTA k-mers AND
    annotates in-CTA-protein membership.  It must call the hitlist
    primitive twice — once per partition — with the right gene_ids
    each time."""
    calls: list[dict] = []

    def _spy(release, lengths, gene_ids=None, **kw):
        calls.append({"release": release, "lengths": lengths, "gene_ids": gene_ids})
        if gene_ids == frozenset({"ENSG_B", "ENSG_C"}):  # non-CTA
            return frozenset({"HUMANNCTA"})
        if gene_ids == frozenset({"ENSG_A"}):  # CTA
            return frozenset({"SHAREDCTA"})
        return frozenset()

    monkeypatch.setattr("hitlist.proteome.proteome_kmer_set", _spy, raising=True)
    monkeypatch.setattr(
        "tsarina.viral.viral_peptides",
        lambda **kw: pd.DataFrame(
            {
                "peptide": ["UNIQUEVIR", "HUMANNCTA", "SHAREDCTA"],
                "length": [9, 9, 9],
                "virus": ["test", "test", "test"],
                "protein_id": ["P1", "P1", "P1"],
            }
        ),
        raising=True,
    )

    out = viral.cancer_specific_viral_peptides(virus="test", lengths=(9,), ensembl_release=112)

    # Two calls — one per partition
    assert len(calls) == 2
    gene_ids_passed = {c["gene_ids"] for c in calls}
    assert gene_ids_passed == {
        frozenset({"ENSG_B", "ENSG_C"}),
        frozenset({"ENSG_A"}),
    }

    # HUMANNCTA dropped (in non-CTA); UNIQUEVIR and SHAREDCTA kept
    kept = set(out["peptide"])
    assert "HUMANNCTA" not in kept
    assert "UNIQUEVIR" in kept
    assert "SHAREDCTA" in kept
    # in_cta_protein flag correctly set
    shared_row = out[out["peptide"] == "SHAREDCTA"].iloc[0]
    unique_row = out[out["peptide"] == "UNIQUEVIR"].iloc[0]
    assert bool(shared_row["in_cta_protein"]) is True
    assert bool(unique_row["in_cta_protein"]) is False


def test_cancer_specific_viral_peptides_reuses_non_cta_frozenset(monkeypatch):
    """The viral path still passes the memoized non-CTA frozenset to hitlist."""
    received_gene_ids: list[frozenset] = []

    def _spy(release, lengths, gene_ids=None, **kw):
        if gene_ids is not None and gene_ids == frozenset({"ENSG_B", "ENSG_C"}):
            received_gene_ids.append(gene_ids)
        return frozenset()

    monkeypatch.setattr("hitlist.proteome.proteome_kmer_set", _spy, raising=True)
    monkeypatch.setattr(
        "tsarina.viral.viral_peptides",
        lambda **kw: pd.DataFrame(
            {
                "peptide": ["Y" * 9],
                "length": [9],
                "virus": ["test"],
                "protein_id": ["P1"],
            }
        ),
        raising=True,
    )

    viral.cancer_specific_viral_peptides(virus="test", lengths=(9,), ensembl_release=112)

    assert len(received_gene_ids) == 1
    assert received_gene_ids[0] is peptides._non_cta_gene_ids(112)
