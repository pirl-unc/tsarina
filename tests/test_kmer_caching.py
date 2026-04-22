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

"""Tests for tsarina's delegation to ``hitlist.proteome.proteome_kmer_set``.

As of v1.1.0 the proteome-walk primitive lives in hitlist (shipped in
hitlist 1.14.2).  Tsarina owns only:

1. The tsarina-specific CTA / non-CTA gene-ID partition frozensets
   (``_non_cta_gene_ids`` / ``_cta_gene_ids``) — memoized wrappers
   around :func:`tsarina.partition.CTA_partition_gene_ids`.
2. Thin call-site glue in ``cta_exclusive_peptides`` /
   ``human_exclusive_viral_peptides`` / ``cancer_specific_viral_peptides``
   that invokes the hitlist primitive with the right gene filter.

These tests exercise the partition memoization and pin down the
dispatch contract (hitlist gets the right ``gene_ids`` frozenset per
caller).  The proteome walk itself is tested in hitlist's own test
suite — tsarina doesn't re-test the primitive.
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


# ── Dispatch to hitlist.proteome.proteome_kmer_set ────────────────────


def test_cta_exclusive_peptides_delegates_to_hitlist_with_non_cta_gene_ids(monkeypatch):
    """cta_exclusive_peptides must call hitlist.proteome.proteome_kmer_set
    with gene_ids=_non_cta_gene_ids(release), not rebuild the walk
    locally."""
    calls: list[dict] = []

    def _spy(release, lengths, gene_ids=None, **kw):
        calls.append({"release": release, "lengths": lengths, "gene_ids": gene_ids})
        return frozenset({"UNWANTEDK"})

    monkeypatch.setattr("hitlist.proteome.proteome_kmer_set", _spy, raising=True)
    monkeypatch.setattr(
        "tsarina.peptides.cta_peptides",
        lambda **kw: pd.DataFrame(
            {
                "peptide": ["GOODPEPTI", "UNWANTEDK"],
                "length": [9, 9],
                "gene_name": ["MAGEA4", "MAGEA4"],
                "gene_id": ["E1", "E1"],
            }
        ),
        raising=True,
    )

    out = peptides.cta_exclusive_peptides(ensembl_release=112, lengths=(9,))

    assert len(calls) == 1
    assert calls[0]["release"] == 112
    assert calls[0]["lengths"] == (9,)
    assert calls[0]["gene_ids"] == frozenset({"ENSG_B", "ENSG_C"})
    # The peptide that appeared in the hitlist-returned background set
    # must be filtered out.
    assert list(out["peptide"]) == ["GOODPEPTI"]


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


def test_non_cta_gene_ids_frozenset_is_reused_across_callers(monkeypatch):
    """Both cta_exclusive_peptides and cancer_specific_viral_peptides
    query the non-CTA partition.  They should pass the SAME frozenset
    instance so hitlist's cache keying is stable across callers."""
    received_gene_ids: list[frozenset] = []

    def _spy(release, lengths, gene_ids=None, **kw):
        if gene_ids is not None and gene_ids == frozenset({"ENSG_B", "ENSG_C"}):
            received_gene_ids.append(gene_ids)
        return frozenset()

    monkeypatch.setattr("hitlist.proteome.proteome_kmer_set", _spy, raising=True)
    monkeypatch.setattr(
        "tsarina.peptides.cta_peptides",
        lambda **kw: pd.DataFrame(
            {
                "peptide": ["X" * 9],
                "length": [9],
                "gene_name": ["MAGEA4"],
                "gene_id": ["E1"],
            }
        ),
        raising=True,
    )
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

    peptides.cta_exclusive_peptides(ensembl_release=112, lengths=(9,))
    viral.cancer_specific_viral_peptides(virus="test", lengths=(9,), ensembl_release=112)

    assert len(received_gene_ids) == 2
    # Both callers receive the SAME frozenset object (not just equal —
    # identical), so hitlist's cache treats them as one lookup.
    assert received_gene_ids[0] is received_gene_ids[1]
