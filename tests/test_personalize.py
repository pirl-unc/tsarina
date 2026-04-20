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

"""Unit tests for :mod:`tsarina.personalize`.

Exercises the modernized ``personalize()`` path: CTA/viral exclusivity gates,
explicit tier assignment, deterministic sort, tumor-specificity filter,
predictor warning, and mandatory-prediction failure propagation.
"""

from __future__ import annotations

import warnings

import pandas as pd
import pytest

from tsarina import personalize as personalize_module
from tsarina.personalize import _assign_tiers, personalize

# ── Tier assignment unit tests ──────────────────────────────────────────


def test_tier1_requires_cancer_ms_for_cta():
    df = pd.DataFrame(
        {
            "peptide": ["P1", "P2"],
            "category": ["cta", "cta"],
            "presentation_percentile": [0.1, 0.1],
            "ms_in_cancer": [True, False],
            "ms_hit_count": [5, 5],
        }
    )
    out = _assign_tiers(df)
    assert list(out["tier"]) == [1, 2]
    assert out.loc[0, "tier_label"] == "STRONG"


def test_tier1_automatic_for_mutant_and_viral_at_strong_percentile():
    df = pd.DataFrame(
        {
            "peptide": ["M1", "V1"],
            "category": ["mutant", "viral"],
            "presentation_percentile": [0.2, 0.2],
            "ms_in_cancer": [False, False],
            "ms_hit_count": [0, 0],
        }
    )
    out = _assign_tiers(df)
    assert list(out["tier"]) == [1, 1]


def test_tier3_when_presented_but_no_evidence_and_self():
    df = pd.DataFrame(
        {
            "peptide": ["C1"],
            "category": ["cta"],
            "presentation_percentile": [0.9],
            "ms_in_cancer": [False],
            "ms_hit_count": [0],
        }
    )
    out = _assign_tiers(df)
    assert out.loc[0, "tier"] == 3
    assert out.loc[0, "tier_reason"] == "presentation_only"


def test_tier4_for_nan_percentile():
    df = pd.DataFrame(
        {
            "peptide": ["X1"],
            "category": ["cta"],
            "presentation_percentile": [float("nan")],
            "ms_in_cancer": [True],
            "ms_hit_count": [3],
        }
    )
    out = _assign_tiers(df)
    assert out.loc[0, "tier"] == 4
    assert out.loc[0, "tier_reason"] == "unscored"


def test_tier4_for_high_percentile():
    df = pd.DataFrame(
        {
            "peptide": ["X2"],
            "category": ["mutant"],
            "presentation_percentile": [5.0],
            "ms_in_cancer": [False],
            "ms_hit_count": [0],
        }
    )
    out = _assign_tiers(df)
    assert out.loc[0, "tier"] == 4
    assert out.loc[0, "tier_reason"] == "below_threshold"


# ── Integration: CTA exclusivity gate ──────────────────────────────────


def test_personalize_uses_cta_exclusive_peptides(monkeypatch):
    """A peptide present in `cta_peptides` but NOT in `cta_exclusive_peptides`
    must be dropped before scoring — confirms the PR-A fix."""
    exclusive = pd.DataFrame(
        {
            "peptide": ["EXCLUSIVE9"],
            "length": [9],
            "gene_name": ["MAGEA4"],
            "gene_id": ["ENSG_MAGEA4"],
        }
    )
    called = {}

    def _fake_cta_exclusive(**kwargs):
        called["ok"] = True
        return exclusive

    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides", _fake_cta_exclusive, raising=True
    )
    monkeypatch.setattr("tsarina.gene_sets.CTA_gene_names", lambda: {"MAGEA4"}, raising=True)
    monkeypatch.setattr("tsarina.gene_sets.CTA_by_axes", lambda **kw: {"MAGEA4"}, raising=True)

    out = personalize(
        hla_alleles=["HLA-A*02:01"],
        cta_expression={"MAGEA4": 10.0},
        score_presentation=False,
        skip_ms_evidence=True,
        drop_weak_tier=False,
    )
    assert called.get("ok") is True
    assert list(out["peptide"]) == ["EXCLUSIVE9"]


# ── Integration: viral exclusivity + enforce_tumor_specificity ─────────


def test_personalize_uses_human_exclusive_viral_by_default(monkeypatch):
    """Confirms PR-B: default viral path routes through the exclusive fn."""
    exclusive = pd.DataFrame(
        {
            "peptide": ["VIRAL9ONLY"],
            "length": [9],
            "virus": ["hpv16"],
            "protein_id": ["E6"],
        }
    )
    touched = {"exclusive": False, "plain": False}

    def _fake_exclusive(**kwargs):
        touched["exclusive"] = True
        return exclusive

    def _fake_plain(**kwargs):
        touched["plain"] = True
        return exclusive

    monkeypatch.setattr(
        "tsarina.viral.human_exclusive_viral_peptides", _fake_exclusive, raising=True
    )
    monkeypatch.setattr("tsarina.viral.viral_peptides", _fake_plain, raising=True)

    out = personalize(
        hla_alleles=["HLA-A*02:01"],
        viruses=["hpv16"],
        score_presentation=False,
        skip_ms_evidence=True,
        drop_weak_tier=False,
    )
    assert touched["exclusive"] is True
    assert touched["plain"] is False
    assert "VIRAL9ONLY" in out["peptide"].tolist()


def test_personalize_no_require_human_exclusive_viral_uses_plain(monkeypatch):
    plain = pd.DataFrame(
        {
            "peptide": ["PEPV"],
            "length": [4],
            "virus": ["hpv16"],
            "protein_id": ["E6"],
        }
    )
    touched = {"exclusive": False, "plain": False}

    def _fake_exclusive(**kwargs):
        touched["exclusive"] = True
        return plain

    def _fake_plain(**kwargs):
        touched["plain"] = True
        return plain

    monkeypatch.setattr(
        "tsarina.viral.human_exclusive_viral_peptides", _fake_exclusive, raising=True
    )
    monkeypatch.setattr("tsarina.viral.viral_peptides", _fake_plain, raising=True)

    personalize(
        hla_alleles=["HLA-A*02:01"],
        viruses=["hpv16"],
        require_human_exclusive_viral=False,
        score_presentation=False,
        skip_ms_evidence=True,
        drop_weak_tier=False,
    )
    assert touched["plain"] is True
    assert touched["exclusive"] is False


# ── enforce_tumor_specificity drops healthy-tissue MS hits ─────────────


def test_enforce_tumor_specificity_drops_healthy_tissue_hit(monkeypatch):
    exclusive = pd.DataFrame(
        {
            "peptide": ["SAFE9PEP", "RISKYPEP9"],
            "length": [9, 9],
            "virus": ["hpv16", "hpv16"],
            "protein_id": ["E6", "E6"],
        }
    )
    monkeypatch.setattr(
        "tsarina.viral.human_exclusive_viral_peptides",
        lambda **kw: exclusive,
        raising=True,
    )

    def _fake_load_ms_evidence(peptides, **kwargs):
        return pd.DataFrame(
            {
                "peptide": ["RISKYPEP9"],
                "mhc_restriction": ["HLA-A*02:01"],
                "src_cancer": [False],
                "src_healthy_tissue": [True],
            }
        )

    monkeypatch.setattr("tsarina.indexing.load_ms_evidence", _fake_load_ms_evidence, raising=True)

    out = personalize(
        hla_alleles=["HLA-A*02:01"],
        viruses=["hpv16"],
        score_presentation=False,
        enforce_tumor_specificity=True,
        drop_weak_tier=False,
    )
    assert "RISKYPEP9" not in out["peptide"].tolist()
    assert "SAFE9PEP" in out["peptide"].tolist()


# ── Predictor warning ──────────────────────────────────────────────────


def test_non_mhcflurry_predictor_emits_warning(monkeypatch):
    """PR-H: warn when --predictor ≠ mhcflurry since tier cutoffs are
    mhcflurry-calibrated."""

    def _fake_score(**kwargs):
        return pd.DataFrame(
            columns=[
                "peptide",
                "allele",
                "presentation_score",
                "presentation_percentile",
                "affinity_nm",
            ]
        )

    monkeypatch.setattr(personalize_module, "_OUTPUT_COLUMNS", personalize_module._OUTPUT_COLUMNS)
    monkeypatch.setattr("tsarina.scoring.score_presentation", _fake_score, raising=True)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        personalize(
            hla_alleles=["HLA-A*02:01"],
            score_presentation=True,
            predictor="netmhcpan",
            skip_ms_evidence=True,
        )
    messages = [str(m.message) for m in w if issubclass(m.category, UserWarning)]
    assert any("netmhcpan" in msg for msg in messages)


# ── Mandatory MHC gate: topiary ImportError propagates ─────────────────


def test_score_presentation_import_error_propagates(monkeypatch):
    """PR-C: when scoring is requested and the backend is missing, raise
    (do not silently continue)."""
    exclusive = pd.DataFrame(
        {
            "peptide": ["P1"],
            "length": [9],
            "virus": ["hpv16"],
            "protein_id": ["E6"],
        }
    )
    monkeypatch.setattr(
        "tsarina.viral.human_exclusive_viral_peptides",
        lambda **kw: exclusive,
        raising=True,
    )

    def _boom(**kwargs):
        raise ImportError("topiary is not installed in this test env")

    monkeypatch.setattr("tsarina.scoring.score_presentation", _boom, raising=True)

    with pytest.raises(ImportError, match="topiary"):
        personalize(
            hla_alleles=["HLA-A*02:01"],
            viruses=["hpv16"],
            score_presentation=True,
            skip_ms_evidence=True,
        )


# ── drop_weak_tier removes T4 rows ─────────────────────────────────────


def test_drop_weak_tier_removes_tier4(monkeypatch):
    exclusive = pd.DataFrame(
        {
            "peptide": ["GOODPEP9", "WEAKPEP99"],
            "length": [9, 9],
            "virus": ["hpv16", "hpv16"],
            "protein_id": ["E6", "E6"],
        }
    )
    monkeypatch.setattr(
        "tsarina.viral.human_exclusive_viral_peptides",
        lambda **kw: exclusive,
        raising=True,
    )

    def _fake_score(**kwargs):
        return pd.DataFrame(
            {
                "peptide": ["GOODPEP9", "WEAKPEP99"],
                "allele": ["HLA-A*02:01", "HLA-A*02:01"],
                "presentation_score": [0.95, 0.05],
                "presentation_percentile": [0.1, 20.0],
                "affinity_nm": [10.0, 30000.0],
            }
        )

    monkeypatch.setattr("tsarina.scoring.score_presentation", _fake_score, raising=True)

    out_kept = personalize(
        hla_alleles=["HLA-A*02:01"],
        viruses=["hpv16"],
        score_presentation=True,
        skip_ms_evidence=True,
        drop_weak_tier=True,
    )
    assert "WEAKPEP99" not in out_kept["peptide"].tolist()
    assert "GOODPEP9" in out_kept["peptide"].tolist()

    out_all = personalize(
        hla_alleles=["HLA-A*02:01"],
        viruses=["hpv16"],
        score_presentation=True,
        skip_ms_evidence=True,
        drop_weak_tier=False,
    )
    assert set(out_all["peptide"]) == {"GOODPEP9", "WEAKPEP99"}


# ── Deterministic sort ─────────────────────────────────────────────────


def test_sort_is_deterministic(monkeypatch):
    """Equal tier + hit count + percentile must still produce a stable
    order (peptide ASC, best_allele ASC)."""
    exclusive = pd.DataFrame(
        {
            "peptide": ["ZZZZZ9PEP", "AAAAA9PEP", "MMMMM9PEP"],
            "length": [9, 9, 9],
            "virus": ["hpv16", "hpv16", "hpv16"],
            "protein_id": ["E6", "E7", "E6"],
        }
    )
    monkeypatch.setattr(
        "tsarina.viral.human_exclusive_viral_peptides",
        lambda **kw: exclusive,
        raising=True,
    )

    def _fake_score(peptides, alleles, **kwargs):
        # All tied at 0.1 percentile → sort falls through to peptide/allele
        return pd.DataFrame(
            {
                "peptide": list(peptides),
                "allele": ["HLA-A*02:01"] * len(peptides),
                "presentation_score": [0.95] * len(peptides),
                "presentation_percentile": [0.1] * len(peptides),
                "affinity_nm": [10.0] * len(peptides),
            }
        )

    monkeypatch.setattr("tsarina.scoring.score_presentation", _fake_score, raising=True)

    a = personalize(
        hla_alleles=["HLA-A*02:01"],
        viruses=["hpv16"],
        score_presentation=True,
        skip_ms_evidence=True,
    )
    b = personalize(
        hla_alleles=["HLA-A*02:01"],
        viruses=["hpv16"],
        score_presentation=True,
        skip_ms_evidence=True,
    )
    assert list(a["peptide"]) == list(b["peptide"])
    # Expect peptide-ASC ordering when ties extend through percentile.
    assert list(a["peptide"]) == sorted(a["peptide"].tolist())


# ── Empty input paths return the canonical column set ─────────────────


def test_empty_inputs_return_canonical_columns():
    out = personalize(
        hla_alleles=["HLA-A*02:01"],
        score_presentation=False,
        skip_ms_evidence=True,
    )
    from tsarina.personalize import _OUTPUT_COLUMNS

    assert list(out.columns) == list(_OUTPUT_COLUMNS)
    assert out.empty


# ── enforce_tumor_specificity=False keeps healthy-tissue hits ──────────


def test_enforce_tumor_specificity_false_keeps_healthy_tissue_hit(monkeypatch):
    """Opt-out path: user explicitly wants healthy-tissue hits retained
    (e.g. for diagnostic inspection)."""
    exclusive = pd.DataFrame(
        {
            "peptide": ["SAFE9PEP", "RISKYPEP9"],
            "length": [9, 9],
            "virus": ["hpv16", "hpv16"],
            "protein_id": ["E6", "E6"],
        }
    )
    monkeypatch.setattr(
        "tsarina.viral.human_exclusive_viral_peptides",
        lambda **kw: exclusive,
        raising=True,
    )

    def _fake_load_ms_evidence(peptides, **kwargs):
        return pd.DataFrame(
            {
                "peptide": ["RISKYPEP9"],
                "mhc_restriction": ["HLA-A*02:01"],
                "src_cancer": [False],
                "src_healthy_tissue": [True],
            }
        )

    monkeypatch.setattr("tsarina.indexing.load_ms_evidence", _fake_load_ms_evidence, raising=True)

    out = personalize(
        hla_alleles=["HLA-A*02:01"],
        viruses=["hpv16"],
        score_presentation=False,
        enforce_tumor_specificity=False,
        drop_weak_tier=False,
    )
    peptides = set(out["peptide"])
    assert "RISKYPEP9" in peptides
    assert "SAFE9PEP" in peptides
    # Safety flag is still surfaced on the row even when not filtered.
    risky_row = out[out["peptide"] == "RISKYPEP9"].iloc[0]
    assert bool(risky_row["ms_in_healthy_tissue"]) is True


# ── mtec_matrix_path gates CTAs by thymic expression ──────────────────


def test_mtec_matrix_path_filters_ctas(tmp_path, monkeypatch):
    """PR-F: when an mTEC matrix is provided, CTAs with mean mTEC TPM above
    mtec_max_tpm must be dropped. Confirms filter_by_mtec is actually wired
    in (personalize.py:222-225 was uncovered)."""
    mtec_tsv = tmp_path / "mtec.tsv"
    mtec_tsv.write_text(
        "gene_symbol\tsample1\tsample2\nLOWMTEC_GENE\t0.1\t0.2\nHIGHMTEC_GENE\t50.0\t60.0\n"
    )

    exclusive = pd.DataFrame(
        {
            "peptide": ["LOWPEP9LEN", "HIGHPEP9LN"],
            "length": [9, 9],
            "gene_name": ["LOWMTEC_GENE", "HIGHMTEC_GENE"],
            "gene_id": ["ENSG_L", "ENSG_H"],
        }
    )
    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides", lambda **kw: exclusive, raising=True
    )
    monkeypatch.setattr(
        "tsarina.gene_sets.CTA_gene_names",
        lambda: {"LOWMTEC_GENE", "HIGHMTEC_GENE"},
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.gene_sets.CTA_by_axes",
        lambda **kw: {"LOWMTEC_GENE", "HIGHMTEC_GENE"},
        raising=True,
    )

    out = personalize(
        hla_alleles=["HLA-A*02:01"],
        cta_expression={"LOWMTEC_GENE": 10.0, "HIGHMTEC_GENE": 10.0},
        mtec_matrix_path=str(mtec_tsv),
        mtec_max_tpm=1.0,
        score_presentation=False,
        skip_ms_evidence=True,
        drop_weak_tier=False,
    )
    peptides = set(out["peptide"])
    assert "LOWPEP9LEN" in peptides
    assert "HIGHPEP9LN" not in peptides


# ── min_restriction_confidence=None bypasses the confidence gate ──────


def test_min_restriction_confidence_none_bypasses_gate(monkeypatch):
    """PR-G opt-out: passing None must skip CTA_by_axes entirely and admit
    genes regardless of restriction_confidence (LOW-confidence CTAs allowed)."""
    exclusive = pd.DataFrame(
        {
            "peptide": ["LOWCONF9PEP"],
            "length": [9],
            "gene_name": ["LOWCONF_CTA"],
            "gene_id": ["ENSG_LC"],
        }
    )

    axes_calls: list[dict] = []

    def _fake_by_axes(**kwargs):
        axes_calls.append(kwargs)
        # Pretend this gene would NOT pass a HIGH/MODERATE filter if asked.
        return set()

    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides", lambda **kw: exclusive, raising=True
    )
    monkeypatch.setattr("tsarina.gene_sets.CTA_gene_names", lambda: {"LOWCONF_CTA"}, raising=True)
    monkeypatch.setattr("tsarina.gene_sets.CTA_by_axes", _fake_by_axes, raising=True)

    out = personalize(
        hla_alleles=["HLA-A*02:01"],
        cta_expression={"LOWCONF_CTA": 10.0},
        min_restriction_confidence=None,
        score_presentation=False,
        skip_ms_evidence=True,
        drop_weak_tier=False,
    )
    assert axes_calls == [], "CTA_by_axes must not be called when confidence gate is disabled"
    assert "LOWCONF9PEP" in out["peptide"].tolist()


# ── Non-mhcflurry predictor still populates tiers correctly ───────────


def test_non_mhcflurry_predictor_populates_tiers(monkeypatch):
    """Combines PR-H warning with end-to-end tier correctness: the warning
    should not compromise tier assignment for an otherwise-strong candidate."""
    exclusive = pd.DataFrame(
        {
            "peptide": ["STRONGVIR9"],
            "length": [9],
            "virus": ["hpv16"],
            "protein_id": ["E6"],
        }
    )
    monkeypatch.setattr(
        "tsarina.viral.human_exclusive_viral_peptides",
        lambda **kw: exclusive,
        raising=True,
    )

    def _fake_score(**kwargs):
        return pd.DataFrame(
            {
                "peptide": ["STRONGVIR9"],
                "allele": ["HLA-A*02:01"],
                "presentation_score": [0.98],
                "presentation_percentile": [0.1],
                "affinity_nm": [8.0],
            }
        )

    monkeypatch.setattr("tsarina.scoring.score_presentation", _fake_score, raising=True)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out = personalize(
            hla_alleles=["HLA-A*02:01"],
            viruses=["hpv16"],
            predictor="netmhcpan_el",
            score_presentation=True,
            skip_ms_evidence=True,
        )

    assert any("netmhcpan_el" in str(m.message) for m in w if issubclass(m.category, UserWarning))
    assert len(out) == 1
    row = out.iloc[0]
    assert int(row["tier"]) == 1  # viral @ 0.1 percentile → STRONG
    assert row["tier_label"] == "STRONG"
    assert row["tier_reason"] == "strong_presentation+viral"
    assert row["best_allele"] == "HLA-A*02:01"
    assert float(row["presentation_percentile"]) == 0.1
