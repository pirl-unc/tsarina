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

"""Unit tests for ``tsarina.spanning.spanning_pmhc_set``.

The library function pulls peptides through ``cta_exclusive_peptides``,
loads public MS evidence, scores them via ``score_presentation``, and
pivots the best result per (CTA, allele) cell. Tests stub each upstream
call so they run without Ensembl / hitlist data / mhcflurry installed.
"""

from __future__ import annotations

import io

import pandas as pd
import pytest

from tsarina.spanning import _DEFAULT_RANK_COLUMN, _resolve_ctas, panel_summary, spanning_pmhc_set


def _stub_peptides() -> pd.DataFrame:
    """Three CTAs, two peptides each at length 9."""
    return pd.DataFrame(
        {
            "peptide": [
                "MAGEAPEP1",
                "MAGEAPEP2",
                "PRAMEPEP1",
                "PRAMEPEP2",
                "NYESPEPT1",
                "NYESPEPT2",
                "NYESPEPT1",
            ],
            "length": [9, 9, 9, 9, 9, 9, 9],
            "gene_name": [
                "MAGEA4",
                "MAGEA4",
                "PRAME",
                "PRAME",
                "CTAG1A",
                "CTAG1B",
                "CTAG1B",
            ],
            "gene_id": ["E1", "E1", "E2", "E2", "E3A", "E3B", "E3B"],
        }
    )


def _stub_scores(peptides, alleles, **kwargs) -> pd.DataFrame:
    """Score every (peptide, allele) pair with deterministic percentiles
    derived from the peptide and allele strings — lets tests assert which
    peptide should win each cell."""
    rows = []
    for pep in peptides:
        for allele in alleles:
            # Lowest percentile = most-presented; PEP1 always beats PEP2
            # for the same gene; A*02:01 beats A*24:02 for all peptides.
            base = 0.1 if pep.endswith("1") else 1.5
            allele_bump = 0.0 if allele == "HLA-A*02:01" else 0.5
            rows.append(
                {
                    "peptide": pep,
                    "allele": allele,
                    "presentation_percentile": base + allele_bump,
                    "presentation_score": 0.95 - (base + allele_bump) / 10,
                    "affinity_nm": 50.0 + (base + allele_bump) * 100,
                    "affinity_percentile": (base + allele_bump) / 2,
                }
            )
    return pd.DataFrame(rows)


def _stub_ms_hits(peptides, **kwargs) -> pd.DataFrame:
    """Default test evidence: every peptide has unrestricted MS support."""
    return pd.DataFrame(
        {
            "peptide": list(peptides),
            "mhc_restriction": ["HLA class I"] * len(peptides),
            "mhc_allele_provenance": ["unmatched"] * len(peptides),
            "mhc_allele_set": [""] * len(peptides),
            "is_monoallelic": [False] * len(peptides),
            "src_cancer": [True] * len(peptides),
            "src_healthy_tissue": [False] * len(peptides),
            "src_healthy_reproductive": [False] * len(peptides),
            "src_healthy_thymus": [False] * len(peptides),
            "src_ebv_lcl": [False] * len(peptides),
            "pmid": ["1"] * len(peptides),
            "cell_line_name": [""] * len(peptides),
            "cell_name": [""] * len(peptides),
            "source_tissue": [""] * len(peptides),
        }
    )


def _stub_gene_ms_evidence(**kwargs) -> pd.DataFrame:
    """Gene-level healthy-MS rows used by automatic CTA safety selection."""
    return pd.DataFrame(
        {
            "peptide": ["MAGEA1UNIQ", "MAGEFAMSHARED", "NONVITAL"],
            "gene_names": ["MAGEA1", "MAGEA4;MAGEA8;MAGEA1", "PRAME"],
            "source_tissue": ["Central nervous system (CNS)", "Heart", "Blood"],
            "src_healthy_tissue": [True, True, True],
            "is_binding_assay": [False, False, False],
        }
    )


@pytest.fixture(autouse=True)
def _stub_pipeline(monkeypatch):
    """Wire stubs for every upstream call spanning_pmhc_set makes."""
    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides",
        lambda **kw: _stub_peptides(),
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.peptides.cta_peptides",
        lambda **kw: _stub_peptides(),
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.scoring.score_presentation",
        _stub_scores,
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.ms_evidence.load_public_ms_hits",
        _stub_ms_hits,
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.indexing.load_ms_evidence",
        _stub_gene_ms_evidence,
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.gene_sets.CTA_gene_names",
        lambda: {
            "MAGEA4",
            "PRAME",
            "CTAG1A",
            "CTAG1B",
            "MAGEA1",
            "MAGEB2",
            "VITALRNA",
            "BAGE",
        },
        raising=True,
    )
    # Stub the CTA dataframe so cta_count / rank_by paths work.
    cta_csv = pd.DataFrame(
        {
            "Symbol": [
                "MAGEA1",
                "MAGEB2",
                "MAGEA4",
                "PRAME",
                "VITALRNA",
                "CTAG1A",
                "CTAG1B",
                "BAGE",
            ],
            "Ensembl_Gene_ID": ["E0", "E0B", "E1", "E2", "E5", "E3A", "E3B", "E4"],
            "filtered": ["true"] * 8,
            "never_expressed": [
                "false",
                "false",
                "false",
                "false",
                "false",
                "false",
                "false",
                "true",
            ],
            "restriction_confidence": [
                "HIGH",
                "HIGH",
                "HIGH",
                "LOW",
                "HIGH",
                "HIGH",
                "MODERATE",
                "LOW",
            ],
            "restriction": ["TESTIS"] * 8,
            "ms_cancer_peptide_count": [60, 55, 50, 30, 20, 10, 10, 0],
            "ms_cta_exclusive_cancer_peptide_count": [60, 55, 50, 30, 20, 10, 10, 0],
            "rna_brain_max_ntpm": [0, 0, 0, 0.2, 1.5, 0.1, 0.7, 0],
            "rna_heart_max_ntpm": [0, 0, 0, 0.3, 0, 0, 0, 0],
            "rna_lung_max_ntpm": [0, 0, 0, 0.2, 0, 0, 0, 0],
            "rna_liver_max_ntpm": [0, 0, 0, 0.1, 0, 0, 0, 0],
            "rna_pancreas_max_ntpm": [0, 0, 0, 0.1, 0, 0, 0, 0],
            "ms_healthy_somatic_tissues": ["heart", "", "heart", "blood", "", "", "", ""],
        }
    )
    cancer_features = pd.DataFrame(
        {
            "gene_id": ["E0", "E0B", "E1", "E2", "E5", "E3A", "E3B", "E4"],
            "symbol": [
                "MAGEA1",
                "MAGEB2",
                "MAGEA4",
                "PRAME",
                "VITALRNA",
                "CTAG1A",
                "CTAG1B",
                "BAGE",
            ],
            "tumor_prevalence_panel_score": [100, 90, 80, 70, 60, 50, 50, 40],
            "cancer_rna_prevalent_cancer_type_count": [20, 19, 18, 17, 16, 15, 15, 14],
            "cancer_rna_sample_prevalence": [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.4, 0.3],
        }
    )
    monkeypatch.setattr("tsarina.loader.cta_dataframe", lambda: cta_csv, raising=True)
    monkeypatch.setattr(
        "tsarina.cancer_expression.cta_cancer_expression_features",
        lambda **kw: cancer_features,
        raising=True,
    )


# ── CTA selection ──────────────────────────────────────────────────────


def test_top_n_ranking_by_default_column():
    """MAGEA1 has the highest raw count but is dropped by the vital-tissue gate.
    PRAME is LOW-confidence but retained by the default clinical allowlist."""
    df = spanning_pmhc_set(
        cta_count=2,
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
    )
    assert list(df["cta"]) == ["MAGEA4", "PRAME"]


def test_selection_allowlist_is_pinned_into_automatic_top_n():
    ctas = _resolve_ctas(
        ctas=None,
        cta_count=3,
        cta_rank_by="ms_cancer_peptide_count",
        min_restriction_confidence=("HIGH", "MODERATE"),
        restriction_levels=None,
        selection_allowlist=["PRAME", "CTAG1A/CTAG1B", "MAGEA4"],
    )
    assert ctas == ["MAGEA4", "PRAME", "CTAG1A/CTAG1B"]


def test_removed_packaged_ms_rank_column_raises(monkeypatch):
    cta_csv = pd.DataFrame(
        {
            "Symbol": ["BROADONLY", "EXCLUSIVE"],
            "Ensembl_Gene_ID": ["E7", "E8"],
            "filtered": ["true", "true"],
            "never_expressed": ["false", "false"],
            "restriction": ["TESTIS", "TESTIS"],
        }
    )
    monkeypatch.setattr("tsarina.loader.cta_dataframe", lambda: cta_csv, raising=True)
    monkeypatch.setattr(
        "tsarina.gene_sets.CTA_gene_names",
        lambda: {"BROADONLY", "EXCLUSIVE"},
        raising=True,
    )

    with pytest.raises(ValueError, match="runtime-derived public-MS column"):
        _resolve_ctas(
            ctas=None,
            cta_count=2,
            cta_rank_by="ms_cta_exclusive_cancer_peptide_count",
            min_restriction_confidence=None,
            restriction_levels=None,
            selection_allowlist=[],
            exclude_vital_tissue_expression=False,
            exclude_non_magea4_mage_family=False,
        )


def test_default_rank_column_prefers_tumor_prevalence(monkeypatch):
    cta_csv = pd.DataFrame(
        {
            "Symbol": ["MSHIGH", "TUMORHIGH"],
            "Ensembl_Gene_ID": ["E7", "E8"],
            "filtered": ["true", "true"],
            "never_expressed": ["false", "false"],
            "restriction_confidence": ["HIGH", "HIGH"],
            "restriction": ["TESTIS", "TESTIS"],
            "ms_cta_exclusive_cancer_peptide_count": [100, 1],
            "ms_healthy_somatic_tissues": ["", ""],
        }
    )
    cancer_features = pd.DataFrame(
        {
            "gene_id": ["E7", "E8"],
            "symbol": ["MSHIGH", "TUMORHIGH"],
            "tumor_prevalence_panel_score": [1.0, 5.0],
        }
    )
    monkeypatch.setattr("tsarina.loader.cta_dataframe", lambda: cta_csv, raising=True)
    monkeypatch.setattr(
        "tsarina.gene_sets.CTA_gene_names",
        lambda: {"MSHIGH", "TUMORHIGH"},
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.cancer_expression.cta_cancer_expression_features",
        lambda **kw: cancer_features,
        raising=True,
    )

    ctas = _resolve_ctas(
        ctas=None,
        cta_count=2,
        cta_rank_by=_DEFAULT_RANK_COLUMN,
        min_restriction_confidence=("HIGH", "MODERATE"),
        restriction_levels=None,
        selection_allowlist=[],
        exclude_vital_tissue_expression=False,
        exclude_non_magea4_mage_family=False,
    )

    assert ctas == ["TUMORHIGH", "MSHIGH"]


def test_live_ms_gate_drops_zero_ms_candidate_after_prevalence_rank(monkeypatch):
    cta_csv = pd.DataFrame(
        {
            "Symbol": ["MSLOW", "ZEROMSHIGH"],
            "Ensembl_Gene_ID": ["E7", "E8"],
            "filtered": ["true", "true"],
            "never_expressed": ["false", "false"],
            "restriction": ["TESTIS", "TESTIS"],
        }
    )
    cancer_features = pd.DataFrame(
        {
            "gene_id": ["E7", "E8"],
            "symbol": ["MSLOW", "ZEROMSHIGH"],
            "tumor_prevalence_panel_score": [1.0, 10.0],
        }
    )
    monkeypatch.setattr("tsarina.loader.cta_dataframe", lambda: cta_csv, raising=True)
    monkeypatch.setattr(
        "tsarina.gene_sets.CTA_gene_names",
        lambda: {"MSLOW", "ZEROMSHIGH"},
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.cancer_expression.cta_cancer_expression_features",
        lambda **kw: cancer_features,
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides",
        lambda **kw: pd.DataFrame(
            {
                "peptide": ["MSLOWPEP", "ZEROPEP"],
                "length": [9, 9],
                "gene_name": ["MSLOW", "ZEROMSHIGH"],
                "gene_id": ["E7", "E8"],
            }
        ),
        raising=True,
    )

    def _live_hits(peptides, **kwargs):
        return pd.DataFrame(
            {
                "peptide": ["MSLOWPEP"],
                "mhc_restriction": ["HLA class I"],
                "mhc_allele_provenance": ["unmatched"],
                "mhc_allele_set": [""],
                "is_monoallelic": [False],
                "src_cancer": [True],
                "src_healthy_tissue": [False],
                "src_healthy_reproductive": [False],
                "src_healthy_thymus": [False],
                "pmid": ["1"],
            }
        )

    monkeypatch.setattr("tsarina.ms_evidence.load_public_ms_hits", _live_hits, raising=True)

    df = spanning_pmhc_set(
        cta_count=1,
        alleles=["HLA-A*02:01"],
        cta_rank_by=_DEFAULT_RANK_COLUMN,
        selection_allowlist=[],
        exclude_vital_tissue_expression=False,
        exclude_non_magea4_mage_family=False,
        max_percentile=10.0,
    )

    assert list(df["cta"]) == ["MSLOW"]
    assert df.attrs["input_cta_order"] == ["ZEROMSHIGH", "MSLOW"]
    assert df.attrs["empty_ctas"] == ["ZEROMSHIGH"]


def test_panel_requires_live_hitlist_ms_evidence(monkeypatch):
    def _missing_hitlist(peptides, **kwargs):
        raise FileNotFoundError("missing observations.parquet")

    monkeypatch.setattr("tsarina.ms_evidence.load_public_ms_hits", _missing_hitlist, raising=True)

    with pytest.raises(RuntimeError, match="hitlist observations index"):
        spanning_pmhc_set(
            cta_count=1,
            alleles=["HLA-A*02:01"],
            max_percentile=10.0,
        )


def test_automatic_selection_hides_empty_ctas_by_default():
    df = spanning_pmhc_set(
        cta_count=4,
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
    )
    assert list(df["cta"]) == ["MAGEA4", "PRAME", "CTAG1A/CTAG1B"]
    assert df.attrs["input_cta_order"] == [
        "MAGEA4",
        "PRAME",
        "CTAG1A/CTAG1B",
        "VITALRNA",
    ]
    assert df.attrs["empty_ctas"] == ["VITALRNA"]
    assert df.attrs["panel_summary"]["empty_cta_count"] == 1


def test_automatic_selection_backfills_empty_ctas_by_default(monkeypatch):
    all_peptides = pd.concat(
        [
            _stub_peptides(),
            pd.DataFrame(
                {
                    "peptide": ["BACKFILLPEP1"],
                    "length": [9],
                    "gene_name": ["BACKFILL"],
                    "gene_id": ["E6"],
                }
            ),
        ],
        ignore_index=True,
    )
    cta_csv = pd.DataFrame(
        {
            "Symbol": ["MAGEA4", "PRAME", "CTAG1A", "CTAG1B", "VITALRNA", "BACKFILL"],
            "Ensembl_Gene_ID": ["E1", "E2", "E3A", "E3B", "E5", "E6"],
            "filtered": ["true"] * 6,
            "never_expressed": ["false"] * 6,
            "restriction_confidence": ["HIGH"] * 6,
            "restriction": ["TESTIS"] * 6,
            "ms_cancer_peptide_count": [50, 30, 10, 10, 20, 15],
            "ms_cta_exclusive_cancer_peptide_count": [50, 30, 10, 10, 20, 15],
            "rna_brain_max_ntpm": [0, 0, 0, 0, 1.5, 0],
            "rna_heart_max_ntpm": [0] * 6,
            "rna_lung_max_ntpm": [0] * 6,
            "rna_liver_max_ntpm": [0] * 6,
            "rna_pancreas_max_ntpm": [0] * 6,
            "ms_healthy_somatic_tissues": [""] * 6,
        }
    )
    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides",
        lambda **kw: all_peptides,
        raising=True,
    )
    monkeypatch.setattr("tsarina.loader.cta_dataframe", lambda: cta_csv, raising=True)
    monkeypatch.setattr(
        "tsarina.gene_sets.CTA_gene_names",
        lambda: {"MAGEA4", "PRAME", "CTAG1A", "CTAG1B", "VITALRNA", "BACKFILL"},
        raising=True,
    )

    df = spanning_pmhc_set(
        cta_count=4,
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
    )

    assert list(df["cta"]) == ["MAGEA4", "PRAME", "CTAG1A/CTAG1B", "BACKFILL"]
    assert df.attrs["input_cta_order"] == [
        "MAGEA4",
        "PRAME",
        "CTAG1A/CTAG1B",
        "VITALRNA",
        "BACKFILL",
    ]
    assert df.attrs["empty_ctas"] == ["VITALRNA"]


def test_include_empty_ctas_preserves_automatic_failures():
    df = spanning_pmhc_set(
        cta_count=4,
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
        include_empty_ctas=True,
    )
    assert list(df["cta"]) == ["MAGEA4", "PRAME", "CTAG1A/CTAG1B", "VITALRNA"]


def test_explicit_ctas_overrides_ranking():
    """When --ctas is supplied, --cta-count is ignored and the order is
    preserved after alias/group normalization."""
    df = spanning_pmhc_set(
        ctas=["CTAG1B", "PRAME"],
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
    )
    assert list(df["cta"]) == ["CTAG1A/CTAG1B", "PRAME"]


def test_explicit_ctas_accepts_clinical_aliases():
    df = spanning_pmhc_set(
        ctas=["MAGE-A4", "NY-ESO-1"],
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
        peptides_per_cell=1,
    )
    assert list(df["cta"]) == ["MAGEA4", "CTAG1A/CTAG1B"]
    assert df.set_index("cta").loc["CTAG1A/CTAG1B", "HLA-A*02:01"] == "NYESPEPT1"


def test_explicit_ctas_bypass_mage_family_gate():
    df = spanning_pmhc_set(
        ctas=["MAGEA1"],
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
    )
    assert list(df["cta"]) == ["MAGEA1"]


def test_ctag1_group_expands_ctag1a_and_ctag1b_for_peptide_resolution(monkeypatch):
    calls: list[list[str]] = []

    def _exclusive(**kw):
        calls.append(kw["gene_names"])
        return _stub_peptides()

    monkeypatch.setattr("tsarina.peptides.cta_exclusive_peptides", _exclusive, raising=True)

    df = spanning_pmhc_set(
        ctas=["NY-ESO-1"],
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
    )

    assert calls == [["CTAG1A", "CTAG1B"]]
    assert list(df["cta"]) == ["CTAG1A/CTAG1B"]


def test_min_restriction_confidence_filters_low():
    """LOW-confidence CTAs drop out unless allowlisted."""
    df = spanning_pmhc_set(
        cta_count=10,
        min_restriction_confidence=("HIGH", "MODERATE"),
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
    )
    assert "BAGE" not in df["cta"].tolist()
    assert "PRAME" in df["cta"].tolist()


def test_min_restriction_confidence_none_admits_low():
    """Disabling the confidence gate does not disable the vital-tissue gate
    or never-expressed filtering."""
    df = spanning_pmhc_set(
        cta_count=10,
        min_restriction_confidence=None,
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
        include_empty_ctas=True,
    )
    assert "BAGE" not in df["cta"].tolist()
    assert "MAGEA1" not in df["cta"].tolist()
    assert "MAGEB2" not in df["cta"].tolist()
    assert set(df["cta"]) == {"MAGEA4", "PRAME", "VITALRNA", "CTAG1A/CTAG1B"}


def test_restriction_levels_filter():
    """restriction_levels=('TESTIS',) keeps the grouped CTAG1 target."""
    df = spanning_pmhc_set(
        cta_count=10,
        restriction_levels=("TESTIS",),
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
    )
    assert "CTAG1A/CTAG1B" in df["cta"].tolist()


def test_vital_tissue_gate_can_be_disabled():
    ctas = _resolve_ctas(
        ctas=None,
        cta_count=1,
        cta_rank_by="ms_cancer_peptide_count",
        min_restriction_confidence=("HIGH", "MODERATE"),
        restriction_levels=None,
        selection_allowlist=[],
        exclude_vital_tissue_expression=False,
        exclude_non_magea4_mage_family=False,
    )
    assert ctas == ["MAGEA1"]


def test_default_mage_family_gate_allows_only_magea4():
    ctas = _resolve_ctas(
        ctas=None,
        cta_count=10,
        cta_rank_by="ms_cancer_peptide_count",
        min_restriction_confidence=None,
        restriction_levels=None,
        selection_allowlist=[],
        exclude_vital_tissue_expression=False,
    )
    assert "MAGEA4" in ctas
    assert "MAGEA1" not in ctas
    assert "MAGEB2" not in ctas


def test_mage_family_gate_allows_selection_allowlist():
    ctas = _resolve_ctas(
        ctas=None,
        cta_count=10,
        cta_rank_by="ms_cancer_peptide_count",
        min_restriction_confidence=None,
        restriction_levels=None,
        selection_allowlist=["MAGEB2"],
        exclude_vital_tissue_expression=False,
    )
    assert "MAGEB2" in ctas
    assert "MAGEA1" not in ctas


def test_default_vital_rna_gate_allows_sub_2_ntpm():
    ctas = _resolve_ctas(
        ctas=None,
        cta_count=10,
        cta_rank_by="ms_cancer_peptide_count",
        min_restriction_confidence=("HIGH", "MODERATE"),
        restriction_levels=None,
        selection_allowlist=[],
    )
    assert "VITALRNA" in ctas


def test_unique_vital_ms_gate_ignores_shared_family_peptide_evidence():
    ctas = _resolve_ctas(
        ctas=None,
        cta_count=10,
        cta_rank_by="ms_cancer_peptide_count",
        min_restriction_confidence=("HIGH", "MODERATE"),
        restriction_levels=None,
        selection_allowlist=[],
    )
    assert "MAGEA4" in ctas
    assert "MAGEA1" not in ctas


def test_vital_rna_gate_threshold_is_parameterizable():
    ctas = _resolve_ctas(
        ctas=None,
        cta_count=10,
        cta_rank_by="ms_cancer_peptide_count",
        min_restriction_confidence=("HIGH", "MODERATE"),
        restriction_levels=None,
        selection_allowlist=[],
        vital_tissue_max_ntpm=1.0,
    )
    assert "VITALRNA" not in ctas


def test_panel_summary_sorts_ctas_by_selected_peptides():
    selected = pd.DataFrame(
        {
            "cta": ["LOW", "HIGH", "HIGH"],
            "allele": ["HLA-A*02:01", "HLA-A*02:01", "HLA-A*24:02"],
            "peptide": ["LOWPEP", "HIGHPEP1", "HIGHPEP2"],
            "evidence_tier": ["unrestricted_ms", "unrestricted_ms", "unrestricted_ms"],
        }
    )
    summary = panel_summary(
        selected=selected,
        cta_list=["ZERO", "LOW", "HIGH"],
        allele_list=["HLA-A*02:01", "HLA-A*24:02"],
        allele_frequencies={"HLA-A*02:01": 0.2, "HLA-A*24:02": 0.1},
    )
    assert [row["cta"] for row in summary["cta_coverage"]] == ["HIGH", "LOW", "ZERO"]
    by_cta = {row["cta"]: row for row in summary["cta_coverage"]}
    assert by_cta["HIGH"]["estimated_population_coverage"] == pytest.approx(0.51)
    assert by_cta["LOW"]["estimated_population_coverage"] == pytest.approx(0.36)


def test_panel_summary_counts_ms_tiers_per_cta():
    selected = pd.DataFrame(
        {
            "cta": ["MAGEA4", "MAGEA4", "MAGEA4", "PRAME"],
            "allele": ["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02", "HLA-A*02:01"],
            "peptide": ["PEP1", "PEP2", "PEP3", "PEP4"],
            "evidence_tier": [
                "monoallelic_ms",
                "sample_allele_ms",
                "unrestricted_ms",
                "sample_allele_ms",
            ],
        }
    )
    summary = panel_summary(
        selected=selected,
        cta_list=["MAGEA4", "PRAME"],
        allele_list=["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02"],
        allele_frequencies={"HLA-A*02:01": 0.2, "HLA-A*24:02": 0.1, "HLA-B*07:02": 0.1},
    )
    by_cta = {row["cta"]: row for row in summary["cta_coverage"]}

    assert by_cta["MAGEA4"]["monoallelic_ms_pmhc_count"] == 1
    assert by_cta["MAGEA4"]["sample_allele_ms_pmhc_count"] == 1
    assert by_cta["MAGEA4"]["unrestricted_ms_pmhc_count"] == 1
    assert by_cta["MAGEA4"]["estimated_population_coverage"] == pytest.approx(0.6031)
    assert by_cta["PRAME"]["monoallelic_ms_pmhc_count"] == 0
    assert by_cta["PRAME"]["sample_allele_ms_pmhc_count"] == 1


# ── Allele resolution ──────────────────────────────────────────────────


def test_panel_default_resolves_via_get_panel():
    """Default 53-allele panel should produce 54 columns including 'cta'."""
    df = spanning_pmhc_set(
        cta_count=2,
        max_percentile=10.0,
    )
    assert df.shape[1] == 54  # 1 cta + 53 alleles


def test_explicit_alleles_override_panel():
    df = spanning_pmhc_set(
        cta_count=2,
        alleles=["HLA-A*02:01", "HLA-A*24:02"],
        max_percentile=10.0,
    )
    assert list(df.columns) == ["cta", "HLA-A*02:01", "HLA-A*24:02"]


def test_no_alleles_no_panel_raises():
    with pytest.raises(ValueError, match="alleles or panel"):
        spanning_pmhc_set(panel=None, alleles=None)


# ── Best-per-cell selection ────────────────────────────────────────────


def test_lowest_percentile_peptide_wins_each_cell():
    """Stub gives PEP1 = 0.1, PEP2 = 1.5 baseline. Both pass max=2.0;
    PEP1 should win every cell."""
    df = spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01"],
        max_percentile=2.0,
        peptides_per_cell=1,
    )
    assert df.set_index("cta").loc["MAGEA4", "HLA-A*02:01"] == "MAGEAPEP1"
    assert df.set_index("cta").loc["PRAME", "HLA-A*02:01"] == "PRAMEPEP1"


def test_max_percentile_filters_weak_cells():
    """Stub: A*24:02 percentiles are 0.6 (PEP1) and 2.0 (PEP2). A cutoff
    of 0.5 leaves every A*24:02 cell empty; A*02:01 cells (0.1) survive."""
    df = spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01", "HLA-A*24:02"],
        max_percentile=0.5,
        peptides_per_cell=1,
    )
    indexed = df.set_index("cta")
    assert indexed.loc["MAGEA4", "HLA-A*02:01"] == "MAGEAPEP1"
    assert pd.isna(indexed.loc["MAGEA4", "HLA-A*24:02"])


def test_monoallelic_ms_tier_and_threshold(monkeypatch):
    hits = pd.DataFrame(
        {
            "peptide": ["MAGEAPEP1"],
            "mhc_restriction": ["HLA-A*02:01"],
            "mhc_allele_provenance": ["exact"],
            "mhc_allele_set": ["HLA-A*02:01"],
            "is_monoallelic": [True],
            "pmid": ["111"],
            "cell_line_name": ["K562-A*02:01"],
        }
    )
    monkeypatch.setattr("tsarina.ms_evidence.load_public_ms_hits", lambda peptides, **kw: hits)

    long = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01"],
        monoallelic_ms_max_percentile=0.2,
        output_format="long",
    )
    assert len(long) == 1
    row = long.iloc[0]
    assert row["evidence_tier"] == "monoallelic_ms"
    assert row["ms_hit_count"] == 1
    assert row["ms_pmids"] == "111"

    blocked = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01"],
        monoallelic_ms_max_percentile=0.05,
        output_format="long",
    )
    assert blocked.empty


def test_sample_allele_ms_requires_best_among_sample_alleles(monkeypatch):
    hits = pd.DataFrame(
        {
            "peptide": ["MAGEAPEP1"],
            "mhc_restriction": ["HLA class I"],
            "mhc_allele_provenance": ["sample_allele_match"],
            "mhc_allele_set": ["HLA-A*02:01;HLA-B*07:02"],
            "is_monoallelic": [False],
            "pmid": ["222"],
            "cell_line_name": ["tumor sample"],
        }
    )
    monkeypatch.setattr("tsarina.ms_evidence.load_public_ms_hits", lambda peptides, **kw: hits)

    long = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01", "HLA-B*07:02"],
        sample_allele_ms_max_percentile=1.0,
        output_format="long",
    )
    assert len(long) == 1
    row = long.iloc[0]
    assert row["allele"] == "HLA-A*02:01"
    assert row["evidence_tier"] == "sample_allele_ms"
    assert row["ms_pmids"] == "222"


def test_peptide_attribution_provenance_treated_like_sample_allele_match(monkeypatch):
    """Regression for hitlist#45 / pirl-unc/tsarina#62.  hitlist v1.30.39
    introduced a new ``mhc_allele_provenance`` value, ``peptide_attribution``,
    emitted for class-only IEDB rows whose peptide appears in a paper
    supplement (Sarkizova 2020 patient cohort etc.) — the candidate set
    is narrowed to the matched donor's typed alleles.  Same shape as
    ``sample_allele_match`` but strictly more specific.

    Spanning's ``_attribute_panel_evidence_with_observations`` previously
    filtered ``provenance == "sample_allele_match"`` and silently skipped
    these rows, dropping ~36 K Sarkizova patient-tumor MS observations
    from panel-evidence attribution.  After the fix the same row gets
    treated as ``sample_allele_ms`` evidence."""
    hits = pd.DataFrame(
        {
            "peptide": ["MAGEAPEP1"],
            "mhc_restriction": ["HLA-A*02:01;HLA-B*07:02"],
            "mhc_allele_provenance": ["peptide_attribution"],  # NEW provenance
            "mhc_allele_set": ["HLA-A*02:01;HLA-B*07:02"],
            "is_monoallelic": [False],
            "pmid": ["31844290"],  # Sarkizova
            "cell_line_name": ["MEL2"],
        }
    )
    monkeypatch.setattr("tsarina.ms_evidence.load_public_ms_hits", lambda peptides, **kw: hits)

    long = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01", "HLA-B*07:02"],
        sample_allele_ms_max_percentile=1.0,
        output_format="long",
    )
    assert len(long) == 1
    row = long.iloc[0]
    assert row["allele"] == "HLA-A*02:01"
    assert row["evidence_tier"] == "sample_allele_ms"
    assert row["ms_pmids"] == "31844290"


def test_extend_alleles_picks_up_peptide_attribution_donor_set(monkeypatch):
    """``_extend_alleles_from_observations`` walks observation rows and
    extends the panel allele list with the donor's typed alleles when
    available.  Pre-tsarina#62 it filtered on ``sample_allele_match``
    only; after the fix ``peptide_attribution`` rows are also visited.

    This is a ``best-of-sample`` correctness gate — without the donor
    alleles in the working list, scoring functions can't pick the best
    haplotype for the peptide even when the data is on the row."""
    hits = pd.DataFrame(
        {
            "peptide": ["MAGEAPEP1"],
            "mhc_restriction": ["HLA-A*02:01;HLA-B*07:02;HLA-C*06:02"],
            "mhc_allele_provenance": ["peptide_attribution"],
            "mhc_allele_set": ["HLA-A*02:01;HLA-B*07:02;HLA-C*06:02"],
            "is_monoallelic": [False],
            "pmid": ["31844290"],
            "cell_line_name": ["MEL2"],
        }
    )
    monkeypatch.setattr("tsarina.ms_evidence.load_public_ms_hits", lambda peptides, **kw: hits)

    # Panel only includes HLA-A*02:01.  After best-of-sample extension
    # the working allele list should also include the row's typed
    # B*07:02 and C*06:02 — and the picked allele for the row should
    # win out as the best haplotype across all three.
    long = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01"],
        sample_allele_ms_max_percentile=1.0,
        output_format="long",
    )
    assert len(long) == 1
    row = long.iloc[0]
    assert row["evidence_tier"] == "sample_allele_ms"
    # The picked allele must come from the panel intersect the
    # row's sample alleles.
    assert row["allele"] == "HLA-A*02:01"


def test_sample_allele_ms_exact_sample_restrictions_still_use_best_haplotype(monkeypatch):
    hits = pd.DataFrame(
        {
            "peptide": ["MAGEAPEP1"],
            "mhc_restriction": ["HLA-A*02:01;HLA-B*07:02"],
            "mhc_allele_provenance": ["sample_allele_match"],
            "mhc_allele_set": ["HLA-A*02:01;HLA-B*07:02"],
            "is_monoallelic": [False],
            "pmid": ["222"],
            "cell_line_name": ["tumor sample"],
        }
    )
    monkeypatch.setattr("tsarina.ms_evidence.load_public_ms_hits", lambda peptides, **kw: hits)

    long = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01", "HLA-B*07:02"],
        sample_allele_ms_max_percentile=1.0,
        output_format="long",
    )
    assert len(long) == 1
    row = long.iloc[0]
    assert row["allele"] == "HLA-A*02:01"
    assert row["evidence_tier"] == "sample_allele_ms"


def test_unrestricted_ms_tier_uses_processing_evidence(monkeypatch):
    hits = pd.DataFrame(
        {
            "peptide": ["MAGEAPEP1"],
            "mhc_restriction": ["HLA class I"],
            "mhc_allele_provenance": ["unmatched"],
            "mhc_allele_set": [""],
            "is_monoallelic": [False],
            "pmid": ["333"],
            "source_tissue": ["tumor"],
        }
    )
    monkeypatch.setattr("tsarina.ms_evidence.load_public_ms_hits", lambda peptides, **kw: hits)

    long = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01", "HLA-A*24:02"],
        unrestricted_ms_max_percentile=0.5,
        output_format="long",
    )
    assert len(long) == 1
    row = long.iloc[0]
    assert row["allele"] == "HLA-A*02:01"
    assert row["evidence_tier"] == "unrestricted_ms"
    assert row["ms_samples"] == "tumor"


def test_predicted_only_is_optional_and_last_priority(monkeypatch):
    empty_hits = pd.DataFrame(
        columns=[
            "peptide",
            "mhc_restriction",
            "mhc_allele_provenance",
            "mhc_allele_set",
            "is_monoallelic",
        ]
    )
    monkeypatch.setattr(
        "tsarina.ms_evidence.load_public_ms_hits", lambda peptides, **kw: empty_hits
    )

    default = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01"],
        output_format="long",
    )
    assert default.empty

    predicted = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01"],
        include_predicted_only=True,
        predicted_only_max_percentile=0.2,
        output_format="long",
    )
    assert len(predicted) == 1
    assert predicted.iloc[0]["evidence_tier"] == "predicted_only"

    ms_hits = pd.DataFrame(
        {
            "peptide": ["MAGEAPEP2"],
            "mhc_restriction": ["HLA-A*02:01"],
            "mhc_allele_provenance": ["exact"],
            "mhc_allele_set": ["HLA-A*02:01"],
            "is_monoallelic": [True],
            "pmid": ["444"],
        }
    )
    monkeypatch.setattr("tsarina.ms_evidence.load_public_ms_hits", lambda peptides, **kw: ms_hits)
    ms_first = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01"],
        include_predicted_only=True,
        predicted_only_max_percentile=0.2,
        monoallelic_ms_max_percentile=2.0,
        output_format="long",
    )
    assert len(ms_first) == 1
    assert ms_first.iloc[0]["peptide"] == "MAGEAPEP2"
    assert ms_first.iloc[0]["evidence_tier"] == "monoallelic_ms"


def test_peptides_per_cell_keeps_ranked_top_n():
    long = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
        peptides_per_cell=2,
        output_format="long",
    )
    assert list(long["peptide"]) == ["MAGEAPEP1", "MAGEAPEP2"]
    assert list(long["peptide_rank_in_cell"]) == [1, 2]

    wide = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
        peptides_per_cell=2,
    )
    assert wide.set_index("cta").loc["MAGEA4", "HLA-A*02:01"] == "MAGEAPEP1; MAGEAPEP2"


def test_peptides_per_cell_ranks_ms_support_before_prediction(monkeypatch):
    hits = pd.DataFrame(
        {
            "peptide": ["MAGEAPEP1", "MAGEAPEP2", "MAGEAPEP2"],
            "mhc_restriction": ["HLA class I", "HLA class I", "HLA class I"],
            "mhc_allele_provenance": ["unmatched", "unmatched", "unmatched"],
            "mhc_allele_set": ["", "", ""],
            "is_monoallelic": [False, False, False],
            "pmid": ["111", "222", "333"],
        }
    )
    monkeypatch.setattr("tsarina.ms_evidence.load_public_ms_hits", lambda peptides, **kw: hits)

    long = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
        peptides_per_cell=2,
        output_format="long",
    )
    assert list(long["peptide"]) == ["MAGEAPEP2", "MAGEAPEP1"]
    assert list(long["ms_hit_count"]) == [2, 1]
    assert list(long["ms_source_count"]) == [2, 1]


# ── Output formats ─────────────────────────────────────────────────────


def test_long_format_has_one_row_per_filled_cell():
    long = spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01", "HLA-A*24:02"],
        max_percentile=2.0,
        peptides_per_cell=1,
        output_format="long",
    )
    # 2 CTAs * 2 alleles = 4 cells, all pass with max=2.0
    assert len(long) == 4
    assert set(long.columns) == {
        "cta",
        "cta_members",
        "allele",
        "peptide",
        "length",
        "evidence_tier",
        "ms_hit_count",
        "ms_source_count",
        "ms_alleles",
        "ms_pmids",
        "ms_samples",
        "presentation_percentile",
        "presentation_score",
        "affinity_nm",
        "affinity_percentile",
        "netmhcpan_affinity_nm",
        "netmhcpan_affinity_percentile",
        "peptide_rank_in_cell",
    }
    # Sorted by cta then allele in the supplied order
    assert list(long["cta"]) == ["MAGEA4", "MAGEA4", "PRAME", "PRAME"]


def test_long_format_drops_filtered_cells_entirely():
    """In long output, cells above max_percentile are absent (not NaN)."""
    long = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01", "HLA-A*24:02"],
        max_percentile=0.5,
        peptides_per_cell=1,
        output_format="long",
    )
    assert len(long) == 1
    assert long.iloc[0]["allele"] == "HLA-A*02:01"


def test_unrecognized_output_format_raises():
    with pytest.raises(ValueError, match="output_format"):
        spanning_pmhc_set(output_format="grid")


def test_output_attrs_include_summary_and_order_metadata():
    long = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01", "HLA-A*24:02"],
        max_percentile=10.0,
        peptides_per_cell=2,
        output_format="long",
    )
    summary = long.attrs["panel_summary"]
    assert long.attrs["cta_order"] == ["MAGEA4"]
    assert long.attrs["allele_order"] == ["HLA-A*02:01", "HLA-A*24:02"]
    assert summary["hla_allele_count"] == 2
    assert summary["cta_count"] == 1
    assert summary["filled_cell_count"] == 2
    assert summary["selected_peptide_count"] == 2
    assert summary["cta_coverage"][0]["covered_hla_count"] == 2
    assert summary["hla_coverage"][0]["covered_cta_fraction"] == 1.0


def test_identical_selected_pmhc_ctas_are_grouped(monkeypatch):
    peptides = pd.DataFrame(
        {
            "peptide": ["SHAREDPEP", "SHAREDPEP", "PRAMEPEP1"],
            "length": [9, 9, 9],
            "gene_name": ["SPANXD", "SPANXA1", "PRAME"],
            "gene_id": ["E10", "E11", "E2"],
        }
    )
    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides",
        lambda **kw: peptides,
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.gene_sets.CTA_gene_names",
        lambda: {"SPANXD", "SPANXA1", "PRAME"},
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.loader.cta_dataframe",
        lambda: pd.DataFrame(
            {
                "Symbol": ["SPANXD", "SPANXA1", "PRAME"],
                "ms_cta_exclusive_cancer_peptide_count": [3, 2, 1],
            }
        ),
        raising=True,
    )

    long = spanning_pmhc_set(
        ctas=["SPANXD", "SPANXA1", "PRAME"],
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
        peptides_per_cell=1,
        output_format="long",
    )

    assert long.attrs["cta_order"] == ["SPANXD/SPANXA1", "PRAME"]
    assert long.attrs["cta_groups"] == [
        {"cta": "SPANXD/SPANXA1", "members": ["SPANXD", "SPANXA1"]},
        {"cta": "PRAME", "members": ["PRAME"]},
    ]
    grouped = long[long["cta"] == "SPANXD/SPANXA1"].iloc[0]
    assert grouped["cta_members"] == "SPANXD;SPANXA1"
    assert grouped["peptide"] == "SHAREDPEP"
    assert long.attrs["panel_summary"]["grouped_cta_member_count"] == 1


def test_identical_peptide_set_ctas_are_grouped_by_combined_name(monkeypatch):
    peptides = pd.DataFrame(
        {
            "peptide": [
                "XAGEPEP1",
                "XAGEPEP2",
                "XAGEPEP1",
                "XAGEPEP2",
                "PAGE2PEP1",
            ],
            "length": [9, 9, 9, 9, 9],
            "gene_name": ["XAGE1A", "XAGE1A", "XAGE1B", "XAGE1B", "PAGE2"],
            "gene_id": ["E10", "E10", "E11", "E11", "E12"],
        }
    )
    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides",
        lambda **kw: peptides,
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.gene_sets.CTA_gene_names",
        lambda: {"XAGE1A", "XAGE1B", "PAGE2"},
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.loader.cta_dataframe",
        lambda: pd.DataFrame(
            {
                "Symbol": ["XAGE1A", "XAGE1B", "PAGE2"],
                "ms_cta_exclusive_cancer_peptide_count": [3, 2, 1],
            }
        ),
        raising=True,
    )
    messages: list[str] = []

    long = spanning_pmhc_set(
        ctas=["XAGE1A", "XAGE1B", "PAGE2"],
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
        peptides_per_cell=1,
        output_format="long",
        on_progress=messages.append,
    )

    assert long.attrs["cta_order"] == ["XAGE1A/XAGE1B", "PAGE2"]
    assert long.attrs["cta_groups"] == [
        {"cta": "XAGE1A/XAGE1B", "members": ["XAGE1A", "XAGE1B"]},
        {"cta": "PAGE2", "members": ["PAGE2"]},
    ]
    grouped = long[long["cta"] == "XAGE1A/XAGE1B"].iloc[0]
    assert grouped["cta_members"] == "XAGE1A;XAGE1B"
    assert grouped["peptide"] == "XAGEPEP1"
    assert any("XAGE1A/XAGE1B" in msg for msg in messages)


def test_automatic_backfill_counts_distinct_cta_pmhc_groups(monkeypatch):
    peptides = pd.DataFrame(
        {
            "peptide": ["SHAREDPEP", "SHAREDPEP", "DISTINCT1"],
            "length": [9, 9, 9],
            "gene_name": ["DUP1", "DUP2", "DISTINCT"],
            "gene_id": ["E10", "E11", "E12"],
        }
    )
    cta_csv = pd.DataFrame(
        {
            "Symbol": ["DUP1", "DUP2", "DISTINCT"],
            "filtered": ["true", "true", "true"],
            "never_expressed": ["false", "false", "false"],
            "restriction_confidence": ["HIGH", "HIGH", "HIGH"],
            "restriction": ["TESTIS", "TESTIS", "TESTIS"],
            "ms_cta_exclusive_cancer_peptide_count": [30, 20, 10],
            "ms_healthy_somatic_tissues": ["", "", ""],
        }
    )
    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides",
        lambda **kw: peptides,
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.gene_sets.CTA_gene_names",
        lambda: {"DUP1", "DUP2", "DISTINCT"},
        raising=True,
    )
    monkeypatch.setattr("tsarina.loader.cta_dataframe", lambda: cta_csv, raising=True)

    long = spanning_pmhc_set(
        cta_count=2,
        cta_rank_by="ms_cta_exclusive_cancer_peptide_count",
        alleles=["HLA-A*02:01"],
        selection_allowlist=(),
        exclude_vital_tissue_expression=False,
        max_percentile=10.0,
        peptides_per_cell=1,
        output_format="long",
    )

    assert long.attrs["cta_order"] == ["DUP1/DUP2", "DISTINCT"]
    assert set(long["cta"]) == {"DUP1/DUP2", "DISTINCT"}
    assert long.attrs["panel_summary"]["input_cta_count"] == 3


def test_netmhcpan_affinity_annotation_scores_selected_pmhcs(monkeypatch):
    calls: list[dict[str, object]] = []

    def _score_with_netmhcpan_annotation(peptides, alleles, predictor="mhcflurry", **kwargs):
        calls.append(
            {
                "predictor": predictor,
                "peptides": tuple(peptides),
                "alleles": tuple(alleles),
            }
        )
        rows = []
        for pep in peptides:
            for allele in alleles:
                is_netmhcpan = predictor == "netmhcpan"
                rows.append(
                    {
                        "peptide": pep,
                        "allele": allele,
                        "presentation_percentile": 0.1,
                        "presentation_score": 0.95,
                        "affinity_nm": 25.0 if is_netmhcpan else 150.0,
                        "affinity_percentile": 0.03 if is_netmhcpan else pd.NA,
                    }
                )
        return pd.DataFrame(rows)

    monkeypatch.setattr(
        "tsarina.scoring.score_presentation",
        _score_with_netmhcpan_annotation,
        raising=True,
    )

    long = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01"],
        max_percentile=10.0,
        peptides_per_cell=1,
        output_format="long",
        annotate_netmhcpan_affinity=True,
    )

    assert [call["predictor"] for call in calls] == ["mhcflurry", "netmhcpan"]
    assert calls[1]["peptides"] == ("MAGEAPEP1",)
    row = long.iloc[0]
    assert row["affinity_nm"] == 150.0
    assert row["netmhcpan_affinity_nm"] == 25.0
    assert row["netmhcpan_affinity_percentile"] == 0.03


def test_panel_alleles_are_ordered_by_weighted_frequency(monkeypatch):
    monkeypatch.setattr(
        "tsarina.alleles.get_panel",
        lambda panel: ["HLA-A*24:02", "HLA-A*02:01"],
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.regions.REGION_POPULATIONS",
        {"test-region": 1.0},
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.regions.region_allele_frequencies",
        lambda: pd.DataFrame(
            {
                "region": ["test-region", "test-region"],
                "allele": ["HLA-A*24:02", "HLA-A*02:01"],
                "frequency": [0.1, 0.3],
            }
        ),
        raising=True,
    )

    df = spanning_pmhc_set(
        ctas=["MAGEA4"],
        panel="test-panel",
        max_percentile=10.0,
        peptides_per_cell=1,
    )
    assert list(df.columns) == ["cta", "HLA-A*02:01", "HLA-A*24:02"]
    assert df.attrs["allele_order"] == ["HLA-A*02:01", "HLA-A*24:02"]


# ── Edge cases ─────────────────────────────────────────────────────────


def test_no_ctas_pass_filter_returns_canonical_empty_frame():
    """If filters eliminate every CTA, we still return a well-formed
    (empty) wide table with the requested allele columns."""
    df = spanning_pmhc_set(
        ctas=[],
        alleles=["HLA-A*02:01"],
        max_percentile=2.0,
    )
    assert list(df.columns) == ["cta", "HLA-A*02:01"]
    assert df.empty


# ── on_progress callback ───────────────────────────────────────────────


def test_on_progress_default_is_silent():
    """Library contract: when on_progress is not supplied, the function
    runs without calling any progress callback (or printing anywhere).

    Confirmed indirectly by the function completing without error; the
    explicit-callback test below asserts the callback IS called on the
    opt-in path."""
    df = spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01"],
        max_percentile=2.0,
    )
    assert not df.empty  # sanity — stub produces rows


def test_on_progress_callback_receives_three_stages():
    messages: list[str] = []
    spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01", "HLA-A*24:02"],
        max_percentile=2.0,
        on_progress=messages.append,
    )
    assert messages[0].startswith("tsarina v")
    assert any("Panel inputs" in msg for msg in messages)
    assert any("Enumerating CTA-exclusive peptides" in msg for msg in messages)
    assert any("Loading public MS evidence" in msg for msg in messages)
    assert any("Scoring" in msg and "mhcflurry" in msg for msg in messages)
    assert any("Scored" in msg and msg.endswith("s.") for msg in messages)
    assert any("Built MS evidence tiers" in msg for msg in messages)
    assert any("Panel selected" in msg for msg in messages)


def test_on_progress_threads_into_peptide_resolution(monkeypatch):
    messages: list[str] = []

    def _exclusive(**kw):
        kw["on_progress"]("inner-peptide-progress")
        assert kw["progress_bar"] is False
        assert kw["gene_names"] == ["MAGEA4"]
        return _stub_peptides()

    monkeypatch.setattr("tsarina.peptides.cta_exclusive_peptides", _exclusive, raising=True)

    spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01"],
        max_percentile=2.0,
        on_progress=messages.append,
    )
    assert "inner-peptide-progress" in messages


def test_on_progress_reports_alias_gene_expansion(monkeypatch):
    messages: list[str] = []

    def _exclusive(**kw):
        assert kw["gene_names"] == ["CTAG1A", "CTAG1B"]
        return _stub_peptides()

    monkeypatch.setattr("tsarina.peptides.cta_exclusive_peptides", _exclusive, raising=True)

    spanning_pmhc_set(
        ctas=["NY-ESO-1"],
        alleles=["HLA-A*02:01"],
        max_percentile=2.0,
        on_progress=messages.append,
    )

    assert any("1 CTA target labels to 2 Ensembl genes" in message for message in messages)


def test_on_progress_pre_scoring_message_has_accurate_counts():
    """The pre-scoring message's numbers should match the actual peptide
    and allele pool sizes: stub has 2 peptides for {MAGEA4, PRAME} each
    = 4 unique peptides, crossed with 3 alleles = 12 predictions."""
    messages: list[str] = []
    spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02"],
        max_percentile=10.0,
        on_progress=messages.append,
    )
    pre = next(msg for msg in messages if msg.startswith("Scoring "))
    assert "4 peptides" in pre
    assert "3 alleles" in pre
    assert "12 predictions" in pre


def test_on_progress_summary_reports_filled_cell_count():
    """Summary message's filled-cell count should equal the count the
    output table actually carries (rows with a non-NaN peptide)."""
    messages: list[str] = []
    spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01", "HLA-A*24:02"],
        max_percentile=2.0,
        output_format="long",
        on_progress=messages.append,
    )
    summary = messages[-1]
    assert "4/4 CTA x HLA cells" in summary


def test_progress_bar_scores_in_chunks(monkeypatch):
    calls: list[tuple[str, ...]] = []

    def _chunked_scores(peptides, alleles, **kwargs):
        calls.append(tuple(alleles))
        return _stub_scores(peptides, alleles, **kwargs)

    monkeypatch.setattr("tsarina.scoring.score_presentation", _chunked_scores, raising=True)

    long = spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01", "HLA-A*24:02"],
        max_percentile=10.0,
        peptides_per_cell=1,
        output_format="long",
        progress_bar=True,
        score_chunk_size=1,
        progress_file=io.StringIO(),
    )
    assert not long.empty
    assert calls == [("HLA-A*02:01",), ("HLA-A*24:02",)]


def test_progress_bar_default_mhcflurry_scores_one_full_allele_batch(monkeypatch):
    from tsarina.spanning import _score_presentations

    calls: list[tuple[str, ...]] = []

    def _chunked_scores(peptides, alleles, **kwargs):
        calls.append(tuple(alleles))
        return _stub_scores(peptides, alleles, **kwargs)

    monkeypatch.setattr("tsarina.scoring.score_presentation", _chunked_scores, raising=True)
    alleles = [f"HLA-A*02:{i:02d}" for i in range(1, 10)]

    scores = _score_presentations(
        peptides=["PEPTIDE"],
        alleles=alleles,
        predictor="mhcflurry",
        progress_bar=True,
        score_chunk_size=None,
        progress_file=io.StringIO(),
    )

    assert not scores.empty
    assert [len(call) for call in calls] == [9]


def test_progress_bar_explicit_mhcflurry_chunk_size_is_honored(monkeypatch):
    from tsarina.spanning import _score_presentations

    calls: list[tuple[str, ...]] = []

    def _chunked_scores(peptides, alleles, **kwargs):
        calls.append(tuple(alleles))
        return _stub_scores(peptides, alleles, **kwargs)

    monkeypatch.setattr("tsarina.scoring.score_presentation", _chunked_scores, raising=True)
    alleles = [f"HLA-A*02:{i:02d}" for i in range(1, 6)]

    scores = _score_presentations(
        peptides=["PEPTIDE"],
        alleles=alleles,
        predictor="mhcflurry",
        progress_bar=True,
        score_chunk_size=2,
        progress_file=io.StringIO(),
    )

    assert not scores.empty
    assert [len(call) for call in calls] == [2, 2, 1]


def test_cli_handler_wires_on_progress_to_stderr(monkeypatch, capsys):
    """The CLI handler must pass a stderr-printing callback into
    spanning_pmhc_set.  Verify callback messages land on stderr (not
    stdout) and that the DataFrame itself still serializes to stdout."""
    import argparse

    from tsarina import cli_spanning

    captured_kwargs: dict = {}

    def _fake_spanning(*args, **kwargs):
        captured_kwargs.update(kwargs)
        if kwargs.get("on_progress"):
            kwargs["on_progress"]("fake-progress-message")
        return pd.DataFrame({"cta": ["MAGEA4"], "HLA-A*02:01": ["PEPTIDE9"]})

    monkeypatch.setattr("tsarina.spanning.spanning_pmhc_set", _fake_spanning)

    args = argparse.Namespace(
        cta_count=25,
        cta_rank_by="ms_cancer_peptide_count",
        ctas=None,
        cancer_rna_threshold=2.0,
        cancer_type_prevalence_floor=0.05,
        min_restriction_confidence=["HIGH", "MODERATE"],
        restriction_levels=None,
        selection_allowlist=["PRAME", "CTAG1A/CTAG1B", "MAGEA4"],
        exclude_vital_tissue_expression=True,
        vital_tissue_max_ntpm=2.0,
        exclude_non_magea4_mage_family=True,
        alleles=None,
        panel="global51_abc_ssa",
        lengths=(8, 9, 10, 11),
        ensembl_release=112,
        require_cta_exclusive=True,
        predictor="mhcflurry",
        monoallelic_ms_max_percentile=2.0,
        sample_allele_ms_max_percentile=1.0,
        unrestricted_ms_max_percentile=0.5,
        include_predicted_only=False,
        predicted_only_max_percentile=0.1,
        max_percentile=None,
        peptides_per_cell=3,
        group_identical_cta_pmhcs=True,
        group_identical_cta_peptide_sets=True,
        annotate_netmhcpan_affinity=False,
        iedb_path=None,
        cedar_path=None,
        format="wide",
        output=None,
        summary=True,
        progress=True,
        progress_bars=False,
        score_chunk_size=None,
        show_empty_ctas=False,
    )
    cli_spanning.handle(args)

    assert "on_progress" in captured_kwargs
    assert callable(captured_kwargs["on_progress"])
    assert captured_kwargs["cancer_rna_threshold"] == 2.0
    assert captured_kwargs["cancer_type_prevalence_floor"] == 0.05
    assert captured_kwargs["selection_allowlist"] == ("PRAME", "CTAG1A/CTAG1B", "MAGEA4")
    assert captured_kwargs["exclude_vital_tissue_expression"] is True
    assert captured_kwargs["vital_tissue_max_ntpm"] == 2.0
    assert captured_kwargs["include_empty_ctas"] is False
    assert captured_kwargs["group_identical_cta_pmhcs"] is True
    assert captured_kwargs["group_identical_cta_peptide_sets"] is True
    assert captured_kwargs["annotate_netmhcpan_affinity"] is False

    captured = capsys.readouterr()
    assert "fake-progress-message" in captured.err
    assert "fake-progress-message" not in captured.out
    assert "MAGEA4" in captured.out  # DataFrame serialized to stdout


def test_cli_handler_default_table_report(monkeypatch, capsys):
    import argparse

    from tsarina import cli_spanning

    def _fake_spanning(*args, **kwargs):
        if kwargs.get("on_progress"):
            kwargs["on_progress"]("fake-progress-message")
        df = pd.DataFrame(
            {
                "cta": ["MAGEA4", "MAGEA4"],
                "allele": ["HLA-A*02:01", "HLA-A*02:01"],
                "peptide": ["MAGEAPEP2", "MAGEAPEP1"],
                "length": [9, 9],
                "evidence_tier": ["unrestricted_ms", "unrestricted_ms"],
                "ms_hit_count": [2, 1],
                "ms_source_count": [2, 1],
                "ms_alleles": ["", ""],
                "ms_pmids": ["222;333", "111"],
                "ms_samples": ["", ""],
                "presentation_percentile": [1.5, 0.1],
                "presentation_score": [0.8, 0.9],
                "affinity_nm": [200.0, 60.0],
                "affinity_percentile": [0.9, 0.2],
                "netmhcpan_affinity_nm": [pd.NA, pd.NA],
                "netmhcpan_affinity_percentile": [pd.NA, pd.NA],
                "peptide_rank_in_cell": [1, 2],
            }
        )
        df.attrs["cta_order"] = ["MAGEA4"]
        df.attrs["allele_order"] = ["HLA-A*02:01"]
        df.attrs["cta_rank_values"] = {"MAGEA4": 50}
        df.attrs["allele_frequencies"] = {"HLA-A*02:01": 0.2}
        df.attrs["panel_summary"] = {
            "hla_allele_count": 1,
            "cta_count": 1,
            "selected_peptide_count": 2,
            "selected_row_count": 2,
            "filled_cell_count": 1,
            "possible_cell_count": 1,
            "filled_cell_fraction": 1.0,
            "average_peptides_per_filled_cell": 2.0,
            "evidence_tier_counts": {"unrestricted_ms": 2},
            "coverage_note": "test coverage note",
            "cta_coverage": [
                {
                    "cta": "MAGEA4",
                    "covered_hla_count": 1,
                    "selected_peptide_count": 2,
                    "estimated_population_coverage": 0.36,
                }
            ],
            "hla_coverage": [
                {
                    "allele": "HLA-A*02:01",
                    "weighted_allele_frequency": 0.2,
                    "covered_cta_count": 1,
                    "covered_cta_fraction": 1.0,
                    "selected_peptide_count": 2,
                }
            ],
        }
        return df

    monkeypatch.setattr("tsarina.spanning.spanning_pmhc_set", _fake_spanning)

    args = argparse.Namespace(
        cta_count=25,
        cta_rank_by="ms_cancer_peptide_count",
        ctas=None,
        cancer_rna_threshold=2.0,
        cancer_type_prevalence_floor=0.05,
        min_restriction_confidence=["HIGH", "MODERATE"],
        restriction_levels=None,
        selection_allowlist=["PRAME", "CTAG1A/CTAG1B", "MAGEA4"],
        exclude_vital_tissue_expression=True,
        vital_tissue_max_ntpm=2.0,
        exclude_non_magea4_mage_family=True,
        alleles=None,
        panel="global51_abc_ssa",
        lengths=(8, 9, 10, 11),
        ensembl_release=112,
        require_cta_exclusive=True,
        predictor="mhcflurry",
        monoallelic_ms_max_percentile=2.0,
        sample_allele_ms_max_percentile=1.0,
        unrestricted_ms_max_percentile=0.5,
        include_predicted_only=False,
        predicted_only_max_percentile=0.1,
        max_percentile=None,
        peptides_per_cell=3,
        group_identical_cta_pmhcs=True,
        group_identical_cta_peptide_sets=True,
        annotate_netmhcpan_affinity=False,
        iedb_path=None,
        cedar_path=None,
        format="table",
        output=None,
        summary=True,
        progress=True,
        progress_bars=False,
        score_chunk_size=None,
        show_empty_ctas=False,
    )
    cli_spanning.handle(args)

    captured = capsys.readouterr()
    assert "fake-progress-message" in captured.err
    assert "CTA rank" in captured.out
    assert "Sources" in captured.out
    assert "MAGEAPEP2" in captured.out
    assert "Summary" in captured.out
    assert "Expected Population Coverage Per CTA" in captured.out
    assert "CTA Coverage Per HLA" in captured.out


def test_size_within_target_envelope_for_default_args(monkeypatch):
    """The headline use case: defaults (25 CTAs x 27 alleles, max 2%)
    should produce a wide table whose count of filled cells is in the
    100-1000 spanning-set range. Stub gives exactly 25*27=675 candidate
    cells; with the stub's 0.1 / 1.5 percentile baselines, every cell
    passes the 2.0 cutoff, so all 675 are filled — comfortably in range.

    Uses ``monkeypatch.setattr`` (not raw attribute assignment) so the
    expanded stubs are torn down after this test and don't leak into the
    rest of the test session.
    """
    # Expand the stub CTA pool so cta_count=25 has enough candidates.
    extra_genes = [f"GENE{i:03d}" for i in range(25)]
    extra_peps = pd.DataFrame(
        {
            "peptide": [f"PEPGENE{i:03d}" for i in range(25)],
            "length": [9] * 25,
            "gene_name": extra_genes,
            "gene_id": [f"E{i + 10}" for i in range(25)],
        }
    )
    all_genes = ["MAGEA4", "PRAME", "CTAG1B", *extra_genes]
    cta_csv = pd.DataFrame(
        {
            "Symbol": all_genes,
            "Ensembl_Gene_ID": [f"E{i}" for i in range(28)],
            "filtered": ["true"] * 28,
            "never_expressed": ["false"] * 28,
            "restriction_confidence": ["HIGH"] * 28,
            "restriction": ["TESTIS"] * 28,
            "ms_cancer_peptide_count": list(range(28, 0, -1)),
            "ms_cta_exclusive_cancer_peptide_count": list(range(28, 0, -1)),
        }
    )

    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides",
        lambda **kw: pd.concat([_stub_peptides(), extra_peps], ignore_index=True),
        raising=True,
    )
    monkeypatch.setattr("tsarina.loader.cta_dataframe", lambda: cta_csv, raising=True)
    monkeypatch.setattr("tsarina.gene_sets.CTA_gene_names", lambda: set(all_genes), raising=True)

    df = spanning_pmhc_set(
        cta_count=25,
        panel="iedb27_ab",
        max_percentile=2.0,
        output_format="long",
    )
    assert 100 < len(df) <= 1000


# ── require_cta_exclusive routing ──────────────────────────────────────


def test_require_cta_exclusive_true_calls_exclusive(monkeypatch):
    """Default True path must dispatch to cta_exclusive_peptides, not
    cta_peptides."""
    calls = {"exclusive": 0, "non_exclusive": 0}

    def _exclusive(**kw):
        calls["exclusive"] += 1
        assert kw["gene_names"] == ["MAGEA4"]
        return _stub_peptides()

    def _non_exclusive(**kw):
        calls["non_exclusive"] += 1
        assert kw["gene_names"] == ["MAGEA4"]
        return _stub_peptides()

    monkeypatch.setattr("tsarina.peptides.cta_exclusive_peptides", _exclusive, raising=True)
    monkeypatch.setattr("tsarina.peptides.cta_peptides", _non_exclusive, raising=True)

    spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01"],
        require_cta_exclusive=True,
        max_percentile=10.0,
    )
    assert calls == {"exclusive": 1, "non_exclusive": 0}


def test_require_cta_exclusive_false_calls_cta_peptides(monkeypatch):
    """False path must dispatch to cta_peptides (full CTA k-mers, no
    exclusivity gate), never the exclusive fn."""
    calls = {"exclusive": 0, "non_exclusive": 0}

    def _exclusive(**kw):
        calls["exclusive"] += 1
        assert kw["gene_names"] == ["MAGEA4"]
        return _stub_peptides()

    def _non_exclusive(**kw):
        calls["non_exclusive"] += 1
        assert kw["gene_names"] == ["MAGEA4"]
        return _stub_peptides()

    monkeypatch.setattr("tsarina.peptides.cta_exclusive_peptides", _exclusive, raising=True)
    monkeypatch.setattr("tsarina.peptides.cta_peptides", _non_exclusive, raising=True)

    spanning_pmhc_set(
        ctas=["MAGEA4"],
        alleles=["HLA-A*02:01"],
        require_cta_exclusive=False,
        max_percentile=10.0,
    )
    assert calls == {"exclusive": 0, "non_exclusive": 1}


# ── Multi-length best-per-cell ─────────────────────────────────────────


def _stub_peptides_multi_length() -> pd.DataFrame:
    """Mixed 9-mer and 10-mer peptides so the best-per-cell selection
    must reach across lengths."""
    return pd.DataFrame(
        {
            "peptide": [
                "MAGEAPEP9M",  # 10-mer; stub gives it percentile 0.3
                "MAGEAPEP1",  # 9-mer; stub gives it percentile 0.1 (wins)
                "PRAMEPEP0X",  # 10-mer; stub gives it percentile 0.05 (wins)
                "PRAMEPEP2",  # 9-mer; stub gives it percentile 1.5
            ],
            "length": [10, 9, 10, 9],
            "gene_name": ["MAGEA4", "MAGEA4", "PRAME", "PRAME"],
            "gene_id": ["E1", "E1", "E2", "E2"],
        }
    )


def _stub_scores_multi_length(peptides, alleles, **kwargs) -> pd.DataFrame:
    """Deterministic percentile per peptide — chosen so the 10-mer wins
    for PRAME and the 9-mer wins for MAGEA4, regardless of allele."""
    percentile_by_peptide = {
        "MAGEAPEP9M": 0.3,
        "MAGEAPEP1": 0.1,
        "PRAMEPEP0X": 0.05,
        "PRAMEPEP2": 1.5,
    }
    rows = []
    for pep in peptides:
        pct = percentile_by_peptide.get(pep, 5.0)
        for allele in alleles:
            rows.append(
                {
                    "peptide": pep,
                    "allele": allele,
                    "presentation_percentile": pct,
                    "presentation_score": 0.95 - pct / 10,
                    "affinity_nm": 50.0 + pct * 100,
                    "affinity_percentile": pct / 2,
                }
            )
    return pd.DataFrame(rows)


def test_multi_length_best_per_cell_winner_can_be_any_length(monkeypatch):
    """With lengths=(9, 10), each cell's winner should be the
    lowest-percentile peptide across BOTH lengths — not length-biased."""
    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides",
        lambda **kw: _stub_peptides_multi_length(),
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.scoring.score_presentation", _stub_scores_multi_length, raising=True
    )

    df = spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01"],
        lengths=(9, 10),
        max_percentile=2.0,
        peptides_per_cell=1,
    )
    indexed = df.set_index("cta")
    # MAGEA4: MAGEAPEP1 (9-mer, 0.1%) beats MAGEAPEP9M (10-mer, 0.3%)
    assert indexed.loc["MAGEA4", "HLA-A*02:01"] == "MAGEAPEP1"
    # PRAME: PRAMEPEP0X (10-mer, 0.05%) beats PRAMEPEP2 (9-mer, 1.5%)
    assert indexed.loc["PRAME", "HLA-A*02:01"] == "PRAMEPEP0X"


def test_multi_length_long_format_preserves_winner_length(monkeypatch):
    """In long format, the length column should reflect the winning
    peptide's length, not a default."""
    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides",
        lambda **kw: _stub_peptides_multi_length(),
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.scoring.score_presentation", _stub_scores_multi_length, raising=True
    )

    long = spanning_pmhc_set(
        ctas=["MAGEA4", "PRAME"],
        alleles=["HLA-A*02:01"],
        lengths=(9, 10),
        max_percentile=2.0,
        peptides_per_cell=1,
        output_format="long",
    )
    rows = {r["cta"]: r for _, r in long.iterrows()}
    assert int(rows["MAGEA4"]["length"]) == 9  # MAGEAPEP1
    assert int(rows["PRAME"]["length"]) == 10  # PRAMEPEP0X
