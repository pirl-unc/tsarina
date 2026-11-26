import warnings

import pandas as pd
import pytest

from tsarina import CTA_gene_ids
from tsarina.cancer_expression import (
    cta_cancer_expression_features,
    cta_tcga_expression_features,
    hpa_cancer_ihc_prevalence,
    hpa_cancer_rna_prevalence,
)
from tsarina.vaccine_inputs import Protein, VaccineInputs, rank_proteoforms


def test_cta_cancer_expression_features_count_sample_and_cancer_type_prevalence():
    rna = pd.DataFrame(
        {
            "gene_id": ["E1", "E1", "E1", "E2"],
            "symbol": ["CTA1", "CTA1", "CTA1", "CTA2"],
            "cancer": ["Cancer A (TCGA)", "Cancer A (validation)", "Cancer B (TCGA)", "Cancer A"],
            "cancer_type": ["Cancer A", "Cancer A", "Cancer B", "Cancer A"],
            "cohort": ["TCGA", "validation", "TCGA", ""],
            "samples": [100, 50, 100, 100],
            "mean_ptpm": [1.0, 2.0, 3.0, 0.0],
            "max_ptpm": [10.0, 20.0, 30.0, 0.0],
            "expressed_samples_ptpm_ge_2": [10, 20, 1, 0],
            "prevalence_ptpm_ge_2": [0.1, 0.4, 0.01, 0.0],
        }
    )
    ihc = pd.DataFrame(
        {
            "gene_id": ["E1", "E1"],
            "symbol": ["CTA1", "CTA1"],
            "cancer": ["Cancer A", "Cancer B"],
            "total": [10, 10],
            "detected": [5, 0],
            "medium_high": [2, 0],
            "prevalence_detected": [0.5, 0.0],
            "prevalence_medium_high": [0.2, 0.0],
        }
    )

    features = cta_cancer_expression_features(
        rna_prevalence=rna,
        ihc_prevalence=ihc,
        rna_threshold=2.0,
        cancer_type_prevalence_floor=0.05,
    ).set_index("symbol")

    assert features.loc["CTA1", "cancer_rna_total_samples"] == 250
    assert features.loc["CTA1", "cancer_rna_expressed_samples"] == 31
    assert features.loc["CTA1", "cancer_rna_sample_prevalence"] == pytest.approx(31 / 250)
    assert features.loc["CTA1", "cancer_rna_prevalent_cancer_type_count"] == 1
    assert features.loc["CTA1", "cancer_rna_max_cancer_type_prevalence"] == pytest.approx(0.2)
    assert features.loc["CTA1", "cancer_ihc_medium_high_cancer_type_count"] == 1
    assert (
        features.loc["CTA1", "tumor_prevalence_panel_score"]
        > features.loc["CTA2", "tumor_prevalence_panel_score"]
    )


def test_cta_tcga_expression_features_aggregates_tcga_cohort_only():
    rna = pd.DataFrame(
        {
            "gene_id": ["E1", "E1", "E1", "E2"],
            "symbol": ["CTA1", "CTA1", "CTA1", "CTA2"],
            "cancer": [
                "Cancer A (TCGA)",
                "Cancer B (TCGA)",
                "Cancer A (validation)",
                "Cancer A (TCGA)",
            ],
            "cancer_type": ["Cancer A", "Cancer B", "Cancer A", "Cancer A"],
            "cohort": ["TCGA", "TCGA", "validation", "TCGA"],
            "samples": [100, 50, 200, 80],
            "mean_ptpm": [1.0, 1.0, 5.0, 0.5],
            "max_ptpm": [50.0, 20.0, 200.0, 4.0],
            "expressed_samples_ptpm_ge_1": [40, 5, 150, 8],
            "expressed_samples_ptpm_ge_5": [25, 1, 100, 0],
        }
    )

    features = cta_tcga_expression_features(rna_prevalence=rna).set_index("symbol")

    # validation cohort row for CTA1 must be excluded (TCGA-only).
    assert features.loc["CTA1", "tcga_sample_count"] == 150
    assert features.loc["CTA1", "tcga_cancer_type_count"] == 2
    assert features.loc["CTA1", "tcga_expressed_samples_tpm_ge_1"] == 45
    assert features.loc["CTA1", "tcga_expressed_samples_tpm_ge_5"] == 26
    assert features.loc["CTA1", "tcga_pan_prevalence_tpm_ge_1"] == pytest.approx(45 / 150)
    assert features.loc["CTA1", "tcga_pan_prevalence_tpm_ge_5"] == pytest.approx(26 / 150)
    assert features.loc["CTA1", "tcga_max_ptpm"] == pytest.approx(50.0)
    # Cancer A has 40/100 = 0.40 prevalence at TPM>=1; Cancer B has 5/50 = 0.10.
    assert features.loc["CTA1", "tcga_top_cancer_type_tpm_ge_1"] == "Cancer A"
    assert features.loc["CTA1", "tcga_top_cancer_type_prevalence_tpm_ge_1"] == pytest.approx(0.4)


def test_cta_tcga_expression_features_returns_empty_without_tcga_rows():
    rna = pd.DataFrame(
        {
            "gene_id": ["E1"],
            "symbol": ["CTA1"],
            "cancer_type": ["Cancer A"],
            "cohort": ["validation"],
            "samples": [100],
            "expressed_samples_ptpm_ge_1": [40],
            "expressed_samples_ptpm_ge_5": [10],
            "max_ptpm": [50.0],
        }
    )

    features = cta_tcga_expression_features(rna_prevalence=rna)
    assert features.empty


def test_cta_cancer_expression_features_requires_requested_threshold():
    rna = pd.DataFrame(
        {
            "gene_id": ["E1"],
            "symbol": ["CTA1"],
            "cancer_type": ["Cancer A"],
            "samples": [1],
            "expressed_samples_ptpm_ge_1": [1],
        }
    )

    with pytest.raises(ValueError, match="threshold 2"):
        cta_cancer_expression_features(
            rna_prevalence=rna,
            ihc_prevalence=pd.DataFrame(),
            rna_threshold=2.0,
        )


def test_live_oncoref_hpa_prevalence_covers_current_ctas():
    """Issue #258: HPA cancer RNA/IHC prevalence delegates to OncoRef v25.1."""
    cta_ids = CTA_gene_ids()
    rna = hpa_cancer_rna_prevalence()
    ihc = hpa_cancer_ihc_prevalence()

    assert "symbol" in rna.columns
    assert "symbol" in ihc.columns
    assert len(cta_ids & set(rna["gene_id"])) >= 350
    assert len(cta_ids & set(ihc["gene_id"])) >= 190

    tcga = cta_tcga_expression_features().set_index("symbol")
    for symbol in ("PLAC1", "DUX4", "NUTM1", "PRM1", "LALBA", "FOXR2", "CTAG2", "MAGEA4"):
        assert symbol in tcga.index
        assert tcga.loc[symbol, "tcga_max_ptpm"] > 0.0
    assert tcga.loc["PLAC1", "tcga_pan_prevalence_tpm_ge_1"] > 0.15

    # Filtering by gene_ids and cohort works cleanly.
    magea4_tcga = hpa_cancer_rna_prevalence(gene_ids=["ENSG00000147381.10"], cohort="TCGA")
    assert not magea4_tcga.empty
    assert set(magea4_tcga["gene_id"]) == {"ENSG00000147381"}
    assert set(magea4_tcga["cohort"]) == {"TCGA"}
    assert set(magea4_tcga["symbol"]) == {"MAGEA4"}


def test_rank_proteoforms_has_no_pandas4_warning_and_rejects_duplicate_burden():
    base_kwargs = {
        "cta_gene_ids": {"G1"},
        "gene_keys": {"G1": "P1"},
        "proteins": [Protein("G1", "P1", "SYM1", "MAAAA")],
        "prevalence": pd.DataFrame(
            {
                "proteoform_key": ["P1"],
                "cancer_code": ["LUAD"],
                "prevalence_p95": [0.25],
                "n_samples": [100],
            }
        ),
        "cohorts": {"lung": ["LUAD"]},
    }
    ok_inputs = VaccineInputs(
        burden=pd.DataFrame({"burden_category": ["lung"], "world_mortality_pct": [18.0]}),
        **base_kwargs,
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        rank_proteoforms(ok_inputs)
    pandas4 = [w for w in caught if "verify_integrity" in str(w.message)]
    assert not pandas4

    dup_inputs = VaccineInputs(
        burden=pd.DataFrame(
            {"burden_category": ["lung", "lung"], "world_mortality_pct": [18.0, 12.0]}
        ),
        **base_kwargs,
    )
    with pytest.raises(ValueError, match="Duplicate burden_category"):
        rank_proteoforms(dup_inputs)
