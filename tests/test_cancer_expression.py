import pandas as pd
import pytest

from tsarina.cancer_expression import cta_cancer_expression_features


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
