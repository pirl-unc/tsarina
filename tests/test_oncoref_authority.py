from pathlib import Path

from oncoref import cta_tissues as oncoref_tissues

import tsarina
import tsarina.tissues as tissues
from tsarina import hpa_cancer_ihc_prevalence, hpa_cancer_rna_prevalence


def _package_data_dir() -> Path:
    return Path(tsarina.__file__).parent / "data"


def test_package_contains_no_cta_definition_proteoform_or_hpa_fallback_tables():
    names = {path.name for path in _package_data_dir().iterdir()}
    assert "cancer-testis-antigens.csv" not in names
    assert "proteoform-groups.csv" not in names
    assert "hpa-cancer-rna-prevalence.csv" not in names
    assert "hpa-cancer-ihc-prevalence.csv" not in names


def test_cancer_prevalence_tables_are_features_not_cta_definitions():
    prohibited = {
        "source_databases",
        "passes_filters",
        "filtered",
        "never_expressed",
        "specificity_status",
        "specificity_action",
        "restriction",
        "restriction_confidence",
    }

    for frame in (
        hpa_cancer_rna_prevalence(gene_ids=["ENSG00000147381"]),
        hpa_cancer_ihc_prevalence(gene_ids=["ENSG00000147381"]),
    ):
        assert not frame.empty
        assert not (prohibited & set(frame.columns))


def test_cta_tissue_definitions_are_direct_oncoref_exports():
    names = (
        "CORE_REPRODUCTIVE_TISSUES",
        "EXTENDED_REPRODUCTIVE_TISSUES",
        "PERMISSIVE_REPRODUCTIVE_TISSUES",
        "HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS",
        "HPA_EXPRESSION_FLOOR_NTPM",
        "PROTEIN_RELIABILITY_ORDER",
        "adaptive_rna_threshold",
    )
    for name in names:
        assert getattr(tissues, name) is getattr(oncoref_tissues, name)
