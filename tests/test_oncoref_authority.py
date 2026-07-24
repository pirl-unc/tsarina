from pathlib import Path

import pandas as pd
from oncoref import cta_tissues as oncoref_tissues

import tsarina
import tsarina.tissues as tissues
from tsarina import CTA_unfiltered_gene_ids


def _package_data_dir() -> Path:
    return Path(tsarina.__file__).parent / "data"


def test_package_contains_no_cta_definition_or_proteoform_fallback_table():
    names = {path.name for path in _package_data_dir().iterdir()}
    assert "cancer-testis-antigens.csv" not in names
    assert "proteoform-groups.csv" not in names


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
    canonical_ids = CTA_unfiltered_gene_ids()

    for filename in (
        "hpa-cancer-rna-prevalence.csv",
        "hpa-cancer-ihc-prevalence.csv",
    ):
        frame = pd.read_csv(_package_data_dir() / filename)
        assert not (prohibited & set(frame.columns))
        assert set(frame["gene_id"].astype(str).str.split(".").str[0]) <= canonical_ids


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
