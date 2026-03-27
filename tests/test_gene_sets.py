from perseus import (
    CTA_evidence,
    CTA_excluded_gene_ids,
    CTA_excluded_gene_names,
    CTA_filtered_gene_ids,
    CTA_filtered_gene_names,
    CTA_gene_ids,
    CTA_gene_names,
    CTA_never_expressed_gene_names,
    CTA_unfiltered_gene_ids,
    CTA_unfiltered_gene_names,
)


def test_gene_names_nonempty():
    assert len(CTA_gene_names()) > 200


def test_gene_ids_nonempty():
    assert len(CTA_gene_ids()) > 200


def test_expressed_is_strict_subset_of_filtered():
    assert CTA_gene_names() < CTA_filtered_gene_names()
    assert CTA_gene_ids() < CTA_filtered_gene_ids()


def test_filtered_is_strict_subset_of_unfiltered():
    assert CTA_filtered_gene_names() < CTA_unfiltered_gene_names()
    assert CTA_filtered_gene_ids() < CTA_unfiltered_gene_ids()


def test_expressed_and_never_expressed_partition_filtered():
    expressed = CTA_gene_names()
    never_expr = CTA_never_expressed_gene_names()
    filtered = CTA_filtered_gene_names()
    assert expressed & never_expr == set()
    assert expressed | never_expr == filtered


def test_filtered_and_excluded_partition_unfiltered():
    filtered = CTA_filtered_gene_names()
    excluded = CTA_excluded_gene_names()
    unfiltered = CTA_unfiltered_gene_names()
    assert filtered & excluded == set()
    assert filtered | excluded == unfiltered


def test_excluded_ids_consistent():
    assert CTA_excluded_gene_ids() == CTA_unfiltered_gene_ids() - CTA_filtered_gene_ids()


def test_evidence_row_count_matches_unfiltered():
    df = CTA_evidence()
    assert len(df) == len(CTA_unfiltered_gene_names())


def test_evidence_has_expected_columns():
    df = CTA_evidence()
    expected = [
        "Symbol",
        "Ensembl_Gene_ID",
        "source_databases",
        "protein_reproductive",
        "protein_thymus",
        "protein_reliability",
        "rna_reproductive",
        "rna_thymus",
        "rna_deflated_reproductive_frac",
        "rna_80_pct_filter",
        "rna_90_pct_filter",
        "rna_95_pct_filter",
        "rna_99_pct_filter",
        "filtered",
        "never_expressed",
        "biotype",
        "Canonical_Transcript_ID",
        "rna_max_ntpm",
    ]
    for col in expected:
        assert col in df.columns, f"Missing column: {col}"


def test_magea4_is_expressed_cta():
    assert "MAGEA4" in CTA_gene_names()


def test_magea4_id_is_expressed_cta():
    assert "ENSG00000147381" in CTA_gene_ids()
