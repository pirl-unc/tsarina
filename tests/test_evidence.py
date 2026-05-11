import pandas as pd

from tsarina.evidence import _attach_tcga_prevalence, _compute_ms_restriction


def test_compute_ms_restriction_uses_public_ms_loader(monkeypatch):
    calls = {}
    input_df = pd.DataFrame({"Symbol": ["MAGEA4", "PRAME"]})

    monkeypatch.setattr(
        "tsarina.peptides.cta_peptides",
        lambda **kwargs: pd.DataFrame(
            {
                "peptide": ["MSPEPTIDE", "OTHERPEP"],
                "gene_name": ["MAGEA4", "PRAME"],
            }
        ),
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides",
        lambda **kwargs: pd.DataFrame(
            {
                "peptide": ["MSPEPTIDE"],
                "gene_name": ["MAGEA4"],
            }
        ),
        raising=True,
    )

    def _fake_load_public_ms_hits(peptides, **kwargs):
        calls["peptides"] = peptides
        calls["kwargs"] = kwargs
        return pd.DataFrame(
            {
                "peptide": ["MSPEPTIDE"],
                "mhc_restriction": ["HLA-A*02:01"],
                "src_cancer": [True],
            }
        )

    def _fake_aggregate_gene_ms_safety(hits, gene_map, exclusive_peptide_gene_map=None):
        assert hits["peptide"].tolist() == ["MSPEPTIDE"]
        assert set(gene_map["gene_name"]) == {"MAGEA4", "PRAME"}
        assert exclusive_peptide_gene_map is not None
        assert exclusive_peptide_gene_map.to_dict("records") == [
            {"peptide": "MSPEPTIDE", "gene_name": "MAGEA4"}
        ]
        return pd.DataFrame(
            {
                "gene_name": ["MAGEA4"],
                "ms_restriction": ["CANCER_ONLY"],
                "ms_cancer_count": [1],
            }
        )

    monkeypatch.setattr(
        "tsarina.ms_evidence.load_public_ms_hits",
        _fake_load_public_ms_hits,
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.tiers.aggregate_gene_ms_safety",
        _fake_aggregate_gene_ms_safety,
        raising=True,
    )

    out = _compute_ms_restriction(input_df, iedb_path=None, cedar_path=None)

    assert calls["peptides"] == {"MSPEPTIDE", "OTHERPEP"}
    assert calls["kwargs"]["mhc_class"] == "I"
    assert out.loc[out["Symbol"] == "MAGEA4", "ms_restriction"].iloc[0] == "CANCER_ONLY"
    assert out.loc[out["Symbol"] == "PRAME", "ms_restriction"].iloc[0] == "NO_MS_DATA"


def test_attach_tcga_prevalence_joins_on_stripped_ensembl_id(monkeypatch):
    input_df = pd.DataFrame(
        {
            "Symbol": ["MAGEA4", "PRAME", "UNKNOWN"],
            # Versioned and unversioned Ensembl IDs must both match.
            "Ensembl_Gene_ID": ["ENSG00000147381.10", "ENSG00000185686", "ENSG99999999999"],
        }
    )
    fake_features = pd.DataFrame(
        {
            "gene_id": ["ENSG00000147381", "ENSG00000185686"],
            "symbol": ["MAGEA4", "PRAME"],
            "tcga_sample_count": [6918, 6918],
            "tcga_cancer_type_count": [21, 21],
            "tcga_max_ptpm": [1602.47, 538.025],
            "tcga_expressed_samples_tpm_ge_1": [1178, 3036],
            "tcga_expressed_samples_tpm_ge_5": [912, 2190],
            "tcga_pan_prevalence_tpm_ge_1": [0.1703, 0.4389],
            "tcga_pan_prevalence_tpm_ge_5": [0.1318, 0.3166],
            "tcga_top_cancer_type_tpm_ge_1": [
                "Lung Squamous Cell Carcinoma",
                "Skin Cuteneous Melanoma",
            ],
            "tcga_top_cancer_type_prevalence_tpm_ge_1": [0.6462, 0.9798],
        }
    )

    monkeypatch.setattr(
        "tsarina.cancer_expression.cta_tcga_expression_features",
        lambda: fake_features,
        raising=True,
    )

    out = _attach_tcga_prevalence(input_df)

    by_symbol = out.set_index("Symbol")
    assert (
        by_symbol.loc["MAGEA4", "tcga_top_cancer_type_tpm_ge_1"] == "Lung Squamous Cell Carcinoma"
    )
    assert by_symbol.loc["PRAME", "tcga_pan_prevalence_tpm_ge_5"] == 0.3166
    # Genes without TCGA features get zeros / empty strings, not NaN.
    assert by_symbol.loc["UNKNOWN", "tcga_sample_count"] == 0
    assert by_symbol.loc["UNKNOWN", "tcga_pan_prevalence_tpm_ge_1"] == 0.0
    assert by_symbol.loc["UNKNOWN", "tcga_top_cancer_type_tpm_ge_1"] == ""
