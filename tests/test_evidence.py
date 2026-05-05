import pandas as pd

from tsarina.evidence import _compute_ms_restriction


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
