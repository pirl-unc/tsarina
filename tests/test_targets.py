import pandas as pd

from tsarina.targets import target_peptides, target_summary


def _cta_targets() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "peptide": ["MSPEPTIDE", "NOEVIDENC"],
            "length": [9, 9],
            "gene_name": ["MAGEA4", "MAGEA4"],
            "gene_id": ["ENSG1", "ENSG1"],
        }
    )


def test_target_peptides_require_ms_evidence_uses_default_loader(monkeypatch):
    calls = {}

    monkeypatch.setattr("tsarina.peptides.cta_exclusive_peptides", lambda **kw: _cta_targets())

    def _fake_load_public_ms_hits(peptides, **kwargs):
        calls["peptides"] = peptides
        calls["kwargs"] = kwargs
        return pd.DataFrame(
            {
                "peptide": ["MSPEPTIDE"],
                "mhc_restriction": ["HLA-A*02:01"],
                "src_cancer": [True],
                "src_healthy_tissue": [False],
            }
        )

    monkeypatch.setattr(
        "tsarina.ms_evidence.load_public_ms_hits",
        _fake_load_public_ms_hits,
        raising=True,
    )

    out = target_peptides(
        cta=True,
        viruses=False,
        mutations=False,
        require_ms_evidence=True,
    )

    assert calls["peptides"] == {"MSPEPTIDE", "NOEVIDENC"}
    assert calls["kwargs"]["iedb_path"] is None
    assert calls["kwargs"]["cedar_path"] is None
    assert calls["kwargs"]["mhc_species"] == "Homo sapiens"
    assert out["peptide"].tolist() == ["MSPEPTIDE"]
    assert out["has_ms_evidence"].tolist() == [True]
    assert out["ms_hit_count"].tolist() == [1]


def test_target_peptides_cancer_specific_uses_default_loader(monkeypatch):
    monkeypatch.setattr("tsarina.peptides.cta_exclusive_peptides", lambda **kw: _cta_targets())

    def _fake_load_public_ms_hits(peptides, **kwargs):
        return pd.DataFrame(
            {
                "peptide": ["MSPEPTIDE", "NOEVIDENC"],
                "mhc_restriction": ["HLA-A*02:01", "HLA-A*02:01"],
                "src_cancer": [True, True],
                "src_healthy_tissue": [False, True],
            }
        )

    monkeypatch.setattr(
        "tsarina.ms_evidence.load_public_ms_hits",
        _fake_load_public_ms_hits,
        raising=True,
    )

    out = target_peptides(
        cta=True,
        viruses=False,
        mutations=False,
        cancer_specific=True,
    )

    assert out["peptide"].tolist() == ["MSPEPTIDE"]
    assert out["ms_cancer"].tolist() == [True]
    assert out["ms_healthy_tissue"].tolist() == [False]


def test_target_summary_on_empty():
    df = pd.DataFrame(columns=["peptide", "length", "category", "source", "source_detail"])
    summary = target_summary(df)
    assert isinstance(summary, pd.DataFrame)
    assert len(summary) == 0


def test_target_summary_categories():
    df = pd.DataFrame(
        {
            "peptide": ["AAAA", "BBBB", "CCCC", "DDDD"],
            "length": [4, 4, 4, 4],
            "category": ["cta", "cta", "viral", "mutant"],
            "source": ["MAGEA4", "MAGEA4", "HPV-16", "KRAS G12D"],
            "source_detail": ["ENSG1", "ENSG1", "P1", "G12D"],
            "has_ms_evidence": [True, False, True, False],
            "ms_cancer": [True, False, False, False],
            "ms_healthy_tissue": [False, False, False, False],
            "ms_cell_lines": ["HeLa", "", "A549", ""],
        }
    )
    summary = target_summary(df)
    assert len(summary) == 3
    assert set(summary["category"]) == {"cta", "viral", "mutant"}
