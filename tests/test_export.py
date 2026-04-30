import pandas as pd

from tsarina.export import build_ms_support_maps


def test_build_ms_support_maps_uses_public_ms_loader(monkeypatch):
    calls = {}
    hits = pd.DataFrame(
        {
            "peptide": ["MSPEPTIDE"],
            "mhc_restriction": ["HLA-A*02:01"],
        }
    )

    def _fake_load_public_ms_hits(peptides, **kwargs):
        calls["peptides"] = peptides
        calls["kwargs"] = kwargs
        return hits

    def _fake_aggregate_per_peptide(df):
        assert df is hits
        return pd.DataFrame({"peptide": ["MSPEPTIDE"], "ms_hit_count": [1]})

    def _fake_aggregate_per_pmhc(df):
        assert df is hits
        return pd.DataFrame(
            {
                "peptide": ["MSPEPTIDE"],
                "mhc_restriction": ["HLA-A*02:01"],
                "ms_hit_count": [1],
            }
        )

    monkeypatch.setattr(
        "tsarina.ms_evidence.load_public_ms_hits",
        _fake_load_public_ms_hits,
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.export.aggregate_per_peptide",
        _fake_aggregate_per_peptide,
        raising=True,
    )
    monkeypatch.setattr(
        "tsarina.export.aggregate_per_pmhc",
        _fake_aggregate_per_pmhc,
        raising=True,
    )

    peptide_map, pmhc_map, summary = build_ms_support_maps(
        peptides={"MSPEPTIDE"},
        mhc_class="I",
    )

    assert calls["peptides"] == {"MSPEPTIDE"}
    assert calls["kwargs"]["mhc_class"] == "I"
    assert peptide_map["MSPEPTIDE"]["ms_hit_count"] == 1
    assert pmhc_map[("MSPEPTIDE", "HLA-A*02:01")]["ms_hit_count"] == 1
    assert summary["peptide"].tolist() == ["MSPEPTIDE"]
