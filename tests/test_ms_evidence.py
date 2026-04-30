from pathlib import Path

import pandas as pd

from tsarina.ms_evidence import (
    aggregate_ms_hits_by_peptide,
    aggregate_ms_hits_for_iedb_columns,
    load_public_ms_hits,
)


def test_load_public_ms_hits_raw_path_drops_binding_assays(monkeypatch):
    calls = {}

    def _fake_resolve(iedb_path, cedar_path, require_iedb=True):
        calls["resolve"] = (iedb_path, cedar_path, require_iedb)
        return Path("iedb.csv"), None

    def _fake_scan_public_ms(**kwargs):
        calls["scan"] = kwargs
        return pd.DataFrame(
            {
                "peptide": ["MSPEPTIDE", "BINDPEPTD"],
                "mhc_restriction": ["HLA-A*02:01", "HLA-A*02:01"],
                "is_binding_assay": [False, True],
            }
        )

    monkeypatch.setattr("tsarina.datasources.resolve_dataset_paths", _fake_resolve, raising=True)
    monkeypatch.setattr("tsarina.iedb.scan_public_ms", _fake_scan_public_ms, raising=True)

    out = load_public_ms_hits(
        peptides={"MSPEPTIDE", "BINDPEPTD"},
        iedb_path="iedb.csv",
        mhc_class="I",
    )

    assert calls["resolve"] == ("iedb.csv", None, True)
    assert calls["scan"]["mhc_species"] == "Homo sapiens"
    assert out["peptide"].tolist() == ["MSPEPTIDE"]


def test_aggregate_ms_hits_by_peptide_carries_source_flags_and_cell_lines():
    hits = pd.DataFrame(
        {
            "peptide": ["MSPEPTIDE", "MSPEPTIDE"],
            "mhc_restriction": ["HLA-A*02:01", "HLA-B*07:02"],
            "src_cancer": [True, False],
            "src_healthy_tissue": [False, True],
            "cell_line_name": ["A375", ""],
        }
    )

    out = aggregate_ms_hits_by_peptide(hits)

    row = out.iloc[0]
    assert row["ms_hit_count"] == 2
    assert row["ms_alleles"] == "HLA-A*02:01;HLA-B*07:02"
    assert row["ms_allele_count"] == 2
    assert bool(row["ms_cancer"]) is True
    assert bool(row["ms_healthy_tissue"]) is True
    assert row["ms_cell_lines"] == "A375"


def test_aggregate_ms_hits_by_peptide_ignores_blank_and_missing_strings():
    hits = pd.DataFrame(
        {
            "peptide": ["MSPEPTIDE", "MSPEPTIDE", "MSPEPTIDE"],
            "mhc_restriction": ["HLA-A*02:01", None, ""],
            "cell_line_name": ["A375", None, ""],
        }
    )

    out = aggregate_ms_hits_by_peptide(hits)

    row = out.iloc[0]
    assert row["ms_hit_count"] == 3
    assert row["ms_alleles"] == "HLA-A*02:01"
    assert row["ms_allele_count"] == 1
    assert row["ms_cell_lines"] == "A375"


def test_aggregate_ms_hits_for_iedb_columns_keeps_legacy_shape():
    hits = pd.DataFrame(
        {
            "peptide": ["MSPEPTIDE"],
            "mhc_restriction": ["HLA-A*02:01"],
        }
    )

    out = aggregate_ms_hits_for_iedb_columns(hits)

    assert out.to_dict(orient="records") == [
        {
            "peptide": "MSPEPTIDE",
            "iedb_hit_count": 1,
            "iedb_alleles": "HLA-A*02:01",
        }
    ]
