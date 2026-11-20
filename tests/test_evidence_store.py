"""Tests for Hitlist 2.0 verified EvidenceStore MS adapter and CLI integration."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from hitlist.mhc_evidence import NATIVE_FIELDS, write_mhc_evidence_package
from hitlist.provenance import CONTRIBUTOR_COLUMNS, file_digest

from tsarina.evidence_store import (
    clear_evidence_store_cache,
    evidence_store_provenance,
    load_evidence_store_ms_hits,
    open_evidence_store,
    resolve_evidence_package,
)
from tsarina.indexing import load_ms_evidence
from tsarina.ms_evidence import cta_healthy_tissue_ms_hits, load_public_ms_hits


def _native_row(
    *,
    provenance_id: str,
    peptide: str,
    mhc_restriction: str = "HLA-A*02:01",
    mhc_class: str = "I",
    mhc_species: str = "Homo sapiens",
    assay_method: str = "mass spectrometry",
    assay_modality: str = "ms",
    ms_outcome_status: str = "positive_result",
    is_binding_assay: bool = False,
    is_ms_observation: bool = True,
    is_non_peptide_ligand: bool = False,
    src_cancer: bool = True,
    src_healthy_tissue: bool = False,
    source_tissue: str = "melanoma",
    pmid: str = "30000001",
    qualitative_measurement: str = "Positive",
    quantitative_measurement: str = "",
    quantitative_value: float | None = None,
    measurement_units: str = "",
    measurement_inequality: str = "",
) -> dict:
    row = {
        "peptide": peptide,
        "peptide_modifications": "",
        "has_ptm": False,
        "peptide_extended": peptide,
        "mhc_restriction": mhc_restriction,
        "mhc_species": mhc_species,
        "mhc_species_source": "reported",
        "mhc_species_context_disagrees": False,
        "mhc_class": mhc_class,
        "mhc_class_reported": mhc_class,
        "mhc_class_source": "reported",
        "mhc_class_corrected": False,
        "allele_resolution": "4-digit",
        "serotype": "HLA-A2",
        "serotypes": "HLA-A2",
        "serotype_source": "derived",
        "assay_iri": f"http://iedb.org/assay/{provenance_id}",
        "reference_iri": f"http://iedb.org/ref/{pmid}",
        "pmid": pmid,
        "submission_id": "",
        "reference_title": "Test study",
        "source_organism": "Homo sapiens",
        "species": mhc_species,
        "host": "Homo sapiens",
        "host_age": "",
        "host_mhc_types": mhc_restriction,
        "process_type": "Natural processing",
        "disease": "melanoma" if src_cancer else "healthy",
        "disease_stage": "",
        "source_tissue": source_tissue,
        "cell_name": "",
        "culture_condition": "",
        "antigen_processing_comments": "",
        "assay_comments": "",
        "qualitative_measurement": qualitative_measurement,
        "assay_method": assay_method,
        "response_measured": "ligand presentation"
        if is_ms_observation
        else "half maximal inhibitory concentration (IC50)",
        "measurement_units": measurement_units,
        "measurement_inequality": measurement_inequality,
        "quantitative_measurement": quantitative_measurement,
        "quantitative_value": quantitative_value,
        "assay_modality": assay_modality,
        "assay_modality_source": "method",
        "is_binding_assay": is_binding_assay,
        "is_ms_observation": is_ms_observation,
        "ms_outcome_status": ms_outcome_status,
        "is_non_peptide_ligand": is_non_peptide_ligand,
        "src_cancer": src_cancer,
        "src_adjacent_to_tumor": False,
        "src_activated_apc": False,
        "src_healthy_tissue": src_healthy_tissue,
        "src_healthy_thymus": False,
        "src_healthy_reproductive": False,
        "src_healthy_reproductive_female": False,
        "src_healthy_reproductive_male": False,
        "src_cell_line": False,
        "src_ebv_lcl": False,
        "src_ex_vivo": False,
        "cell_line_name": "",
        "cell_type": "",
        "is_monoallelic": True,
        "monoallelic_host": "",
        "restriction_evidence": "monoallelic" if is_ms_observation else "experimental",
        "mhc_allele_set": mhc_restriction,
        "mhc_allele_provenance": "exact",
        "mhc_allele_set_size": 1,
        "attributed_sample_label": "",
        "provenance_id": provenance_id,
        "source": "iedb",
        "supplementary_file": "",
        "is_potential_contaminant": False,
    }
    assert tuple(row) == NATIVE_FIELDS
    return row


def _write_native_table(path: Path, rows: list[dict]) -> None:
    if rows:
        df = pd.DataFrame(rows)
        df["quantitative_value"] = pd.to_numeric(df["quantitative_value"], errors="coerce")
        table = pa.Table.from_pandas(df, preserve_index=False)
    else:
        schema = pa.schema(
            [
                (
                    col,
                    pa.bool_()
                    if col.startswith(("is_", "has_", "src_"))
                    or col in ("mhc_species_context_disagrees", "mhc_class_corrected")
                    else pa.float64()
                    if col == "quantitative_value"
                    else pa.int64()
                    if col == "mhc_allele_set_size"
                    else pa.string(),
                )
                for col in NATIVE_FIELDS
            ]
        )
        table = pa.Table.from_pylist([], schema=schema)
    pq.write_table(table, path)


@pytest.fixture
def verified_mhc_package(tmp_path: Path) -> Path:
    """Build a real, hash-verified Hitlist 2.x MHC evidence package."""
    clear_evidence_store_cache()
    src = tmp_path / "source_corpus"
    src.mkdir()

    obs_rows = [
        _native_row(
            provenance_id="obs-1",
            peptide="SLLQHLIGL",
            src_cancer=True,
            src_healthy_tissue=False,
            source_tissue="melanoma",
            pmid="30000001",
        ),
        _native_row(
            provenance_id="obs-2",
            peptide="ALYVDSLFFL",
            src_cancer=False,
            src_healthy_tissue=True,
            source_tissue="heart muscle",
            pmid="30000002",
        ),
        _native_row(
            provenance_id="obs-3",
            peptide="SLLQHLIGL",
            mhc_class="II",
            mhc_restriction="HLA-DRB1*01:01",
            src_cancer=True,
            src_healthy_tissue=False,
            source_tissue="melanoma",
            pmid="30000003",
        ),
        _native_row(
            provenance_id="obs-4",
            peptide="GLC_LIPID",
            mhc_restriction="CD1d",
            is_non_peptide_ligand=True,
            src_cancer=False,
            src_healthy_tissue=False,
            source_tissue="blood",
            pmid="30000004",
        ),
        _native_row(
            provenance_id="obs-5",
            peptide="SLLQHLIGL",
            mhc_species="Mus musculus",
            mhc_restriction="H2-Kb",
            src_cancer=True,
            src_healthy_tissue=False,
            source_tissue="melanoma",
            pmid="30000005",
        ),
    ]
    binding_rows = [
        _native_row(
            provenance_id="bind-1",
            peptide="SLLQHLIGL",
            assay_method="purified MHC/competitive/fluorescence",
            assay_modality="binding",
            ms_outcome_status="",
            is_binding_assay=True,
            is_ms_observation=False,
            src_cancer=False,
            src_healthy_tissue=False,
            source_tissue="",
            pmid="30000010",
            quantitative_measurement="45",
            quantitative_value=45.0,
            measurement_units="nM",
            measurement_inequality="=",
        ),
    ]
    other_rows = [
        _native_row(
            provenance_id="other-1",
            peptide="SLLQHLIGL",
            assay_method="x-ray crystallography",
            assay_modality="structural",
            ms_outcome_status="",
            is_binding_assay=False,
            is_ms_observation=False,
            src_cancer=False,
            src_healthy_tissue=False,
            source_tissue="",
            pmid="30000020",
        ),
    ]

    _write_native_table(src / "observations.parquet", obs_rows)
    _write_native_table(src / "binding.parquet", binding_rows)
    _write_native_table(src / "other_assays.parquet", other_rows)

    all_ids = [r["provenance_id"] for r in obs_rows + binding_rows + other_rows]
    contrib_rows = [
        {
            "provenance_id": pid,
            "source_record_id": f"rec-{idx}",
            "source_dataset": "iedb",
            "source_row": idx + 1,
            "source_file": "iedb.csv",
            "source_key": f"key-{idx}",
            "original_fields": '["peptide"]',
            "source_row_values": '["SLLQHLIGL"]',
            "relationship_status": "exact",
        }
        for idx, pid in enumerate(all_ids)
    ]
    # Keep exact CONTRIBUTOR_COLUMNS ordering.
    contrib_df = pd.DataFrame(
        [{col: row.get(col, "x") for col in CONTRIBUTOR_COLUMNS} for row in contrib_rows]
    )
    pq.write_table(
        pa.Table.from_pandas(contrib_df, preserve_index=False),
        src / "observation_contributors.parquet",
    )

    mappings_df = pd.DataFrame(
        [
            {
                "peptide": "SLLQHLIGL",
                "gene_name": "PRAME",
                "gene_id": "ENSG00000185686",
                "protein_id": "ENSP00000381762",
            },
            {
                "peptide": "ALYVDSLFFL",
                "gene_name": "PRAME",
                "gene_id": "ENSG00000185686",
                "protein_id": "ENSP00000381762",
            },
        ]
    )
    pq.write_table(
        pa.Table.from_pandas(mappings_df, preserve_index=False),
        src / "peptide_mappings.parquet",
    )

    files_meta = {
        name: file_digest(src / name)
        for name in (
            "observations.parquet",
            "binding.parquet",
            "other_assays.parquet",
            "observation_contributors.parquet",
            "peptide_mappings.parquet",
        )
    }
    (src / "observations_meta.json").write_text(
        json.dumps(
            {
                "ci_corpus_manifest_version": 2,
                "artifact_version": 10,
                "files": files_meta,
                "provenance": {
                    "sources": {
                        "iedb": {"sha256": "0" * 64, "size": 123},
                    }
                },
            }
        )
    )

    pkg = tmp_path / "verified_pkg"
    write_mhc_evidence_package(src, pkg)
    return pkg


def test_evidence_store_loads_filtered_positive_ms_with_provenance(verified_mhc_package: Path):
    hits = load_evidence_store_ms_hits(
        verified_mhc_package,
        peptides={"SLLQHLIGL", "ALYVDSLFFL", "GLC_LIPID"},
        mhc_class="I",
        mhc_species="Homo sapiens",
    )
    # Class II (obs-3), non-peptide lipid (obs-4), mouse (obs-5) and binding (bind-1) are excluded.
    assert sorted(hits["peptide"]) == ["ALYVDSLFFL", "SLLQHLIGL"]
    assert set(hits["evidence_type"]) == {"mhc_positive_ms"}
    assert set(hits["evidence_kind"]) == {"ms"}
    assert hits.loc[hits["peptide"] == "SLLQHLIGL", "src_cancer"].all()
    assert hits.loc[hits["peptide"] == "ALYVDSLFFL", "src_healthy_tissue"].all()
    assert hits.attrs["evidence_backend"] == "hitlist_evidence_store"
    assert hits.attrs["evidence_package_manifest_sha256"]
    assert hits.attrs["evidence_package_schema_version"] == 2


def test_evidence_store_include_binding_and_length_filters(verified_mhc_package: Path):
    hits = load_evidence_store_ms_hits(
        verified_mhc_package,
        peptides={"SLLQHLIGL", "ALYVDSLFFL"},
        mhc_class="I",
        mhc_species="Homo sapiens",
        include_binding=True,
        drop_binding_assays=False,
        length_min=9,
        length_max=9,
    )
    assert set(hits["peptide"]) == {"SLLQHLIGL"}
    assert set(hits["evidence_kind"]) == {"ms", "binding"}


def test_load_ms_evidence_and_public_hits_route_through_evidence_package(
    verified_mhc_package: Path, monkeypatch
):
    direct = load_ms_evidence(
        peptides={"SLLQHLIGL"},
        evidence_package=verified_mhc_package,
    )
    assert len(direct) == 1
    assert direct["pmid"].iloc[0] == "30000001"
    assert direct.attrs["evidence_backend"] == "hitlist_evidence_store"

    # Environment variable routing avoids touching observations.parquet.
    monkeypatch.setenv("TSARINA_EVIDENCE_PACKAGE", str(verified_mhc_package))
    assert resolve_evidence_package() == verified_mhc_package
    via_env = load_public_ms_hits({"ALYVDSLFFL"})
    assert len(via_env) == 1
    assert via_env["source_tissue"].iloc[0] == "heart muscle"
    assert via_env.attrs["evidence_backend"] == "hitlist_evidence_store"


def test_cta_healthy_tissue_ms_hits_uses_verified_package_gene_mappings(verified_mhc_package: Path):
    healthy = cta_healthy_tissue_ms_hits("PRAME", evidence_package=verified_mhc_package)
    assert len(healthy) == 1
    assert healthy["peptide"].iloc[0] == "ALYVDSLFFL"
    assert healthy["tissue"].iloc[0] == "heart muscle"
    assert bool(healthy["vital_organ"].iloc[0]) is True


def test_cli_hits_supports_evidence_package_flag(verified_mhc_package: Path, monkeypatch, capsys):
    from tsarina.cli import main

    monkeypatch.setattr(
        "sys.argv",
        [
            "tsarina",
            "hits",
            "--gene",
            "PRAME",
            "--evidence-package",
            str(verified_mhc_package),
            "--format",
            "raw",
        ],
    )
    main()
    out = capsys.readouterr().out
    assert "SLLQHLIGL" in out
    assert "ALYVDSLFFL" in out


def test_tampered_evidence_package_fails_verification(verified_mhc_package: Path):
    clear_evidence_store_cache()
    target = verified_mhc_package / "typed" / "mhc_observations.parquet"
    target.write_bytes(target.read_bytes() + b"corrupt")
    with pytest.raises(ValueError, match=r"digest|mismatch|hash"):
        open_evidence_store(verified_mhc_package)


def test_evidence_store_provenance_reports_counts_and_hash(verified_mhc_package: Path):
    prov = evidence_store_provenance(verified_mhc_package)
    assert prov["evidence_backend"] == "hitlist_evidence_store"
    assert prov["mhc_positive_ms_rows"] == 5
    assert prov["mhc_binding_rows"] == 1
    assert len(prov["evidence_package_manifest_sha256"]) == 64


def test_full_evidence_store_open_path_is_invoked_for_multimodality_packages(
    tmp_path: Path, monkeypatch
):
    clear_evidence_store_cache()
    pkg = tmp_path / "full_pkg"
    pkg.mkdir()
    (pkg / "manifest.json").write_text(
        json.dumps({"schema_version": 10, "hitlist_version": "2.1.22"})
    )
    sentinel = SimpleNamespace(directory=pkg, manifest={"schema_version": 10}, bindings={})
    monkeypatch.setattr("hitlist.EvidenceStore.open", lambda path: sentinel)
    assert open_evidence_store(pkg) is sentinel
