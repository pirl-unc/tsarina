from __future__ import annotations

import json
from pathlib import Path

import pytest

from tsarina import reference_data


@pytest.fixture
def isolated_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("TSARINA_DATA_DIR", str(tmp_path))
    # cache_dir() = get_data_dir(subdir="tsarina", envkey=...) / "sources"
    return tmp_path / "tsarina" / "sources"


def _stub_downloader(monkeypatch, content: bytes, counter: dict | None = None):
    """Replace datacache's fetch/decompress with one that writes *content*.

    datacache normally downloads + decompresses to ``full_path``; the stub
    simulates the post-decompress result so the tests stay offline.
    """

    def _impl(full_path, download_url, timeout=None, use_wget_if_available=False):
        if counter is not None:
            counter["n"] += 1
        Path(full_path).write_bytes(content)

    monkeypatch.setattr(
        reference_data.datacache.download,
        "_download_and_decompress_if_necessary",
        _impl,
    )


def test_resolve_version_defaults_and_errors():
    # HPA datasets default to the pinned version; NCBI to the rolling sentinel.
    assert reference_data.resolve_version("hpa_rna_consensus") == reference_data.DEFAULT_HPA_VERSION
    assert reference_data.resolve_version("ncbi_gene_info") == reference_data.ROLLING
    with pytest.raises(reference_data.ReferenceDataError):
        reference_data.resolve_version("hpa_normal_tissue", "v99")
    with pytest.raises(reference_data.ReferenceDataError):
        reference_data.resolve_version("does_not_exist")


def test_download_writes_file_and_records_manifest(isolated_cache, monkeypatch):
    _stub_downloader(monkeypatch, b"Gene\tTissue\tnTPM\n")

    path = reference_data.download("hpa_rna_consensus", version="v23")
    assert path.exists()
    assert path.read_bytes() == b"Gene\tTissue\tnTPM\n"
    assert path == isolated_cache / "hpa_rna_consensus" / "v23" / "rna_tissue_consensus.tsv"

    manifest = json.loads((isolated_cache / "manifest.json").read_text())
    rec = manifest["hpa_rna_consensus"]
    assert rec["version"] == "v23"
    assert rec["bytes"] == path.stat().st_size
    assert len(rec["sha256"]) == 64
    assert rec["url"].endswith(".zip")


def test_download_reuses_cache_unless_forced(isolated_cache, monkeypatch):
    counter = {"n": 0}
    _stub_downloader(monkeypatch, b"x", counter)

    reference_data.download("hpa_normal_tissue")
    reference_data.ensure("hpa_normal_tissue")  # cached -> no new fetch
    assert counter["n"] == 1
    reference_data.download("hpa_normal_tissue", force=True)
    assert counter["n"] == 2


def test_status_reflects_cache_state(isolated_cache, monkeypatch):
    _stub_downloader(monkeypatch, b"y")

    before = {r["name"]: r for r in reference_data.status()}
    assert before["hpa_rna_consensus"]["cached"] is False

    reference_data.download("hpa_rna_consensus")
    after = {r["name"]: r for r in reference_data.status()}
    assert after["hpa_rna_consensus"]["cached"] is True
    assert after["hpa_rna_consensus"]["cached_version"] == reference_data.DEFAULT_HPA_VERSION
    assert "v23" in after["hpa_rna_consensus"]["available_versions"]


def test_download_failure_raises_reference_error(isolated_cache, monkeypatch):
    def _boom(full_path, download_url, timeout=None, use_wget_if_available=False):
        raise OSError("network down")

    monkeypatch.setattr(
        reference_data.datacache.download, "_download_and_decompress_if_necessary", _boom
    )
    with pytest.raises(reference_data.ReferenceDataError):
        reference_data.download("ncbi_gene_info")


def test_gz_payload_kept_verbatim(isolated_cache, monkeypatch):
    # gene_info is a .gz whose local name keeps the .gz suffix (datacache leaves
    # it compressed); the stub mirrors that by writing the raw bytes.
    _stub_downloader(monkeypatch, b"\x1f\x8bRAWGZ")
    path = reference_data.download("ncbi_gene_info")
    assert path.read_bytes() == b"\x1f\x8bRAWGZ"
    assert path.name == "Homo_sapiens.gene_info.gz"
