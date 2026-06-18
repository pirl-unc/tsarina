from __future__ import annotations

import io
import json
import zipfile
from contextlib import contextmanager
from pathlib import Path

import pytest
from hitlist import downloads as hl_downloads

from tsarina import reference_data


@pytest.fixture
def isolated_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("TSARINA_DATA_DIR", str(tmp_path))
    # cache_dir() = get_data_dir(subdir="tsarina", envkey=...) / "sources"
    return tmp_path / "tsarina" / "sources"


def _stub_downloader(monkeypatch, content: bytes, counter: dict | None = None):
    """Replace hitlist's ``download_to_file`` with one that writes *content*.

    Mirrors the real helper's cache semantics: a cached *dest* (without ``force``)
    is left untouched, otherwise the post-decompress bytes are written to *dest* --
    so the offline tests exercise the same reuse behaviour.
    """

    def _impl(url, dest, *, label="", verbose=True, force=False, decompress=False):
        dest = Path(dest)
        if dest.exists() and not force:
            return dest
        if counter is not None:
            counter["n"] += 1
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(content)
        return dest

    # The registry lives in hitlist and calls download_to_file in *its* module
    # namespace, so patch there (not on reference_data).
    monkeypatch.setattr(hl_downloads, "download_to_file", _impl)


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
    def _boom(url, dest, *, label="", verbose=True, force=False, decompress=False):
        raise OSError("network down")

    monkeypatch.setattr(hl_downloads, "download_to_file", _boom)
    with pytest.raises(reference_data.ReferenceDataError):
        reference_data.download("ncbi_gene_info")


def test_reference_data_delegates_to_hitlist_registry():
    # The versioned-download machinery lives in hitlist; reference_data only
    # supplies the dataset defs + cache namespace.
    assert isinstance(reference_data._REGISTRY, hl_downloads.VersionedDatasetRegistry)
    assert issubclass(reference_data.ReferenceDataError, hl_downloads.VersionedDatasetError)


def test_download_through_real_helper_decompresses_zip(isolated_cache, monkeypatch):
    # End-to-end through the *real* hitlist download_to_file (urlopen stubbed):
    # an HPA-style .zip is streamed and expanded into the .tsv dest, proving the
    # tsarina <- hitlist reuse, not just the stub contract.
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("rna_tissue_consensus.tsv", b"Gene\tTissue\tnTPM\nTSPY1\ttestis\t99\n")

    @contextmanager
    def fake_urlopen(url, timeout=None):
        yield io.BytesIO(buf.getvalue())

    monkeypatch.setattr(hl_downloads.urllib.request, "urlopen", fake_urlopen)

    path = reference_data.download("hpa_rna_consensus", version="v23", verbose=False)
    assert path.name == "rna_tissue_consensus.tsv"
    assert path.read_bytes() == b"Gene\tTissue\tnTPM\nTSPY1\ttestis\t99\n"
    # The compressed archive temp is cleaned up after decompression.
    assert not path.with_name(path.name + ".zip").exists()

    # Manifest records the decompressed payload.
    manifest = json.loads((isolated_cache / "manifest.json").read_text())
    assert manifest["hpa_rna_consensus"]["bytes"] == path.stat().st_size


def test_gz_payload_kept_verbatim(isolated_cache, monkeypatch):
    # gene_info is a .gz whose local name keeps the .gz suffix, so the downloader
    # leaves it compressed (no decompress); the stub mirrors that with raw bytes.
    _stub_downloader(monkeypatch, b"\x1f\x8bRAWGZ")
    path = reference_data.download("ncbi_gene_info")
    assert path.read_bytes() == b"\x1f\x8bRAWGZ"
    assert path.name == "Homo_sapiens.gene_info.gz"
