import io
import json
import zipfile

import pytest

from tsarina import reference_data


@pytest.fixture
def isolated_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("TSARINA_DATA_DIR", str(tmp_path))
    return tmp_path / "sources"


def _fake_zip(member_name: str, content: bytes) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr(member_name, content)
    return buf.getvalue()


class _FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def read(self):
        return self._payload

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def test_resolve_version_defaults_and_errors():
    # HPA datasets default to the pinned version; NCBI to the rolling sentinel.
    assert reference_data.resolve_version("hpa_rna_consensus") == reference_data.DEFAULT_HPA_VERSION
    assert reference_data.resolve_version("ncbi_gene_info") == reference_data.ROLLING
    with pytest.raises(reference_data.ReferenceDataError):
        reference_data.resolve_version("hpa_normal_tissue", "v99")
    with pytest.raises(reference_data.ReferenceDataError):
        reference_data.resolve_version("does_not_exist")


def test_download_extracts_zip_and_records_manifest(isolated_cache, monkeypatch):
    payload = _fake_zip("rna_tissue_consensus.tsv", b"Gene\tTissue\tnTPM\n")
    monkeypatch.setattr(
        reference_data.urllib.request, "urlopen", lambda *a, **k: _FakeResponse(payload)
    )

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
    calls = {"n": 0}

    def _count(*a, **k):
        calls["n"] += 1
        return _FakeResponse(_fake_zip("normal_tissue.tsv", b"x"))

    monkeypatch.setattr(reference_data.urllib.request, "urlopen", _count)

    reference_data.download("hpa_normal_tissue")
    reference_data.ensure("hpa_normal_tissue")  # cached -> no new fetch
    assert calls["n"] == 1
    reference_data.download("hpa_normal_tissue", force=True)
    assert calls["n"] == 2


def test_status_reflects_cache_state(isolated_cache, monkeypatch):
    monkeypatch.setattr(
        reference_data.urllib.request,
        "urlopen",
        lambda *a, **k: _FakeResponse(_fake_zip("rna_tissue_consensus.tsv", b"y")),
    )
    before = {r["name"]: r for r in reference_data.status()}
    assert before["hpa_rna_consensus"]["cached"] is False

    reference_data.download("hpa_rna_consensus")
    after = {r["name"]: r for r in reference_data.status()}
    assert after["hpa_rna_consensus"]["cached"] is True
    assert after["hpa_rna_consensus"]["cached_version"] == reference_data.DEFAULT_HPA_VERSION
    assert "v23" in after["hpa_rna_consensus"]["available_versions"]


def test_non_zip_payload_written_verbatim(isolated_cache, monkeypatch):
    monkeypatch.setattr(
        reference_data.urllib.request, "urlopen", lambda *a, **k: _FakeResponse(b"\x1f\x8bRAWGZ")
    )
    path = reference_data.download("ncbi_gene_info")
    assert path.read_bytes() == b"\x1f\x8bRAWGZ"
    assert path.name == "Homo_sapiens.gene_info.gz"
