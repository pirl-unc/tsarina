"""Downloaded saved reports must retain their original bytes after MkDocs."""

import json
from hashlib import sha256

import pytest

from tsarina.docs_reports import on_post_build


def test_preserve_report_markdown_and_manifest(tmp_path):
    docs, site = tmp_path / "docs", tmp_path / "site"
    source = docs / "example"
    source.mkdir(parents=True)
    original = b"# Original report\n\nDLA-88*001:01\n"
    (source / "report.md").write_bytes(original)
    manifest = {"artifact_sha256": {"report.md": sha256(original).hexdigest()}}
    (source / "manifest.json").write_text(json.dumps(manifest))
    on_post_build({"docs_dir": str(docs), "site_dir": str(site)})
    assert (site / "example/report.md").read_bytes() == original
    assert (site / "example/manifest.json").read_bytes() == (source / "manifest.json").read_bytes()
    (source / "report.md").write_text("changed after hashing")
    with pytest.raises(ValueError, match="failed verification"):
        on_post_build({"docs_dir": str(docs), "site_dir": str(site)})


def test_report_artifact_cannot_escape_directory(tmp_path):
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "manifest.json").write_text(json.dumps({"artifact_sha256": {"../other": "bad"}}))
    with pytest.raises(ValueError, match="escapes"):
        on_post_build({"docs_dir": str(docs), "site_dir": str(tmp_path / "site")})
