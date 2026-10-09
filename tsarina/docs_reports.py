"""MkDocs hook preserving downloadable, hash-verified saved report artifacts."""

import json
import shutil
from hashlib import sha256
from pathlib import Path


def on_post_build(config, **kwargs):
    """Keep original Markdown downloads alongside MkDocs' rendered pages."""
    docs = Path(config["docs_dir"]).resolve()
    site = Path(config["site_dir"]).resolve()
    for manifest in docs.rglob("manifest.json"):
        hashes = json.loads(manifest.read_text()).get("artifact_sha256")
        if hashes is None:
            continue
        source_dir = manifest.parent
        destination = site / source_dir.relative_to(docs)
        files = []
        for name, digest in hashes.items():
            source = (source_dir / name).resolve()
            target = (destination / name).resolve()
            if not source.is_relative_to(source_dir) or not target.is_relative_to(destination):
                raise ValueError(f"Saved report artifact escapes its directory: {name}")
            if not source.is_file() or sha256(source.read_bytes()).hexdigest() != digest:
                raise ValueError(f"Saved report artifact failed verification: {name}")
            files.append((source, target))
        files.append((manifest, destination / "manifest.json"))
        for source, target in files:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
