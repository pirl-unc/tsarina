# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Versioned download + cache for the upstream curation reference data.

The bundled CTA table (``cancer-testis-antigens.csv``) is derived from large
upstream reference files -- HPA RNA consensus, HPA IHC ``normal_tissue``, and
NCBI ``gene_info``.  This module downloads and caches those files **with an
explicit version stamp**, so the regeneration of the bundled table is
reproducible and can't silently mix HPA releases (the ``www`` "latest" mirror,
for example, serves the RNA consensus but not ``normal_tissue`` -- only the
pinned ``v23`` archive carries both, so the two must be fetched as a matched
pair).

Surfaced on the CLI as ``tsarina reference {list,fetch,path}`` -- the
HPA/NCBI counterpart of the hitlist-backed ``tsarina data`` (IEDB/CEDAR/viral).

Cache layout (one subdir per dataset+version, plus a JSON manifest)::

    <cache>/sources/<name>/<version>/<filename>
    <cache>/sources/manifest.json

The cache dir is ``$TSARINA_DATA_DIR`` if set, else ``<repo>/.cache`` when run
from a source checkout, else ``~/.cache/tsarina``.
"""

from __future__ import annotations

import hashlib
import io
import json
import os
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

#: Default HPA release.  ``v23`` is the most recent HPA version whose download
#: mirror serves *both* ``rna_tissue_consensus`` and ``normal_tissue`` as a
#: matched pair; newer ``www``/``vNN`` mirrors drop ``normal_tissue``, which
#: would force a release-mismatched (RNA-new / protein-old) table.
DEFAULT_HPA_VERSION = "v23"

#: Sentinel "version" for continuously-updated, unversioned upstreams (NCBI).
ROLLING = "latest"

#: Registry of reference datasets.  ``urls`` maps a version label to a download
#: URL; ``.zip`` URLs are transparently extracted to ``filename``.
REFERENCE_DATASETS: dict[str, dict] = {
    "hpa_rna_consensus": {
        "description": "HPA RNA consensus tissue nTPM (per-tissue expression)",
        "filename": "rna_tissue_consensus.tsv",
        "kind": "hpa",
        "urls": {
            "v23": "https://v23.proteinatlas.org/download/rna_tissue_consensus.tsv.zip",
            "latest": "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip",
        },
    },
    "hpa_normal_tissue": {
        "description": "HPA IHC protein expression (normal_tissue, per tissue/cell type)",
        "filename": "normal_tissue.tsv",
        "kind": "hpa",
        "urls": {
            "v23": "https://v23.proteinatlas.org/download/normal_tissue.tsv.zip",
        },
    },
    "ncbi_gene_info": {
        "description": "NCBI human gene_info (official symbols + synonyms)",
        "filename": "Homo_sapiens.gene_info.gz",
        "kind": "ncbi",
        "urls": {
            "latest": (
                "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/"
                "Mammalia/Homo_sapiens.gene_info.gz"
            ),
        },
    },
}


class ReferenceDataError(RuntimeError):
    """Raised for unknown datasets/versions or download failures."""


# ── Cache location ───────────────────────────────────────────────────────────


def cache_dir() -> Path:
    """Return the reference-data cache directory (created on demand)."""
    env = os.environ.get("TSARINA_DATA_DIR")
    if env:
        base = Path(env).expanduser()
    else:
        repo_cache = Path(__file__).resolve().parent.parent / ".cache"
        # Use the repo's .cache only from a source checkout (where scripts/ lives).
        if (repo_cache.parent / "scripts").is_dir():
            base = repo_cache
        else:
            base = Path.home() / ".cache" / "tsarina"
    sources = base / "sources"
    sources.mkdir(parents=True, exist_ok=True)
    return sources


def _manifest_path() -> Path:
    return cache_dir() / "manifest.json"


def _read_manifest() -> dict:
    path = _manifest_path()
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return {}


def _write_manifest(manifest: dict) -> None:
    _manifest_path().write_text(json.dumps(manifest, indent=2, sort_keys=True))


# ── Dataset / version resolution ─────────────────────────────────────────────


def _dataset(name: str) -> dict:
    try:
        return REFERENCE_DATASETS[name]
    except KeyError:
        known = ", ".join(sorted(REFERENCE_DATASETS))
        raise ReferenceDataError(f"unknown dataset {name!r}; known: {known}") from None


def resolve_version(name: str, version: str | None = None) -> str:
    """Return the concrete version for *name*, applying per-kind defaults."""
    spec = _dataset(name)
    if version is None:
        version = DEFAULT_HPA_VERSION if spec["kind"] == "hpa" else ROLLING
    if version not in spec["urls"]:
        avail = ", ".join(sorted(spec["urls"]))
        raise ReferenceDataError(f"{name!r} has no version {version!r}; available: {avail}")
    return version


def local_path(name: str, version: str | None = None) -> Path:
    """Return the expected cache path for *name*/*version* (may not exist yet)."""
    version = resolve_version(name, version)
    spec = _dataset(name)
    return cache_dir() / name / version / spec["filename"]


def is_cached(name: str, version: str | None = None) -> bool:
    return local_path(name, version).exists()


# ── Download ─────────────────────────────────────────────────────────────────


def _extract_to(raw: bytes, url: str, dest: Path) -> None:
    """Write *raw* (a .zip or plain payload) to *dest*, unzipping if needed."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if url.endswith(".zip"):
        with zipfile.ZipFile(io.BytesIO(raw)) as zf:
            member = next(n for n in zf.namelist() if n.endswith(dest.suffix))
            dest.write_bytes(zf.read(member))
    else:
        dest.write_bytes(raw)


def download(name: str, version: str | None = None, *, force: bool = False) -> Path:
    """Download *name*/*version* into the cache and record it in the manifest.

    Returns the local path.  A cached copy is reused unless ``force=True``.
    """
    version = resolve_version(name, version)
    spec = _dataset(name)
    dest = local_path(name, version)
    url = spec["urls"][version]

    if dest.exists() and not force:
        return dest

    try:
        with urllib.request.urlopen(url, timeout=600) as response:
            raw = response.read()
    except Exception as e:  # surface any network/HTTP failure uniformly
        raise ReferenceDataError(f"failed to download {name} ({url}): {e}") from e

    _extract_to(raw, url, dest)

    manifest = _read_manifest()
    manifest[name] = {
        "version": version,
        "url": url,
        "path": str(dest),
        "bytes": dest.stat().st_size,
        "sha256": hashlib.sha256(dest.read_bytes()).hexdigest(),
        "downloaded_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    _write_manifest(manifest)
    return dest


def ensure(name: str, version: str | None = None) -> Path:
    """Return a local path to *name*/*version*, downloading if absent."""
    path = local_path(name, version)
    return path if path.exists() else download(name, version)


# ── Status (for the CLI) ─────────────────────────────────────────────────────


def status() -> list[dict]:
    """Return one status row per dataset for ``tsarina reference list``."""
    manifest = _read_manifest()
    rows = []
    for name, spec in sorted(REFERENCE_DATASETS.items()):
        default_v = DEFAULT_HPA_VERSION if spec["kind"] == "hpa" else ROLLING
        path = cache_dir() / name / default_v / spec["filename"]
        record = manifest.get(name, {})
        rows.append(
            {
                "name": name,
                "description": spec["description"],
                "default_version": default_v,
                "available_versions": sorted(spec["urls"]),
                "cached": path.exists(),
                "cached_version": record.get("version") if record else None,
                "bytes": record.get("bytes") if path.exists() else None,
                "downloaded_at": record.get("downloaded_at") if path.exists() else None,
                "path": str(path),
            }
        )
    return rows
