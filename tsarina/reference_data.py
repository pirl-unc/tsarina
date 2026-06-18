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

Cache dir resolution is delegated to **datacache** (the same openvax layer
pyensembl uses); the download + ``.zip``/``.gz`` decompress go through
**hitlist**'s ``download_to_file`` -- a streaming, progress-reporting downloader
shared with ``tsarina data fetch``, so reference fetches show download progress
and cache status instead of silently buffering the whole file in memory.  The
base dir is ``$TSARINA_DATA_DIR`` if set, else a platform-appropriate
``appdirs`` dir (e.g.
``~/Library/Caches/tsarina`` on macOS).  Layout (one subdir per dataset+version,
plus a JSON manifest this module adds on top of datacache)::

    <cache>/sources/<name>/<version>/<filename>
    <cache>/sources/manifest.json
"""

from __future__ import annotations

from pathlib import Path

import datacache
from hitlist.downloads import VersionedDatasetError, VersionedDatasetRegistry

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
        "default_version": DEFAULT_HPA_VERSION,
        "urls": {
            "v23": "https://v23.proteinatlas.org/download/rna_tissue_consensus.tsv.zip",
            "latest": "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip",
        },
    },
    "hpa_normal_tissue": {
        "description": "HPA IHC protein expression (normal_tissue, per tissue/cell type)",
        "filename": "normal_tissue.tsv",
        "kind": "hpa",
        "default_version": DEFAULT_HPA_VERSION,
        "urls": {
            "v23": "https://v23.proteinatlas.org/download/normal_tissue.tsv.zip",
        },
    },
    "ncbi_gene_info": {
        "description": "NCBI human gene_info (official symbols + synonyms)",
        "filename": "Homo_sapiens.gene_info.gz",
        "kind": "ncbi",
        "default_version": ROLLING,
        "urls": {
            "latest": (
                "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/"
                "Mammalia/Homo_sapiens.gene_info.gz"
            ),
        },
    },
}


class ReferenceDataError(VersionedDatasetError):
    """Raised for unknown datasets/versions or download failures."""


# ── Cache location ───────────────────────────────────────────────────────────


#: Env var pointing the cache at a custom location (mirrors pyensembl's
#: ``PYENSEMBL_CACHE_DIR``).  When unset, datacache picks a platform-appropriate
#: dir via ``appdirs`` (e.g. ``~/Library/Caches/tsarina`` on macOS).
CACHE_DIR_ENV_KEY = "TSARINA_DATA_DIR"


def cache_dir() -> Path:
    """Return the reference-data cache directory (created on demand).

    Delegates dir resolution to ``datacache.get_data_dir`` -- the same
    download/cache layer pyensembl uses -- so tsarina, pyensembl, and other
    openvax tools share one convention instead of hand-rolling cache paths.
    """
    base = datacache.get_data_dir(subdir="tsarina", envkey=CACHE_DIR_ENV_KEY)
    sources = Path(base) / "sources"
    datacache.ensure_dir(str(sources))
    return sources


# ── Dataset registry (machinery lives in hitlist) ────────────────────────────
#
# The versioned download + provenance-manifest + cache machinery is hitlist's
# ``VersionedDatasetRegistry`` (shared with ``tsarina data fetch``'s downloader),
# so tsarina supplies only the dataset *definitions* + its cache namespace
# instead of re-implementing a parallel data layer.  The on-disk layout is
# unchanged (``<cache>/<name>/<version>/<filename>`` + ``<cache>/manifest.json``),
# so existing caches keep working.

_REGISTRY = VersionedDatasetRegistry(
    REFERENCE_DATASETS, cache_dir=cache_dir, error_cls=ReferenceDataError
)


def resolve_version(name: str, version: str | None = None) -> str:
    """Return the concrete version for *name*, applying its default."""
    return _REGISTRY.resolve_version(name, version)


def is_hpa_dataset(name: str) -> bool:
    """True if *name* is version-pinned to an HPA release (honors ``--hpa-version``).

    Non-HPA datasets (e.g. the rolling NCBI gene_info) ignore an HPA version, so
    callers should not forward ``--hpa-version`` to them. This stays tsarina-side
    because ``--hpa-version`` forwarding is a tsarina CLI concern.
    """
    spec = REFERENCE_DATASETS.get(name)
    return bool(spec and spec["kind"] == "hpa")


def local_path(name: str, version: str | None = None) -> Path:
    """Return the expected cache path for *name*/*version* (may not exist yet)."""
    return _REGISTRY.local_path(name, version)


def is_cached(name: str, version: str | None = None) -> bool:
    return _REGISTRY.is_cached(name, version)


def download(
    name: str, version: str | None = None, *, force: bool = False, verbose: bool = True
) -> Path:
    """Download *name*/*version* into the cache and record it in the manifest.

    A cached copy is reused unless ``force=True``. Delegates to hitlist's
    ``VersionedDatasetRegistry`` (streaming + progress + decompress + manifest).
    """
    return _REGISTRY.download(name, version, force=force, verbose=verbose)


def ensure(name: str, version: str | None = None) -> Path:
    """Return a local path to *name*/*version*, downloading if absent."""
    return _REGISTRY.ensure(name, version)


def status() -> list[dict]:
    """Return one status row per dataset for ``tsarina reference list``."""
    return _REGISTRY.status()
