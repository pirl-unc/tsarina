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

"""Reproducible source-data fetchers for the CTA curation scripts.

Thin wrappers over :mod:`tsarina.reference_data` (the versioned download/cache
layer also exposed as ``tsarina data sources``), so the curation scripts and the
CLI share one version-stamped cache.  Each helper accepts an optional explicit
local path/URL override; otherwise it returns the cached pinned-version copy,
downloading on demand.

Inspect or pre-fetch from the CLI::

    tsarina data sources list
    tsarina data sources fetch
"""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
import sys  # noqa: E402

sys.path.insert(0, str(REPO_ROOT))

from tsarina import reference_data  # noqa: E402


def _resolve(name: str, source: str | None, version: str | None = None) -> Path:
    """Return a local path for *name*, honoring an explicit local/URL override."""
    if source is None:
        return reference_data.ensure(name, version=version)
    if not source.startswith(("http://", "https://")):
        # Explicit local file (e.g. a pre-downloaded .tsv) — use as-is.
        return Path(source)
    return _fetch_url(name, source)


def _fetch_url(name: str, url: str) -> Path:
    """Fetch an explicit URL override into a side cache slot (no manifest)."""
    import urllib.request

    spec = reference_data.REFERENCE_DATASETS[name]
    dest = reference_data.cache_dir() / name / "override" / spec["filename"]
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        raw = urllib.request.urlopen(url, timeout=600).read()
        reference_data._extract_to(raw, url, dest)
    return dest


def rna_tissue_consensus_path(source: str | None = None) -> Path:
    """Local path to HPA ``rna_tissue_consensus.tsv`` (pinned HPA version)."""
    return _resolve("hpa_rna_consensus", source)


def normal_tissue_path(source: str | None = None) -> Path:
    """Local path to HPA ``normal_tissue.tsv`` IHC table (pinned HPA version)."""
    return _resolve("hpa_normal_tissue", source)


def ncbi_gene_info_path(source: str | None = None) -> Path:
    """Local path to NCBI human ``gene_info.gz``."""
    return _resolve("ncbi_gene_info", source)
