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

Each curation script (``add_xage2.py``, ``add_cta_gene.py``,
``backfill_aliases.py``) depends on a large upstream reference file (HPA RNA
tissue consensus, NCBI ``gene_info``).  These helpers download and cache those
files on demand so the scripts are self-contained, durable artifacts — re-run
any of them later and the source is fetched fresh, no manual ``/tmp`` setup.

Cache lives at ``<repo>/.cache/`` (gitignored).  Pass an explicit local path to
any of the scripts' ``--*`` source flags to use a pre-downloaded copy instead.
"""

from __future__ import annotations

import urllib.request
import zipfile
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
_CACHE_DIR = _REPO_ROOT / ".cache"

RNA_TISSUE_CONSENSUS_URL = "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"
NCBI_GENE_INFO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
)


def _download(url: str, dest: Path) -> None:
    _CACHE_DIR.mkdir(exist_ok=True)
    print(f"  downloading {url} -> {dest} ...")
    with urllib.request.urlopen(url, timeout=300) as response:
        dest.write_bytes(response.read())


def rna_tissue_consensus_path(source: str | None = None) -> Path:
    """Return a local path to HPA ``rna_tissue_consensus.tsv``.

    ``source`` may be a local ``.tsv``/``.zip`` path or a URL (default: HPA).
    The zip is extracted into the cache; an existing extract is reused.
    """
    src = source or RNA_TISSUE_CONSENSUS_URL
    if src.endswith(".tsv") and not src.startswith(("http://", "https://")):
        return Path(src)

    tsv_path = _CACHE_DIR / "rna_tissue_consensus.tsv"
    if tsv_path.exists():
        return tsv_path

    if src.startswith(("http://", "https://")):
        zip_path = _CACHE_DIR / "rna_tissue_consensus.tsv.zip"
        if not zip_path.exists():
            _download(src, zip_path)
    else:
        zip_path = Path(src)

    with zipfile.ZipFile(zip_path) as zf:
        member = next(n for n in zf.namelist() if n.endswith(".tsv"))
        _CACHE_DIR.mkdir(exist_ok=True)
        with zf.open(member) as fsrc:
            tsv_path.write_bytes(fsrc.read())
    return tsv_path


def ncbi_gene_info_path(source: str | None = None) -> Path:
    """Return a local path to NCBI human ``gene_info.gz`` (download + cache)."""
    src = source or NCBI_GENE_INFO_URL
    if not src.startswith(("http://", "https://")):
        return Path(src)

    dest = _CACHE_DIR / "Homo_sapiens.gene_info.gz"
    if not dest.exists():
        _download(src, dest)
    return dest
