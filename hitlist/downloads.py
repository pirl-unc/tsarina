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

"""Data management for external datasets (IEDB, CEDAR, UniProt proteomes).

Provides a registry for tracking downloaded data files across sessions,
auto-fetching for datasets that can be programmatically downloaded (e.g.
UniProt proteomes), and path resolution for use in analysis pipelines.

Storage location: ``~/.hitlist/`` (override with ``HITLIST_DATA_DIR`` env var).

Typical usage::

    from hitlist.downloads import register, get_path, fetch, list_datasets

    # Register a manually downloaded IEDB export
    register("iedb", "/data/mhc_ligand_full.csv")

    # Later, resolve the path
    path = get_path("iedb")

    # Auto-fetch a UniProt viral proteome
    fetch("hpv16")

    # List everything
    list_datasets()

CLI::

    hitlist data register iedb /data/mhc_ligand_full.csv
    hitlist data register cedar /data/cedar-mhc-ligand-full.csv
    hitlist data fetch hpv16
    hitlist data list
    hitlist data path iedb
"""

from __future__ import annotations

import json
import os
import shutil
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

# ── Data directory ──────────────────────────────────────────────────────────

_DEFAULT_DATA_DIR = Path.home() / ".hitlist"


def data_dir() -> Path:
    """Return the hitlist data directory, creating it if needed."""
    d = Path(os.environ.get("HITLIST_DATA_DIR", str(_DEFAULT_DATA_DIR)))
    d.mkdir(parents=True, exist_ok=True)
    return d


def _manifest_path() -> Path:
    return data_dir() / "manifest.json"


def _load_manifest() -> dict:
    p = _manifest_path()
    if p.exists():
        return json.loads(p.read_text())
    return {"datasets": {}}


def _save_manifest(manifest: dict) -> None:
    p = _manifest_path()
    p.write_text(json.dumps(manifest, indent=2, default=str) + "\n")


# ── Known fetchable datasets ───────────────────────────────────────────────

#: Datasets that can be auto-downloaded.  Keys are dataset names,
#: values have ``url`` (download URL), ``filename`` (local name),
#: and ``description``.
FETCHABLE_DATASETS: dict[str, dict[str, str]] = {
    "hpv16": {
        "url": "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000006729&format=fasta",
        "filename": "hpv16.fasta",
        "description": "HPV-16 proteome (UniProt UP000006729)",
    },
    "hpv18": {
        "url": "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000006728&format=fasta",
        "filename": "hpv18.fasta",
        "description": "HPV-18 proteome (UniProt UP000006728)",
    },
    "ebv": {
        "url": "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000153037&format=fasta",
        "filename": "ebv.fasta",
        "description": "EBV/HHV-4 proteome (UniProt UP000153037)",
    },
    "htlv1": {
        "url": "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000002063&format=fasta",
        "filename": "htlv1.fasta",
        "description": "HTLV-1 proteome (UniProt UP000002063)",
    },
    "hbv": {
        "url": "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000126453&format=fasta",
        "filename": "hbv.fasta",
        "description": "HBV proteome (UniProt UP000126453)",
    },
    "hcv": {
        "url": "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000000518&format=fasta",
        "filename": "hcv.fasta",
        "description": "HCV proteome (UniProt UP000000518)",
    },
    "kshv": {
        "url": "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000009113&format=fasta",
        "filename": "kshv.fasta",
        "description": "KSHV/HHV-8 proteome (UniProt UP000009113)",
    },
    "mcpyv": {
        "url": "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000116695&format=fasta",
        "filename": "mcpyv.fasta",
        "description": "MCPyV proteome (UniProt UP000116695)",
    },
    "hiv1": {
        "url": "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000002241&format=fasta",
        "filename": "hiv1.fasta",
        "description": "HIV-1 proteome (UniProt UP000002241)",
    },
}

#: Datasets requiring manual download (terms-of-use gated).
MANUAL_DATASETS: dict[str, dict[str, str]] = {
    "iedb": {
        "download_url": "https://www.iedb.org/downloader.php?file_name=doc/mhc_ligand_full_single_file.zip",
        "description": "IEDB MHC ligand full export (requires IEDB terms acceptance)",
        "expected_filename": "mhc_ligand_full.csv",
    },
    "cedar": {
        "download_url": "https://cedar.iedb.org/downloader.php?file_name=doc/cedar_mhc_ligand_full.zip",
        "description": "CEDAR MHC ligand full export",
        "expected_filename": "cedar-mhc-ligand-full.csv",
    },
}


# ── Core API ────────────────────────────────────────────────────────────────


def register(name: str, path: str | Path, description: str | None = None) -> Path:
    """Register a local file path for a named dataset.

    Parameters
    ----------
    name
        Dataset name (e.g. ``"iedb"``, ``"cedar"``).
    path
        Path to the file on disk.
    description
        Optional human-readable description.

    Returns
    -------
    Path
        The resolved absolute path.

    Raises
    ------
    FileNotFoundError
        If the path does not exist.
    """
    p = Path(path).resolve()
    if not p.exists():
        raise FileNotFoundError(f"File not found: {p}")

    manifest = _load_manifest()
    desc = description
    if desc is None and name in MANUAL_DATASETS:
        desc = MANUAL_DATASETS[name]["description"]
    if desc is None and name in FETCHABLE_DATASETS:
        desc = FETCHABLE_DATASETS[name]["description"]

    manifest["datasets"][name] = {
        "path": str(p),
        "registered": datetime.now(timezone.utc).isoformat(),
        "size_bytes": p.stat().st_size,
        "description": desc or "",
        "source": "registered",
    }
    _save_manifest(manifest)
    return p


def fetch(name: str, force: bool = False) -> Path:
    """Download a fetchable dataset (e.g. UniProt viral proteome).

    Parameters
    ----------
    name
        Dataset name (must be a key in :data:`FETCHABLE_DATASETS`).
    force
        If True, re-download even if already present.

    Returns
    -------
    Path
        Path to the downloaded file.

    Raises
    ------
    ValueError
        If ``name`` is not a fetchable dataset.
    """
    if name not in FETCHABLE_DATASETS:
        if name in MANUAL_DATASETS:
            info = MANUAL_DATASETS[name]
            raise ValueError(
                f"'{name}' requires manual download from:\n"
                f"  {info['download_url']}\n"
                f"Then register with: hitlist data register {name} /path/to/{info['expected_filename']}"
            )
        available = sorted(set(FETCHABLE_DATASETS) | set(MANUAL_DATASETS))
        raise ValueError(f"Unknown dataset '{name}'. Available: {available}")

    info = FETCHABLE_DATASETS[name]
    dest = data_dir() / info["filename"]

    if dest.exists() and not force:
        # Already downloaded, ensure registered
        manifest = _load_manifest()
        if name not in manifest.get("datasets", {}):
            register(name, dest, info["description"])
        return dest

    print(f"Downloading {info['description']}...")
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    try:
        urllib.request.urlretrieve(info["url"], str(tmp))
        shutil.move(str(tmp), str(dest))
    finally:
        if tmp.exists():
            tmp.unlink()

    manifest = _load_manifest()
    manifest["datasets"][name] = {
        "path": str(dest),
        "registered": datetime.now(timezone.utc).isoformat(),
        "size_bytes": dest.stat().st_size,
        "description": info["description"],
        "source": info["url"],
    }
    _save_manifest(manifest)
    print(f"  Saved to {dest} ({dest.stat().st_size:,} bytes)")
    return dest


def get_path(name: str) -> Path:
    """Resolve a dataset name to its local file path.

    Parameters
    ----------
    name
        Dataset name.

    Returns
    -------
    Path
        Path to the file.

    Raises
    ------
    KeyError
        If the dataset is not registered or fetched.
    FileNotFoundError
        If the registered path no longer exists.
    """
    manifest = _load_manifest()
    entry = manifest.get("datasets", {}).get(name)
    if entry is None:
        hint = ""
        if name in FETCHABLE_DATASETS:
            hint = f"\n  Fetch with: hitlist data fetch {name}"
        elif name in MANUAL_DATASETS:
            info = MANUAL_DATASETS[name]
            hint = (
                f"\n  Download from: {info['download_url']}"
                f"\n  Then register: hitlist data register {name} /path/to/file"
            )
        raise KeyError(f"Dataset '{name}' not registered.{hint}")

    p = Path(entry["path"])
    if not p.exists():
        raise FileNotFoundError(
            f"Registered path for '{name}' no longer exists: {p}\nRe-register or re-fetch."
        )
    return p


def list_datasets() -> dict[str, dict]:
    """Return all registered/fetched datasets.

    Returns
    -------
    dict[str, dict]
        Mapping from dataset name to metadata (path, size, date, etc.).
    """
    return dict(_load_manifest().get("datasets", {}))


def available_datasets() -> dict[str, str]:
    """Return all known dataset names with descriptions.

    Returns
    -------
    dict[str, str]
        Mapping from name to description, including both fetchable
        and manual-download datasets.
    """
    result = {}
    for name, info in FETCHABLE_DATASETS.items():
        result[name] = info["description"] + " [auto-fetch]"
    for name, info in MANUAL_DATASETS.items():
        result[name] = info["description"] + " [manual download]"
    return result


def remove(name: str) -> None:
    """Remove a dataset from the registry (does not delete the file).

    Parameters
    ----------
    name
        Dataset name to unregister.
    """
    manifest = _load_manifest()
    manifest.get("datasets", {}).pop(name, None)
    _save_manifest(manifest)
