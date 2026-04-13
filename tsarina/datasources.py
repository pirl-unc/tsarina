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

"""Path resolution for IEDB/CEDAR via the hitlist data registry.

Callers pass ``None`` (the default) to let tsarina pick up whatever the user
registered with ``tsarina data register iedb /path/to/mhc_ligand_full.csv`` or
``hitlist data register iedb ...``.  Explicit paths are honored unchanged.

CEDAR is optional throughout tsarina; a missing CEDAR registration resolves
silently to ``None`` rather than raising.  IEDB is required and surfaces an
actionable error message pointing at the registration command.
"""

from __future__ import annotations

from pathlib import Path

from hitlist.downloads import get_path


class DatasetNotRegisteredError(FileNotFoundError):
    """Raised when a required dataset is not in the hitlist registry."""


def _register_hint(name: str) -> str:
    return (
        f"Dataset '{name}' is not registered with hitlist.  Download the "
        f"{name.upper()} MHC ligand full export from https://www.iedb.org/ "
        f"(or https://cedar.iedb.org/) and register it with:\n"
        f"    tsarina data register {name} /path/to/{name}-mhc-ligand-full.csv"
    )


def resolve_iedb_path(iedb_path: str | Path | None) -> Path:
    """Resolve an IEDB path, falling back to the hitlist registry.

    Parameters
    ----------
    iedb_path
        Explicit path, or ``None`` to auto-resolve.

    Returns
    -------
    Path
        Resolved IEDB CSV path.

    Raises
    ------
    DatasetNotRegisteredError
        If ``iedb_path`` is None and IEDB is not registered.
    """
    if iedb_path is not None:
        return Path(iedb_path)
    try:
        return Path(get_path("iedb"))
    except (KeyError, FileNotFoundError) as e:
        raise DatasetNotRegisteredError(_register_hint("iedb")) from e


def resolve_cedar_path(cedar_path: str | Path | None) -> Path | None:
    """Resolve a CEDAR path, silently returning None if unregistered.

    CEDAR is optional throughout tsarina; a missing registration is not an
    error.  Pass an explicit path to force a specific file.
    """
    if cedar_path is not None:
        return Path(cedar_path)
    try:
        return Path(get_path("cedar"))
    except (KeyError, FileNotFoundError):
        return None


def resolve_dataset_paths(
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    require_iedb: bool = True,
) -> tuple[Path | None, Path | None]:
    """Resolve both IEDB and CEDAR paths from the hitlist registry.

    Parameters
    ----------
    iedb_path, cedar_path
        Explicit overrides. ``None`` means auto-resolve.
    require_iedb
        If True (default), raise if IEDB is unregistered.  If False, IEDB
        resolves silently to None (useful for ``--skip-ms-evidence``-style
        call paths).

    Returns
    -------
    (iedb, cedar)
        Both ``Path`` or ``None``.
    """
    if require_iedb:
        iedb = resolve_iedb_path(iedb_path)
    else:
        try:
            iedb = resolve_iedb_path(iedb_path)
        except DatasetNotRegisteredError:
            iedb = None
    cedar = resolve_cedar_path(cedar_path)
    return iedb, cedar
