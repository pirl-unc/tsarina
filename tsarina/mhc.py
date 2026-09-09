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

"""Shared MHC restriction normalization and matching helpers.

Every designation tsarina reads or accepts goes through :func:`parse_mhc`, so
mhcgnomes owns what a token means and tsarina owns nothing but the caching.
``expect`` lets a caller state what the token is supposed to be — the kind of
expectation a CLI flag or a curated paper record carries — which is what
separates the serological ``HLA-A2`` from the molecular ``HLA-A*02:01`` when a
string could be read either way.
"""

from __future__ import annotations

import re
from functools import lru_cache

#: Shape of a serological designation: a locus prefix, a digit, then optional
#: split / public-epitope characters (``A2``, ``A2.1``, ``Bw4``, ``DR1B``).
#: Guards the fallback in :func:`serotype_key` so an allele or a class-only
#: string is rejected rather than accepted as an unknown serotype.
_SEROTYPE_SHAPE = re.compile(r"^[A-Za-z]{1,3}[0-9][0-9A-Za-z.]*$")

#: What a token is expected to be, mapped to the mhcgnomes result type that
#: :func:`parse_mhc` will require.
_EXPECTED_RESULT_TYPES = {
    "serotype": ("mhcgnomes.serotype", "Serotype"),
    "allele": ("mhcgnomes.allele", "Allele"),
}


@lru_cache(maxsize=1 << 14)
def parse_mhc(value: str, expect: str = ""):
    """Parse one MHC designation through mhcgnomes, or ``None`` if it will not.

    Parameters
    ----------
    value
        Any designation mhcgnomes accepts: ``A2``, ``HLA-A2``, ``hla-a2``,
        ``Bw4``, ``A*02:01``, ``HLA-DRB1*04:01``.  Case and the ``HLA-``
        prefix are mhcgnomes' problem, not the caller's.
    expect
        ``"serotype"`` or ``"allele"`` to require that reading, so a stated
        expectation wins over mhcgnomes' own preference; ``""`` accepts
        whatever the token parses as.  A token that cannot satisfy the
        requirement yields ``None`` rather than a different kind of answer.

    Notes
    -----
    Cached on ``(value, expect)``: the distinct designation vocabulary is a few
    thousand strings against millions of observation rows.
    """
    from importlib import import_module

    from mhcgnomes import parse

    stripped = value.strip()
    if not stripped:
        return None
    required = None
    if expect:
        module, attr = _EXPECTED_RESULT_TYPES[expect]
        required = [getattr(import_module(module), attr)]
    try:
        return parse(stripped, required_result_types=required, raise_on_error=False)
    except Exception:
        return None


def serotype_key(value: object) -> str | None:
    """Return a comparison key for one serotype designation.

    Serotype names live in mhcgnomes' HLA table, and hitlist stores them into
    the ``serotypes`` column as ``HLA-<name>``.  Passing both a query and a
    stored token through this function makes the two agree by construction,
    across case, the optional ``HLA-`` prefix, and split serotypes (``A2403``
    resolves alongside its broad ``A24``).

    A handful of curated names mhcgnomes writes but will not read back
    (``DR1B``, ``DR3A``, ``DR7A``) fall back to the token as given, so a stored
    serotype stays queryable rather than silently matching nothing.  Anything
    that is not serotype-shaped — an allele, a class-only string — yields
    ``None``, so it can be reported instead of quietly matching nothing.
    """
    if not isinstance(value, str):
        return None
    parsed = parse_mhc(value, expect="serotype")
    if parsed is not None:
        return parsed.name.upper()
    bare = _bare_serotype_name(value)
    return bare.upper() if _SEROTYPE_SHAPE.match(bare) else None


def _bare_serotype_name(value: str) -> str:
    stripped = value.strip()
    return stripped[4:] if stripped[:4].upper() == "HLA-" else stripped


def serotype_keys(cell: object) -> frozenset[str]:
    """Comparison keys for a semicolon-joined ``serotypes`` cell."""
    if not isinstance(cell, str) or not cell:
        return frozenset()
    return _serotype_keys_cached(cell)


@lru_cache(maxsize=1 << 14)
def _serotype_keys_cached(cell: str) -> frozenset[str]:
    keys = (serotype_key(token) for token in cell.split(";"))
    return frozenset(key for key in keys if key)


def _parse_hla(value: str):
    return parse_mhc(value)


def normalize_mhc_restriction(value: object) -> str | None:
    """Return a canonical HLA restriction string when mhcgnomes can parse it.

    Examples
    --------
    ``"A*02:01"`` and ``"HLA-A*02:01"`` both normalize to
    ``"HLA-A*02:01"``.  Unparseable non-empty strings are returned stripped so
    exact string filters still work for unusual restrictions.
    """
    if not isinstance(value, str):
        return None
    stripped = value.strip()
    if not stripped:
        return None
    parsed = _parse_hla(stripped)
    if parsed is not None and hasattr(parsed, "to_string"):
        return parsed.to_string()
    return stripped


def split_mhc_restrictions(value: object) -> tuple[str, ...]:
    """Split a semicolon-joined restriction cell into normalized tokens."""
    if not isinstance(value, str):
        return ()
    normalized = []
    for token in value.split(";"):
        restriction = normalize_mhc_restriction(token)
        if restriction is not None:
            normalized.append(restriction)
    return tuple(normalized)


def normalize_mhc_restriction_set(values: list[str] | tuple[str, ...] | set[str]) -> set[str]:
    """Normalize an iterable of user-supplied MHC restrictions."""
    normalized = set()
    for value in values:
        restriction = normalize_mhc_restriction(value)
        if restriction is not None:
            normalized.add(restriction)
    return normalized


def mhc_restriction_matches_any(value: object, wanted: set[str]) -> bool:
    """Return True when a restriction cell contains any normalized target."""
    if not wanted:
        return True
    return any(restriction in wanted for restriction in split_mhc_restrictions(value))
