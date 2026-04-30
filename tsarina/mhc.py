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

"""Shared MHC restriction normalization and matching helpers."""

from __future__ import annotations


def _parse_hla(value: str):
    from mhcgnomes import parse

    query = value if value.upper().startswith("HLA-") else f"HLA-{value}"
    return parse(query)


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
    try:
        parsed = _parse_hla(stripped)
    except Exception:
        return stripped
    if hasattr(parsed, "to_string"):
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
