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

Restriction matching and serotype queries share :func:`parse_mhc`, so
mhcgnomes owns parsing and tsarina caches the results.
``expect`` lets a caller state what the token is supposed to be — the kind of
expectation a CLI flag or a curated paper record carries — which is what
separates the serological ``HLA-A2`` from the molecular ``HLA-A*02:01`` when a
string could be read either way.
"""

from __future__ import annotations

from functools import lru_cache

#: Legacy curated labels retained for existing hitlist indexes even when
#: mhcgnomes cannot parse them. All other labels require a Serotype result;
#: spelling alone cannot distinguish a serotype from an allele like A0201.
_LEGACY_SEROTYPES = frozenset({"DR1B", "DR3A", "DR7A"})

#: Values :func:`parse_mhc`'s ``expect`` parameter recognizes.
_EXPECTED_KINDS = ("serotype", "allele")


@lru_cache(maxsize=1 << 14)
def parse_mhc(value: str, expect: str = "", species: str | None = None):
    """Parse one MHC designation through mhcgnomes, or ``None`` if it will not.

    Never raises. A non-string ``value``, an ``expect`` this function does
    not recognize, or any failure inside mhcgnomes -- including the deferred
    import itself -- all degrade to ``None`` rather than propagating.

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
    species
        A strict species constraint passed straight to mhcgnomes (e.g.
        ``"HLA"``): the result must belong to exactly that species, with no
        fallback to a default or inferred one.  ``None`` accepts any
        species mhcgnomes recognizes.

    Notes
    -----
    Cached on ``(value, expect, species)``: the distinct designation
    vocabulary is a few thousand strings against millions of observation
    rows.
    """
    if not isinstance(value, str):
        return None
    stripped = value.strip()
    if not stripped:
        return None
    if expect and expect not in _EXPECTED_KINDS:
        return None
    try:
        from mhcgnomes import parse

        required = None
        if expect == "serotype":
            from mhcgnomes.serotype import Serotype

            required = [Serotype]
        elif expect == "allele":
            from mhcgnomes.allele import Allele

            required = [Allele]
        return parse(
            stripped,
            species=species,
            required_result_types=required,
            raise_on_error=False,
        )
    except Exception:
        return None


def serotype_key(value: object) -> str | None:
    """Return a comparison key for one serotype designation.

    Serotype names live in mhcgnomes' HLA table, and hitlist stores them into
    the ``serotypes`` column as ``HLA-<name>``.  Passing both a query and a
    stored token through this function makes the two agree by construction,
    across case and the optional ``HLA-`` prefix.

    This keys a single token to itself -- it does not expand a split
    serotype (``A2403``) into its broad parent (``A24``); ``serotype_key
    ("A2403")`` returns ``"A2403"``, not ``"A24"``.  That resolution comes
    from hitlist storing both names in a matching row's ``serotypes`` cell;
    :func:`serotype_keys` set-intersects a query against everything the cell
    actually lists, which is where split/broad matching happens.

    Three legacy curated names (``DR1B``, ``DR3A``, ``DR7A``) remain queryable
    even when mhcgnomes cannot parse them. Every other token must parse as a
    serotype; molecular alleles and unknown labels yield ``None`` so callers
    can report invalid queries instead of quietly matching nothing.
    """
    if not isinstance(value, str):
        return None
    parsed = parse_mhc(value, expect="serotype")
    if parsed is not None:
        return parsed.name.upper()
    bare = strip_hla_prefix(value).upper()
    return bare if bare in _LEGACY_SEROTYPES else None


def strip_hla_prefix(value: str) -> str:
    """Return ``value`` with a leading, case-insensitive ``HLA-`` removed.

    Shared so a second, independently-drifting copy of this stripping rule
    doesn't accumulate elsewhere (:func:`tsarina.spanning._allele_locus`
    uses this one rather than hand-rolling it).
    """
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


def normalize_mhc_restriction(value: object) -> str | None:
    """Return a canonical HLA restriction string when mhcgnomes can parse it.

    Examples
    --------
    ``"A*02:01"`` and ``"HLA-A*02:01"`` both normalize to
    ``"HLA-A*02:01"``.  Unparseable non-empty strings are returned stripped so
    exact string filters still work for unusual restrictions.

    Constrained to the HLA species so this genuinely returns an *HLA*
    restriction string, as documented.  A non-human designation (e.g. a pig
    ``SLA1*01:01`` or a rat ``RT1A*01:01``) is left unchanged rather than
    silently reparsed into that species' own canonical form -- callers on
    the multi-species paths (``--species any`` filtering) rely on getting
    back what they put in for anything this function doesn't own.
    """
    if not isinstance(value, str):
        return None
    stripped = value.strip()
    if not stripped:
        return None
    parsed = parse_mhc(stripped, species="HLA")
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
