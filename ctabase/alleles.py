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

"""Public HLA class I allele panels for population-spanning vaccine design.

Panels are organized as nested tiers: each larger panel is a superset of
the previous one, adding alleles for broader geographic coverage.

Tier structure::

    IEDB-27 (A/B)           27 alleles   Global baseline
    IEDB-36 (A/B/C)         36 alleles   + HLA-C coverage
    Global-44               44 alleles   + East Asia, South Asia, Sub-Saharan Africa
    Global-48               48 alleles   + Latin America, MENA
    Global-51               51 alleles   + additional Sub-Saharan Africa

The IEDB-27 baseline corresponds to the published IEDB/TepiTool panel of
the 27 most frequent global MHC class I alleles.
"""

from __future__ import annotations

from collections import OrderedDict

# ── Allele lists ────────────────────────────────────────────────────────────

IEDB27_AB: list[str] = [
    "HLA-A*01:01",
    "HLA-A*02:01",
    "HLA-A*02:03",
    "HLA-A*02:06",
    "HLA-A*03:01",
    "HLA-A*11:01",
    "HLA-A*23:01",
    "HLA-A*24:02",
    "HLA-A*26:01",
    "HLA-A*30:01",
    "HLA-A*30:02",
    "HLA-A*31:01",
    "HLA-A*32:01",
    "HLA-A*33:01",
    "HLA-A*68:01",
    "HLA-A*68:02",
    "HLA-B*07:02",
    "HLA-B*08:01",
    "HLA-B*15:01",
    "HLA-B*35:01",
    "HLA-B*40:01",
    "HLA-B*44:02",
    "HLA-B*44:03",
    "HLA-B*51:01",
    "HLA-B*53:01",
    "HLA-B*57:01",
    "HLA-B*58:01",
]

HLA_C_ADDON: list[str] = [
    "HLA-C*03:04",
    "HLA-C*04:01",
    "HLA-C*05:01",
    "HLA-C*06:02",
    "HLA-C*07:01",
    "HLA-C*07:02",
    "HLA-C*08:02",
    "HLA-C*12:03",
    "HLA-C*15:02",
]

GLOBAL44_ADDON: list[str] = [
    "HLA-A*02:07",
    "HLA-A*02:11",
    "HLA-A*74:01",
    "HLA-B*46:01",
    "HLA-B*40:06",
    "HLA-B*15:03",
    "HLA-B*58:02",
    "HLA-C*01:02",
]

GLOBAL48_ADDON: list[str] = [
    "HLA-A*29:02",
    "HLA-B*15:02",
    "HLA-B*40:02",
    "HLA-B*50:01",
]

GLOBAL51_SSA_ADDON: list[str] = [
    "HLA-B*57:03",
    "HLA-B*81:01",
    "HLA-C*17:01",
]


def _sorted_alleles(values: set[str]) -> list[str]:
    return sorted(values, key=lambda allele: (allele.split("*", 1)[0], allele))


# ── Panel definitions ───────────────────────────────────────────────────────

PANEL_DEFINITIONS: OrderedDict[str, dict] = OrderedDict(
    [
        (
            "iedb27_ab",
            {
                "label": "IEDB-27 (A/B)",
                "description": (
                    "Published IEDB/TepiTool panel of 27 most frequent global class I alleles."
                ),
                "alleles": IEDB27_AB,
            },
        ),
        (
            "iedb36_abc",
            {
                "label": "IEDB-36 (A/B/C)",
                "description": "IEDB-27 plus a 9-allele HLA-C add-on for broader class I coverage.",
                "alleles": _sorted_alleles(set(IEDB27_AB) | set(HLA_C_ADDON)),
            },
        ),
        (
            "global44_abc",
            {
                "label": "Global-44",
                "description": (
                    "IEDB-36 plus East Asia, South Asia, and Sub-Saharan Africa add-ons."
                ),
                "alleles": _sorted_alleles(set(IEDB27_AB) | set(HLA_C_ADDON) | set(GLOBAL44_ADDON)),
            },
        ),
        (
            "global48_abc",
            {
                "label": "Global-48",
                "description": "Global-44 plus Latin America and MENA-focused add-ons.",
                "alleles": _sorted_alleles(
                    set(IEDB27_AB) | set(HLA_C_ADDON) | set(GLOBAL44_ADDON) | set(GLOBAL48_ADDON)
                ),
            },
        ),
        (
            "global51_abc_ssa",
            {
                "label": "Global-51",
                "description": "Global-48 plus extra Sub-Saharan Africa-enriched alleles.",
                "alleles": _sorted_alleles(
                    set(IEDB27_AB)
                    | set(HLA_C_ADDON)
                    | set(GLOBAL44_ADDON)
                    | set(GLOBAL48_ADDON)
                    | set(GLOBAL51_SSA_ADDON)
                ),
            },
        ),
    ]
)

# ── Per-allele source categories ────────────────────────────────────────────

PANEL_SOURCE_CATEGORIES: dict[str, str] = {
    "HLA-A*01:01": "IEDB-27 baseline",
    "HLA-A*02:01": "IEDB-27 baseline",
    "HLA-A*02:03": "IEDB-27 baseline",
    "HLA-A*02:06": "IEDB-27 baseline",
    "HLA-A*03:01": "IEDB-27 baseline",
    "HLA-A*11:01": "IEDB-27 baseline",
    "HLA-A*23:01": "IEDB-27 baseline",
    "HLA-A*24:02": "IEDB-27 baseline",
    "HLA-A*26:01": "IEDB-27 baseline",
    "HLA-A*30:01": "IEDB-27 baseline",
    "HLA-A*30:02": "IEDB-27 baseline",
    "HLA-A*31:01": "IEDB-27 baseline",
    "HLA-A*32:01": "IEDB-27 baseline",
    "HLA-A*33:01": "IEDB-27 baseline",
    "HLA-A*68:01": "IEDB-27 baseline",
    "HLA-A*68:02": "IEDB-27 baseline",
    "HLA-B*07:02": "IEDB-27 baseline",
    "HLA-B*08:01": "IEDB-27 baseline",
    "HLA-B*15:01": "IEDB-27 baseline",
    "HLA-B*35:01": "IEDB-27 baseline",
    "HLA-B*40:01": "IEDB-27 baseline",
    "HLA-B*44:02": "IEDB-27 baseline",
    "HLA-B*44:03": "IEDB-27 baseline",
    "HLA-B*51:01": "IEDB-27 baseline",
    "HLA-B*53:01": "IEDB-27 baseline",
    "HLA-B*57:01": "IEDB-27 baseline",
    "HLA-B*58:01": "IEDB-27 baseline",
    "HLA-C*03:04": "HLA-C global add-on",
    "HLA-C*04:01": "HLA-C global add-on",
    "HLA-C*05:01": "HLA-C global add-on",
    "HLA-C*06:02": "HLA-C global add-on",
    "HLA-C*07:01": "HLA-C global add-on",
    "HLA-C*07:02": "HLA-C global add-on",
    "HLA-C*08:02": "HLA-C global add-on",
    "HLA-C*12:03": "HLA-C global add-on",
    "HLA-C*15:02": "HLA-C global add-on",
    "HLA-A*02:07": "East Asia / Southeast Asia add-on",
    "HLA-A*02:11": "South Asia add-on",
    "HLA-A*74:01": "Sub-Saharan Africa add-on",
    "HLA-B*46:01": "East Asia / Southeast Asia add-on",
    "HLA-B*40:06": "South Asia / MENA add-on",
    "HLA-B*15:03": "Sub-Saharan Africa add-on",
    "HLA-B*58:02": "Sub-Saharan Africa add-on",
    "HLA-C*01:02": "East Asia / Southeast Asia / Latin America add-on",
    "HLA-A*29:02": "Latin America / IEDB MS add-on",
    "HLA-B*15:02": "Southeast Asia / IEDB MS add-on",
    "HLA-B*40:02": "Latin America / IEDB MS add-on",
    "HLA-B*50:01": "MENA / Arab add-on",
    "HLA-B*57:03": "Sub-Saharan Africa extension",
    "HLA-B*81:01": "Sub-Saharan Africa extension",
    "HLA-C*17:01": "Sub-Saharan Africa / MENA extension",
}


# ── Convenience functions ───────────────────────────────────────────────────


def get_panel(name: str) -> list[str]:
    """Return allele list for a named panel.

    Parameters
    ----------
    name : str
        Panel key: ``"iedb27_ab"``, ``"iedb36_abc"``, ``"global44_abc"``,
        ``"global48_abc"``, or ``"global51_abc_ssa"``.

    Returns
    -------
    list[str]
        Sorted list of HLA allele names (e.g. ``"HLA-A*02:01"``).

    Raises
    ------
    KeyError
        If ``name`` is not a valid panel key.
    """
    return list(PANEL_DEFINITIONS[name]["alleles"])


def panel_names() -> list[str]:
    """Return available panel keys in tier order."""
    return list(PANEL_DEFINITIONS.keys())
