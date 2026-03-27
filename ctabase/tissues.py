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

"""Tissue definitions and HPA confidence thresholds for CTA restriction analysis.

Three tiers of reproductive tissue sets determine how strictly a gene's
expression must be confined to qualify as a cancer-testis antigen:

- **Core**: testis, ovary, placenta -- the classic CT restriction definition.
- **Extended**: core + other reproductive-tract tissues (cervix, endometrium,
  epididymis, fallopian tube, prostate, seminal vesicle, vagina).
- **Permissive**: extended + breast (including lactating breast).

Thymus is excluded from all restriction checks because AIRE-mediated
expression in medullary thymic epithelial cells (mTECs) is expected
for CTAs and does not indicate somatic tissue leakage.

The adaptive protein-RNA thresholds scale the required deflated RNA
reproductive fraction based on how reliable the protein (IHC) evidence
is for a given gene.
"""

from __future__ import annotations

# ── Tissue sets ─────────────────────────────────────────────────────────────

CORE_REPRODUCTIVE_TISSUES: frozenset[str] = frozenset({"ovary", "placenta", "testis"})

EXTENDED_REPRODUCTIVE_TISSUES: frozenset[str] = CORE_REPRODUCTIVE_TISSUES | frozenset(
    {
        "cervix",
        "endometrium",
        "epididymis",
        "fallopian tube",
        "prostate",
        "seminal vesicle",
        "vagina",
    }
)

PERMISSIVE_REPRODUCTIVE_TISSUES: frozenset[str] = EXTENDED_REPRODUCTIVE_TISSUES | frozenset(
    {"breast", "lactating breast"}
)


# ── Adaptive HPA protein → RNA thresholds ──────────────────────────────────

#: Mapping from HPA antibody reliability tier to the minimum deflated
#: RNA reproductive fraction required for a gene to pass the adaptive
#: confidence filter.  Higher-confidence protein data allows a more
#: relaxed RNA threshold.
HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS: dict[str, float] = {
    "Enhanced": 0.80,
    "Supported": 0.90,
    "Approved": 0.95,
    "Uncertain": 0.99,
    "Missing": 0.99,
}

#: Maximum number of tissues with RNA detected (nTPM >= 1) for the
#: strict marker-based filter.
HPA_MARKER_STRICT_MAX_RNA_TISSUES: int = 7

# ── Protein reliability ordering ────────────────────────────────────────────

#: HPA antibody reliability tiers, ordered from strongest to weakest.
PROTEIN_RELIABILITY_ORDER: list[str] = [
    "Enhanced",
    "Supported",
    "Approved",
    "Uncertain",
]


def is_tissue_restricted(
    tissues: set[str] | frozenset[str],
    allowed: frozenset[str] = CORE_REPRODUCTIVE_TISSUES,
    exclude_thymus: bool = True,
) -> bool:
    """Check whether all detected tissues are within an allowed set.

    Parameters
    ----------
    tissues
        Set of tissue names where expression was detected.
    allowed
        Set of tissue names considered acceptable (default: core reproductive).
    exclude_thymus
        If True (default), remove "thymus" from *tissues* before checking.
    """
    check = tissues - {"thymus"} if exclude_thymus else tissues
    return check.issubset(allowed)


def nonreproductive_tissues(
    tissues: set[str] | frozenset[str],
    definition: frozenset[str] = CORE_REPRODUCTIVE_TISSUES,
    exclude_thymus: bool = True,
) -> frozenset[str]:
    """Return tissues that are NOT in the reproductive set.

    Parameters
    ----------
    tissues
        Set of tissue names where expression was detected.
    definition
        Which reproductive definition to use (default: core).
    exclude_thymus
        If True (default), thymus is not counted as non-reproductive.
    """
    excluded = tissues - definition
    if exclude_thymus:
        excluded = excluded - {"thymus"}
    return frozenset(excluded)


def adaptive_rna_threshold(protein_reliability: str) -> float:
    """Return the minimum deflated RNA reproductive fraction for a protein tier.

    Parameters
    ----------
    protein_reliability
        HPA antibody reliability label (Enhanced, Supported, Approved,
        Uncertain) or "Missing" / "no data".

    Returns
    -------
    float
        Minimum deflated reproductive fraction required to pass the
        adaptive confidence filter.
    """
    key = (
        protein_reliability
        if protein_reliability in HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS
        else "Missing"
    )
    return HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS[key]
