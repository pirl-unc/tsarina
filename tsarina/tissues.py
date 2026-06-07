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
    "Uncertain": 0.97,
    "Missing": 0.97,
}

#: Maximum number of tissues with RNA detected (nTPM >= 1) for the
#: strict marker-based filter.
HPA_MARKER_STRICT_MAX_RNA_TISSUES: int = 7

#: RNA expression floor (nTPM) for the ``never_expressed`` flag.  A gene with
#: no HPA protein (IHC) data and a maximum RNA nTPM below this value across all
#: tissues is flagged ``never_expressed`` -- it passes the reproductive-
#: restriction filter only because the ``+1`` deflation pseudocount yields a
#: 1.0 fraction when every tissue is below 1 nTPM, but HPA lacks the signal to
#: confirm restriction.  Kept explicit/parameterized rather than baked as a
#: magic number into the bundled table.  See tsarina#78.
HPA_EXPRESSION_FLOOR_NTPM: float = 2.0

#: Cancer-testis antigens manually rescued into the expressed CTA set despite
#: being flagged ``never_expressed`` (max RNA < ``HPA_EXPRESSION_FLOOR_NTPM``,
#: no protein data).  These are borderline-but-real CTAs where the low HPA
#: signal reflects detection limits rather than true absence -- corroborated by
#: source-database membership and/or tumor MS evidence.  Preferred over a
#: blanket lowering of the global floor, which would admit paralog
#: cross-mapping noise.  Keyed by (unversioned) Ensembl gene ID.  See
#: tsarina#78.
MANUALLY_EXPRESSED_CTA: frozenset[str] = frozenset(
    {
        "ENSG00000171405",  # XAGE5 -- testis 1.1 nTPM; CTpedia/CTexploreR/daSilva2017
    }
)

#: Genes that entered the CTA universe via a source database but are NOT cancer-
#: testis antigens: conserved/multicopy housekeeping families (core histones,
#: alpha-tubulins) and placental chorionic gonadotropin.  These can pass the
#: reproductive-restriction filter on a testis-enriched copy, yet they make poor
#: targets (ubiquitous/essential proteins) and pollute any CTA-keyed sequence or
#: expression analysis.  Dropped from the universe at load time.  Keyed by
#: (unversioned) Ensembl gene ID.  See tsarina#92.  Note: testis-specific
#: histone variants (e.g. H1-6 / HIST1H1T) are deliberately NOT excluded.
NON_CTA_EXCLUDED_GENE_IDS: frozenset[str] = frozenset(
    {
        "ENSG00000274618",  # H4C6   -- core histone H4 (conserved/multicopy)
        "ENSG00000146047",  # H2BC1  -- core histone H2B
        "ENSG00000276410",  # H2BC3  -- core histone H2B
        "ENSG00000124610",  # H1-1   -- somatic linker histone H1.1 (HIST1H1A)
        "ENSG00000213030",  # CGB8   -- chorionic gonadotropin beta (placental)
        "ENSG00000198033",  # TUBA3C -- alpha-tubulin (conserved/essential)
        "ENSG00000152086",  # TUBA3E -- alpha-tubulin
    }
)

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
    # Normalize whitespace/case so a non-canonical label (e.g. " enhanced ")
    # maps to its tier instead of silently falling back to the strict "Missing"
    # threshold.  Unknown labels (incl. "no data") use the "Missing" floor.
    normalized = str(protein_reliability).strip().casefold()
    by_casefold = {k.casefold(): k for k in HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS}
    key = by_casefold.get(normalized, "Missing")
    return HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS[key]
