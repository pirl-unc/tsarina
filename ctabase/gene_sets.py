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

"""CTA gene set functions returning sets of gene symbols or Ensembl IDs."""

from __future__ import annotations

from .loader import cta_dataframe


def _cta_by_column(
    column: str,
    filtered_only: bool = False,
    exclude_never_expressed: bool = False,
) -> set[str]:
    """Internal: retrieve gene values from the CTA CSV with optional masks."""
    df = cta_dataframe()
    mask = True
    if filtered_only and "filtered" in df.columns:
        mask = df["filtered"].astype(str).str.lower() == "true"
    if exclude_never_expressed and "never_expressed" in df.columns:
        mask = mask & ~(df["never_expressed"].astype(str).str.lower() == "true")
    subset = df[mask] if not isinstance(mask, bool) else df
    result: set[str] = set()
    if column in subset.columns:
        for x in subset[column]:
            if isinstance(x, str):
                result.update(xi.strip() for xi in x.split(";"))
    return result


def _all_by_column(column: str) -> set[str]:
    """Return all values for a column across the full CTA universe."""
    df = cta_dataframe()
    result: set[str] = set()
    if column in df.columns:
        for x in df[column]:
            if isinstance(x, str):
                result.update(xi.strip() for xi in x.split(";"))
    return result


# ── Primary gene set functions ──────────────────────────────────────────────


def CTA_gene_names() -> set[str]:
    """CTA gene symbols: filtered AND expressed (>= 2 nTPM somewhere).

    This is the recommended default for pMHC discovery.  Excludes
    never-expressed CTAs (no protein data + max RNA < 2 nTPM).
    For the full filtered set including never-expressed, use
    ``CTA_filtered_gene_names()``.
    """
    return _cta_by_column("Symbol", filtered_only=True, exclude_never_expressed=True)


def CTA_gene_ids() -> set[str]:
    """CTA Ensembl gene IDs: filtered AND expressed."""
    return _cta_by_column("Ensembl_Gene_ID", filtered_only=True, exclude_never_expressed=True)


def CTA_filtered_gene_names() -> set[str]:
    """All CTA gene symbols that pass the HPA filter (including never-expressed)."""
    return _cta_by_column("Symbol", filtered_only=True)


def CTA_filtered_gene_ids() -> set[str]:
    """All CTA Ensembl gene IDs that pass the HPA filter (including never-expressed)."""
    return _cta_by_column("Ensembl_Gene_ID", filtered_only=True)


def CTA_never_expressed_gene_names() -> set[str]:
    """CTA genes that pass filter but have no meaningful HPA expression.

    No protein data AND max RNA nTPM < 2.  These are in source databases
    but lack positive evidence of tissue restriction from HPA.
    """
    return CTA_filtered_gene_names() - CTA_gene_names()


def CTA_never_expressed_gene_ids() -> set[str]:
    """CTA Ensembl gene IDs that pass filter but have no meaningful HPA expression."""
    return CTA_filtered_gene_ids() - CTA_gene_ids()


def CTA_unfiltered_gene_names() -> set[str]:
    """All CTA gene symbols from all source databases (unfiltered).

    This is the full CTA universe -- use for excluding CTA genes from
    a non-CTA comparison set.  Any gene in this set was identified as
    a candidate CTA by at least one source database.
    """
    return _all_by_column("Symbol")


def CTA_unfiltered_gene_ids() -> set[str]:
    """All CTA Ensembl gene IDs from all source databases (unfiltered)."""
    return _all_by_column("Ensembl_Gene_ID")


def CTA_excluded_gene_names() -> set[str]:
    """CTA genes that FAIL the reproductive-tissue filter.

    These are candidate CTAs with evidence of somatic tissue expression.
    Use this set to exclude from a non-CTA comparison set: they should
    not be in the clean CTA set (they leak into healthy tissue) but also
    should not be in a non-CTA set (they are still CTA candidates).
    """
    return CTA_unfiltered_gene_names() - CTA_filtered_gene_names()


def CTA_excluded_gene_ids() -> set[str]:
    """CTA Ensembl gene IDs that FAIL the reproductive-tissue filter."""
    return CTA_unfiltered_gene_ids() - CTA_filtered_gene_ids()
