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


def _extract_values(df, column: str) -> set[str]:
    """Extract unique values from a column, splitting semicolons."""
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


# ── Axis-based gene set functions ──────────────────────────────────────────


def _axis_filter(column: str, axis_col: str, values: str | set[str]) -> set[str]:
    """Internal: filter genes by a single axis column value(s)."""
    df = cta_dataframe()
    if axis_col not in df.columns:
        return set()
    if isinstance(values, str):
        values = {values}
    mask = df[axis_col].isin(values)
    return _extract_values(df[mask], column)


def CTA_testis_restricted_gene_names() -> set[str]:
    """CTA genes with synthesized restriction = TESTIS.

    Includes genes with IHC testis-only protein AND genes where RNA
    indicates testis-only expression.  The blood-testis barrier prevents
    testicular proteins from reaching circulation.
    """
    return _axis_filter("Symbol", "restriction", "TESTIS")


def CTA_testis_restricted_gene_ids() -> set[str]:
    """Ensembl gene IDs for TESTIS restriction CTAs."""
    return _axis_filter("Ensembl_Gene_ID", "restriction", "TESTIS")


def CTA_placental_restricted_gene_names() -> set[str]:
    """CTA genes with IHC protein in placenta ± testis (no ovary).

    Serum detection in non-pregnant individuals implies cancer.
    """
    return _axis_filter("Symbol", "restriction", "PLACENTAL")


def CTA_placental_restricted_gene_ids() -> set[str]:
    """Ensembl gene IDs for PLACENTAL restriction CTAs."""
    return _axis_filter("Ensembl_Gene_ID", "restriction", "PLACENTAL")


def CTA_by_axes(
    restriction: str | set[str] | None = None,
    protein_restriction: str | set[str] | None = None,
    rna_restriction: str | set[str] | None = None,
    rna_restriction_level: str | set[str] | None = None,
    ms_restriction: str | set[str] | None = None,
    restriction_confidence: str | set[str] | None = None,
    column: str = "Symbol",
) -> set[str]:
    """Return CTA gene identifiers matching specified axis values.

    Parameters
    ----------
    restriction
        Synthesized restriction (None = any).
    protein_restriction
        Per-modality protein restriction (None = any).
    rna_restriction
        Per-modality RNA restriction (None = any).
    rna_restriction_level
        RNA restriction quality: STRICT / MODERATE / PERMISSIVE (None = any).
    ms_restriction
        MS restriction classification (None = any).
    restriction_confidence
        Synthesized confidence: HIGH / MODERATE / LOW (None = any).
    column
        Column to return (``"Symbol"`` or ``"Ensembl_Gene_ID"``).
    """
    df = cta_dataframe()
    mask = True

    for axis_col, values in [
        ("restriction", restriction),
        ("protein_restriction", protein_restriction),
        ("rna_restriction", rna_restriction),
        ("rna_restriction_level", rna_restriction_level),
        ("ms_restriction", ms_restriction),
        ("restriction_confidence", restriction_confidence),
    ]:
        if values is not None and axis_col in df.columns:
            if isinstance(values, str):
                values = {values}
            mask = mask & df[axis_col].isin(values)

    subset = df[mask] if not isinstance(mask, bool) else df
    return _extract_values(subset, column)
