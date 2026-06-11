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

from functools import lru_cache

import pandas as pd

from .loader import cta_dataframe, passes_filters_mask
from .tissues import MANUALLY_EXPRESSED_CTA

#: ``protein_reliability`` spellings that mean "no HPA IHC data for this gene".
_NO_PROTEIN_RELIABILITY = {"no data", "nan", ""}


def _normalize_alias(name: object) -> str:
    """Collapse a gene name/alias to its comparison key (uppercase alphanumeric).

    So ``NY-ESO-1``, ``ny eso 1`` and ``NYESO1`` all compare equal.
    """
    return "".join(ch for ch in str(name).upper() if ch.isalnum())


#: Canonical resolution for synonyms NCBI shares across near-identical paralogs,
#: where the colloquial name conventionally denotes one specific gene.  Keyed by
#: normalized alias.  Highest precedence in :func:`_alias_to_symbol`.
_CANONICAL_ALIAS_OVERRIDES: dict[str, str] = {
    "NYESO1": "CTAG1B",  # NY-ESO-1 is the cloned CTAG1B; CTAG1A is the paralog.
    "ESO1": "CTAG1B",
}


@lru_cache(maxsize=1)
def _alias_to_symbol() -> dict[str, str]:
    """Map normalized alias/synonym -> official CTA Symbol.

    Precedence: curated canonical overrides > official symbols > first synonym
    in table order.  Built from the ``Symbol`` and ``Aliases`` columns.
    """
    df = cta_dataframe()
    symbols = set(df["Symbol"].astype(str))
    mapping: dict[str, str] = {}
    if "Aliases" in df.columns:
        for symbol, aliases in zip(df["Symbol"], df["Aliases"]):
            if not isinstance(aliases, str):
                continue
            for alias in aliases.split(";"):
                key = _normalize_alias(alias)
                if key and key not in mapping:
                    mapping[key] = str(symbol)
    # Official symbols take precedence over any synonym collision.
    for symbol in df["Symbol"]:
        key = _normalize_alias(symbol)
        if key:
            mapping[key] = str(symbol)
    # Curated canonical overrides win, but only for symbols actually present.
    for key, symbol in _CANONICAL_ALIAS_OVERRIDES.items():
        if symbol in symbols:
            mapping[key] = symbol
    return mapping


def cta_symbol_for_alias(name: str) -> str | None:
    """Resolve a CTA gene name or synonym to its official Symbol.

    Case- and punctuation-insensitive, so ``"NY-ESO-1"``, ``"ESO1"`` and
    ``"CTAG1B"`` all resolve to ``"CTAG1B"``.  Returns ``None`` if *name* is
    not a known CTA symbol or synonym.  See tsarina#77.

    Examples
    --------
    >>> cta_symbol_for_alias("NY-ESO-1")
    'CTAG1B'
    >>> cta_symbol_for_alias("not-a-gene") is None
    True
    """
    return _alias_to_symbol().get(_normalize_alias(name))


def _cta_by_column(
    column: str,
    filtered_only: bool = False,
    exclude_never_expressed: bool = False,
) -> set[str]:
    """Internal: retrieve gene values from the CTA CSV with optional masks."""
    df = cta_dataframe()
    mask = True
    if filtered_only:
        mask = passes_filters_mask(df)
    if exclude_never_expressed and "never_expressed" in df.columns:
        never_expressed = df["never_expressed"].astype(str).str.lower() == "true"
        # Keep manually rescued borderline-but-real CTAs (e.g. XAGE5) in the
        # expressed set even though HPA flags them never_expressed.  See
        # tsarina#78 and tsarina.tissues.MANUALLY_EXPRESSED_CTA.
        if "Ensembl_Gene_ID" in df.columns:
            rescued = (
                df["Ensembl_Gene_ID"].astype(str).str.split(".").str[0].isin(MANUALLY_EXPRESSED_CTA)
            )
            never_expressed = never_expressed & ~rescued
        mask = mask & ~never_expressed
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
    """CTA gene symbols: passes filters AND expressed (>= 2 nTPM somewhere).

    This is the recommended default for pMHC discovery.  Excludes
    never-expressed CTAs (no protein data + max RNA < 2 nTPM).
    For the full filter-passing set including never-expressed, use
    ``CTA_filtered_gene_names()``.
    """
    return _cta_by_column("Symbol", filtered_only=True, exclude_never_expressed=True)


def CTA_gene_ids() -> set[str]:
    """CTA Ensembl gene IDs: passes filters AND expressed."""
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


def _relaxed_reproductive_mask(df: pd.DataFrame, min_deflated_frac: float) -> pd.Series:
    """Mask for the opt-in relaxed reproductive tier (tsarina#119).

    RNA-only genes (no HPA protein data, so the filter failure is purely the RNA
    fraction rather than somatic protein detection) that FAIL the default gate
    but whose deflated reproductive fraction is still >= *min_deflated_frac*.
    """
    failed = ~passes_filters_mask(df)
    rna_only = (
        df["protein_reliability"].astype(str).str.strip().str.lower().isin(_NO_PROTEIN_RELIABILITY)
    )
    frac = pd.to_numeric(df["rna_deflated_reproductive_frac"], errors="coerce")
    return failed & rna_only & (frac >= float(min_deflated_frac))


def CTA_relaxed_reproductive_gene_names(min_deflated_frac: float = 0.80) -> set[str]:
    """Opt-in less-stringent tier: reproductive-dominant CTAs with somatic leak.

    RNA-only candidates (no HPA protein data) that fail the **default**
    reproductive-restriction gate but keep a dominant testis/placenta signal --
    a deflated reproductive fraction of at least *min_deflated_frac* (default
    0.80), with some somatic leakage. Restricting to RNA-only genes means the
    failure is the RNA fraction, not somatic protein detection, so relaxing the
    fraction threshold is meaningful.

    Includes the placental ERV envelopes syncytin-1 (``ERVW-1``) and syncytin-2
    (``ERVFRD-1``) -- well-known onco-placental antigens -- and ``ERVV-1``.
    Provided for consumers willing to accept more somatic leak without loosening
    the default sets (tsarina#119). Disjoint from ``CTA_filtered_gene_names()``.
    """
    df = cta_dataframe()
    return set(df.loc[_relaxed_reproductive_mask(df, min_deflated_frac), "Symbol"])


def CTA_relaxed_reproductive_gene_ids(min_deflated_frac: float = 0.80) -> set[str]:
    """Ensembl gene IDs for :func:`CTA_relaxed_reproductive_gene_names`."""
    df = cta_dataframe()
    return set(df.loc[_relaxed_reproductive_mask(df, min_deflated_frac), "Ensembl_Gene_ID"])


# ── Axis-based gene set functions ──────────────────────────────────────────


def _axis_filter(
    column: str, axis_col: str, values: str | set[str], filtered_only: bool = True
) -> set[str]:
    """Internal: filter genes by a single axis column value(s)."""
    df = cta_dataframe()
    if axis_col not in df.columns:
        return set()
    if isinstance(values, str):
        values = {values}
    mask = df[axis_col].isin(values)
    if filtered_only:
        mask = mask & passes_filters_mask(df)
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
    filtered_only: bool = True,
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
    filtered_only
        If True (default), restrict to genes passing the HPA filter.
        Set to False to query the full 358-gene CTA universe.
    """
    df = cta_dataframe()
    mask = True

    if filtered_only:
        mask = passes_filters_mask(df)

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
