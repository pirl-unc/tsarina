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

"""CTA helpers.

All CTA membership, specificity, alias, and HPA-axis definitions are direct
aliases of :mod:`oncoref.cta`. Tsarina adds only the optional ``ms_restriction``
filter to ``CTA_by_axes`` over its downstream gene-level MS evidence.
"""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd
from oncoref import cta as _oncoref_cta

from .loader import cta_dataframe

# Direct aliases make ownership executable, not just documentary: the
# foundational Tsarina API objects are the oncoref implementations themselves.
CTA_gene_names = _oncoref_cta.cta_gene_names
CTA_gene_ids = _oncoref_cta.cta_gene_ids
CTA_filtered_gene_names = _oncoref_cta.cta_filtered_gene_names
CTA_filtered_gene_ids = _oncoref_cta.cta_filtered_gene_ids
CTA_never_expressed_gene_names = _oncoref_cta.cta_never_expressed_gene_names
CTA_never_expressed_gene_ids = _oncoref_cta.cta_never_expressed_gene_ids
CTA_unfiltered_gene_names = _oncoref_cta.cta_unfiltered_gene_names
CTA_unfiltered_gene_ids = _oncoref_cta.cta_unfiltered_gene_ids
CTA_excluded_gene_names = _oncoref_cta.cta_excluded_gene_names
CTA_excluded_gene_ids = _oncoref_cta.cta_excluded_gene_ids
CTA_relaxed_reproductive_gene_names = _oncoref_cta.cta_relaxed_reproductive_gene_names
CTA_relaxed_reproductive_gene_ids = _oncoref_cta.cta_relaxed_reproductive_gene_ids
CTA_testis_restricted_gene_names = _oncoref_cta.cta_testis_restricted_gene_names
CTA_testis_restricted_gene_ids = _oncoref_cta.cta_testis_restricted_gene_ids
CTA_placental_restricted_gene_names = _oncoref_cta.cta_placental_restricted_gene_names
CTA_placental_restricted_gene_ids = _oncoref_cta.cta_placental_restricted_gene_ids
cta_symbol_for_alias = _oncoref_cta.cta_symbol_for_alias


def _filter_values(values: str | Iterable[str] | None) -> set[str] | None:
    if values is None:
        return None
    if isinstance(values, str):
        values = {values}
    return {str(value).upper() for value in values}


def _extract_values(df: pd.DataFrame, column: str) -> set[str]:
    result: set[str] = set()
    if column not in df.columns:
        return result
    for value in df[column]:
        if isinstance(value, str):
            result.update(part.strip() for part in value.split(";") if part.strip())
    return result


def CTA_by_axes(
    *,
    restriction: str | Iterable[str] | None = None,
    protein_restriction: str | Iterable[str] | None = None,
    rna_restriction: str | Iterable[str] | None = None,
    rna_restriction_level: str | Iterable[str] | None = None,
    ms_restriction: str | Iterable[str] | None = None,
    restriction_confidence: str | Iterable[str] | None = None,
    column: str = "Symbol",
    filtered_only: bool = True,
) -> set[str]:
    """Return oncoref CTAs matching HPA axes and optional Tsarina MS evidence.

    The row universe and canonical filtered tier come directly from oncoref.
    Tsarina's only extension is the explicit ``ms_restriction`` predicate.
    """
    df = cta_dataframe().copy()
    if column not in df.columns or "Ensembl_Gene_ID" not in df.columns:
        return set()

    allowed_ids = (
        _oncoref_cta.cta_filtered_gene_ids()
        if filtered_only
        else _oncoref_cta.cta_unfiltered_gene_ids()
    )
    gene_ids = df["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
    mask = gene_ids.isin(allowed_ids)

    for axis_col, values in (
        ("restriction", restriction),
        ("protein_restriction", protein_restriction),
        ("rna_restriction", rna_restriction),
        ("rna_restriction_level", rna_restriction_level),
        ("ms_restriction", ms_restriction),
        ("restriction_confidence", restriction_confidence),
    ):
        wanted = _filter_values(values)
        if wanted is None:
            continue
        if axis_col not in df.columns:
            return set()
        mask &= df[axis_col].astype(str).str.upper().isin(wanted)

    return _extract_values(df[mask], column)
