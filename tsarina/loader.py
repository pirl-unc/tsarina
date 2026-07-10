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

from __future__ import annotations

from functools import lru_cache
from os.path import dirname, join

import pandas as pd

_DATA_DIR = join(dirname(__file__), "data")
PASSES_FILTERS_COLUMN = "passes_filters"
LEGACY_FILTERED_COLUMN = "filtered"
SPECIFICITY_ACTION_COLUMN = "specificity_action"
SPECIFICITY_STATUS_COLUMN = "specificity_status"

_MS_EVIDENCE_COLUMNS = (
    "ms_restriction",
    "ms_healthy_somatic_tissues",
    "ms_pmids",
)
_MS_EVIDENCE_DEFAULTS = {
    "ms_restriction": "NO_MS_DATA",
    "ms_healthy_somatic_tissues": "",
    "ms_pmids": "",
}

_TSARINA_ONLY_EXCLUDED_CTA_ROWS = {
    "ENSG00000187475": {
        "specificity_status": "excluded_structural_gene",
        "specificity_action": "exclude_default",
        "specificity_source_anchor": "tsarina#141:oncoref_histone_family_exclusion",
        "specificity_rationale": (
            "Excluded from the canonical oncoref CTA set by histone-family filtering; "
            "retained here as a tsarina-only evidence candidate because its bundled "
            "MS evidence is RECURRENT_HEALTHY across healthy somatic tissues."
        ),
    }
}


def _unversioned_gene_ids(values: pd.Series) -> pd.Series:
    return values.astype(str).str.split(".").str[0]


def _local_cta_dataframe() -> pd.DataFrame:
    """Return the local bundled table used for tsarina-only evidence columns."""
    from .tissues import NON_CTA_EXCLUDED_GENE_IDS

    path = join(_DATA_DIR, "cancer-testis-antigens.csv")
    df = pd.read_csv(path)
    # Keep this local view aligned with the previous tsarina evidence universe:
    # old hard-excluded non-CTA structural genes stay out, but H1-6 remains so
    # tsarina can preserve its MS safety evidence as an excluded candidate.
    if "Ensembl_Gene_ID" in df.columns and NON_CTA_EXCLUDED_GENE_IDS:
        unversioned = _unversioned_gene_ids(df["Ensembl_Gene_ID"])
        df = df[~unversioned.isin(NON_CTA_EXCLUDED_GENE_IDS)].reset_index(drop=True)
    return df


def _oncoref_cta_dataframe() -> pd.DataFrame:
    """Return oncoref's canonical CTA evidence frame."""
    import oncoref

    return oncoref.cta_evidence().copy()


def _attach_local_ms_evidence(df: pd.DataFrame, local: pd.DataFrame) -> pd.DataFrame:
    """Left-join tsarina's bundled MS safety columns onto oncoref CTA evidence."""
    if "Ensembl_Gene_ID" not in df.columns or "Ensembl_Gene_ID" not in local.columns:
        return df

    ms_columns = [
        column
        for column in _MS_EVIDENCE_COLUMNS
        if column in local.columns and column not in df.columns
    ]
    if not ms_columns:
        return df

    join_key = "_tsarina_gene_id"
    out = df.copy()
    out[join_key] = _unversioned_gene_ids(out["Ensembl_Gene_ID"])
    ms = local[[*ms_columns, "Ensembl_Gene_ID"]].copy()
    ms[join_key] = _unversioned_gene_ids(ms["Ensembl_Gene_ID"])
    ms = ms.drop(columns=["Ensembl_Gene_ID"]).drop_duplicates(subset=[join_key], keep="first")
    out = out.merge(ms, on=join_key, how="left").drop(columns=[join_key])
    for column in ms_columns:
        out[column] = out[column].fillna(_MS_EVIDENCE_DEFAULTS[column])
    return out


def _append_tsarina_only_excluded_rows(df: pd.DataFrame, local: pd.DataFrame) -> pd.DataFrame:
    """Append local evidence rows that oncoref excludes from its CTA universe.

    H1-6 is the motivating row: oncoref removes histones before `cta_evidence()`,
    but tsarina still needs to expose the local RECURRENT_HEALTHY MS signal for
    auditability while keeping the gene out of canonical CTA defaults.
    """
    if "Ensembl_Gene_ID" not in df.columns or "Ensembl_Gene_ID" not in local.columns:
        return df

    present = set(_unversioned_gene_ids(df["Ensembl_Gene_ID"]))
    local_ids = _unversioned_gene_ids(local["Ensembl_Gene_ID"])
    rows: list[dict] = []
    for gene_id, specificity in _TSARINA_ONLY_EXCLUDED_CTA_ROWS.items():
        if gene_id in present:
            continue
        matches = local[local_ids == gene_id]
        if matches.empty:
            continue
        row = {column: matches.iloc[0].get(column, "") for column in df.columns}
        row.update(specificity)
        for column, default in _MS_EVIDENCE_DEFAULTS.items():
            row[column] = row.get(column, default) or default
        rows.append(row)

    if not rows:
        return df
    return pd.DataFrame([*df.to_dict("records"), *rows], columns=df.columns)


@lru_cache(maxsize=1)
def _load_cta_dataframe() -> pd.DataFrame:
    local = _local_cta_dataframe()
    df = _oncoref_cta_dataframe()
    df = _attach_local_ms_evidence(df, local)
    df = _append_tsarina_only_excluded_rows(df, local)

    # Recompute tsarina's synthesized axes after the MS evidence join. oncoref
    # owns the HPA tissue evidence and specificity decisions; tsarina keeps the
    # target-selection-layer MS restriction contribution.
    from .tiers import assign_all_axes

    df = assign_all_axes(df)
    # Surface the legacy inclusion flag alongside the canonical one so
    # downstream consumers that still schema-check for ``filtered`` (e.g.
    # pirlygenes pre-#241) can read the live tsarina table without
    # falling back to their own bundled snapshot.  See tsarina#61.
    if PASSES_FILTERS_COLUMN in df.columns and LEGACY_FILTERED_COLUMN not in df.columns:
        df[LEGACY_FILTERED_COLUMN] = df[PASSES_FILTERS_COLUMN]
    return df


def cta_dataframe() -> pd.DataFrame:
    """Return the full CTA evidence DataFrame (cached after first load)."""
    return _load_cta_dataframe()


def _boolish(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    return values.astype(str).str.lower().isin({"true", "1", "yes"})


def passes_filters_mask(df: pd.DataFrame) -> pd.Series:
    """Return a boolean mask for rows passing CTA curation filters.

    The bundled evidence table now uses ``passes_filters``. The legacy
    ``filtered`` column is accepted so tests and callers using older evidence
    tables still behave as before.
    """
    if PASSES_FILTERS_COLUMN in df.columns:
        values = df[PASSES_FILTERS_COLUMN]
    elif LEGACY_FILTERED_COLUMN in df.columns:
        values = df[LEGACY_FILTERED_COLUMN]
    else:
        return pd.Series(True, index=df.index)

    return _boolish(values)


def canonical_filtered_mask(df: pd.DataFrame) -> pd.Series:
    """Rows in the canonical CTA specificity tier.

    oncoref supplies explicit specificity decisions on top of the raw HPA
    ``passes_filters`` column. When those columns are present, use them for
    public CTA set helpers; otherwise fall back to the historical filter mask.
    """
    if {SPECIFICITY_ACTION_COLUMN, SPECIFICITY_STATUS_COLUMN} <= set(df.columns):
        action = df[SPECIFICITY_ACTION_COLUMN].astype(str)
        status = df[SPECIFICITY_STATUS_COLUMN].astype(str)
        return action.eq("include_default") | status.eq("canonical_low_expression")
    return passes_filters_mask(df)


def canonical_default_mask(df: pd.DataFrame) -> pd.Series:
    """Rows in the recommended expressed canonical CTA default set."""
    if SPECIFICITY_ACTION_COLUMN in df.columns:
        return df[SPECIFICITY_ACTION_COLUMN].astype(str).eq("include_default")

    mask = passes_filters_mask(df)
    if "never_expressed" in df.columns:
        from .tissues import MANUALLY_EXPRESSED_CTA

        never_expressed = _boolish(df["never_expressed"])
        if "Ensembl_Gene_ID" in df.columns:
            rescued = _unversioned_gene_ids(df["Ensembl_Gene_ID"]).isin(MANUALLY_EXPRESSED_CTA)
            never_expressed = never_expressed & ~rescued
        mask = mask & ~never_expressed
    return mask
