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

"""Load oncoref-owned CTA evidence with Tsarina's MS safety overlay."""

from __future__ import annotations

from functools import lru_cache
from os.path import dirname, join

import pandas as pd
from oncoref.cta import cta_evidence as oncoref_cta_evidence

_DATA_DIR = join(dirname(__file__), "data")
_MS_EVIDENCE_PATH = join(_DATA_DIR, "gene-ms-safety-evidence.csv")
PASSES_FILTERS_COLUMN = "passes_filters"
LEGACY_FILTERED_COLUMN = "filtered"

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
_MS_EVIDENCE_SCHEMA = ("Ensembl_Gene_ID", *_MS_EVIDENCE_COLUMNS)


def _unversioned_gene_ids(values: pd.Series) -> pd.Series:
    return values.astype(str).str.split(".").str[0]


def _local_ms_evidence() -> pd.DataFrame:
    """Return the gene-level MS overlay.

    This table deliberately contains no CTA symbols, membership flags, HPA
    columns, or specificity decisions. It is a generic gene-evidence table;
    whether a row is a CTA is determined exclusively by oncoref.
    """
    df = pd.read_csv(_MS_EVIDENCE_PATH, dtype=str, keep_default_na=False)
    missing = set(_MS_EVIDENCE_SCHEMA) - set(df.columns)
    if missing:
        raise ValueError(f"{_MS_EVIDENCE_PATH} is missing MS evidence column(s): {sorted(missing)}")
    df = df[list(_MS_EVIDENCE_SCHEMA)].copy()
    df["Ensembl_Gene_ID"] = _unversioned_gene_ids(df["Ensembl_Gene_ID"])
    duplicates = df.loc[df["Ensembl_Gene_ID"].duplicated(), "Ensembl_Gene_ID"]
    if not duplicates.empty:
        raise ValueError(
            f"gene-ms-safety-evidence.csv has duplicate Ensembl gene IDs: {sorted(set(duplicates))}"
        )
    return df


def _oncoref_cta_dataframe() -> pd.DataFrame:
    """Return oncoref's canonical CTA evidence frame."""
    return oncoref_cta_evidence().copy()


def _attach_local_ms_evidence(df: pd.DataFrame, ms: pd.DataFrame) -> pd.DataFrame:
    """Left-join Tsarina's gene-level MS safety evidence onto oncoref CTAs."""
    if "Ensembl_Gene_ID" not in df.columns:
        return df

    overlap = set(_MS_EVIDENCE_COLUMNS) & set(df.columns)
    if overlap:
        raise ValueError(
            "oncoref CTA evidence unexpectedly contains Tsarina-owned MS column(s): "
            f"{sorted(overlap)}"
        )

    join_key = "_tsarina_gene_id"
    out = df.copy()
    out[join_key] = _unversioned_gene_ids(out["Ensembl_Gene_ID"])
    overlay = ms.rename(columns={"Ensembl_Gene_ID": join_key})
    out = out.merge(overlay, on=join_key, how="left").drop(columns=[join_key])
    for column, default in _MS_EVIDENCE_DEFAULTS.items():
        out[column] = out[column].fillna(default)
    return out


@lru_cache(maxsize=1)
def _load_cta_dataframe() -> pd.DataFrame:
    df = _attach_local_ms_evidence(_oncoref_cta_dataframe(), _local_ms_evidence())

    # Compatibility only: the canonical raw-HPA flag is oncoref's
    # ``passes_filters``. Older consumers may still schema-check ``filtered``.
    if PASSES_FILTERS_COLUMN in df.columns and LEGACY_FILTERED_COLUMN not in df.columns:
        df[LEGACY_FILTERED_COLUMN] = df[PASSES_FILTERS_COLUMN]
    return df


def cta_dataframe() -> pd.DataFrame:
    """Return the oncoref CTA evidence frame with Tsarina MS annotations."""
    return _load_cta_dataframe()
