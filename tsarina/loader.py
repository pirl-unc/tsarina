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


@lru_cache(maxsize=1)
def _load_cta_dataframe() -> pd.DataFrame:
    path = join(_DATA_DIR, "cancer-testis-antigens.csv")
    return pd.read_csv(path)


def cta_dataframe() -> pd.DataFrame:
    """Return the full CTA evidence DataFrame (cached after first load)."""
    return _load_cta_dataframe()


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

    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    return values.astype(str).str.lower().isin({"true", "1", "yes"})
