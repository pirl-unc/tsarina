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


@lru_cache(maxsize=1)
def _load_cta_dataframe() -> pd.DataFrame:
    path = join(_DATA_DIR, "cancer-testis-antigens.csv")
    return pd.read_csv(path)


def cta_dataframe() -> pd.DataFrame:
    """Return the full CTA evidence DataFrame (cached after first load)."""
    return _load_cta_dataframe()
