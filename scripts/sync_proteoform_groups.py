#!/usr/bin/env python
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

"""Sync the proteoform-group registry from ``cancerdata`` into tsarina.

The identical-protein CGA grouping (CTAG1A/CTAG1B, the CT47A family, …) used by
the pMHC-collapse logic in ``spanning.py`` is *owned* by ``cancerdata``
(``cancerdata.proteoform_groups()``, derived from Ensembl protein sequences).
tsarina ships a materialized mirror so the grouping is deterministic and needs no
runtime dependency on cancerdata; this script regenerates that mirror, and
``tests/test_spanning.py`` asserts the two stay in sync when cancerdata is
importable.

Run (with cancerdata installed/importable):
    python scripts/sync_proteoform_groups.py
"""

from __future__ import annotations

from pathlib import Path

_OUTPUT = Path(__file__).resolve().parent.parent / "tsarina" / "data" / "proteoform-groups.csv"
_COLUMNS = ["proteoform_id", "member_symbol", "member_gene_id", "protein_length", "n_members"]


def main() -> None:
    from cancerdata import proteoform_groups

    df = proteoform_groups()[_COLUMNS].sort_values(["proteoform_id", "member_symbol"])
    _OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(_OUTPUT, index=False)
    n_groups = df["proteoform_id"].nunique()
    print(f"Wrote {len(df)} member rows across {n_groups} proteoform groups to {_OUTPUT}")


if __name__ == "__main__":
    main()
