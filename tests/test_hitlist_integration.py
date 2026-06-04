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

"""Real, non-mocked hitlist integration smoke test (tsarina#82).

Every other test in the suite mocks hitlist, so CI cannot catch a break in
the *contract* between tsarina and the installed hitlist package — e.g. the
``observations`` / ``peptide_mappings`` schema splits in hitlist 1.30.46
(gene columns) and 1.30.57 (``cell_name`` -> ``cell_line_name`` +
``cell_type``).  The dev env once silently drifted to hitlist 1.30.0 (below
the pin) while the mocked suite stayed green.

This test runs the **real** ``tsarina hits`` cached path against a tiny,
committed, real-schema fixture index (``tests/fixtures/hitlist_mini/``,
sliced from a built hitlist index for a handful of MAGEA4 peptides).  No
mocking, no proteome download, no multi-minute build.  It exercises the
fragile coupling points:

- ``hitlist.observations.load_observations`` (gene pushdown)
- the gene-column re-attach: ``load_peptide_mappings`` (singular
  ``gene_name``/``gene_id``/``protein_id``) + ``annotate_observations_with_genes``
  (-> plural semicolon-joined ``gene_names``/``gene_ids``/``protein_ids``)
- ``hitlist.aggregate.aggregate_per_pmhc_with_refs`` output columns

If hitlist changes a schema or signature tsarina depends on, this fails at
PR time instead of at runtime.
"""

from __future__ import annotations

import io
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

_FIXTURE_DIR = Path(__file__).parent / "fixtures" / "hitlist_mini"


def _run_tsarina_hits(env_data_dir: Path) -> subprocess.CompletedProcess:
    import os

    env = {**os.environ, "HITLIST_DATA_DIR": str(env_data_dir)}
    return subprocess.run(
        [
            sys.executable,
            "-c",
            "from tsarina.cli import main; main()",
            "hits",
            "--gene",
            "MAGEA4",
            "--mhc-class",
            "I",
            "--format",
            "pmhc",
        ],
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
    )


@pytest.mark.skipif(
    not (_FIXTURE_DIR / "observations.parquet").exists(),
    reason="hitlist_mini fixture index missing",
)
def test_tsarina_hits_against_real_hitlist_fixture():
    """The real `tsarina hits` path consumes the real hitlist index schema."""
    result = _run_tsarina_hits(_FIXTURE_DIR)
    assert result.returncode == 0, (
        f"tsarina hits failed against real hitlist:\nSTDERR:\n{result.stderr}"
    )

    df = pd.read_csv(io.StringIO(result.stdout))

    # Gene-column re-attach contract: singular sidecar columns -> plural,
    # semicolon-joined columns that tsarina (and pirlygenes) consume.
    for col in ("peptide", "gene_names", "gene_ids", "protein_ids", "mhc_restriction"):
        assert col in df.columns, f"missing expected column {col!r}; got {list(df.columns)}"

    # aggregate_per_pmhc_with_refs output contract.
    for col in ("ms_pmhc_hit_count", "ms_pmhc_ref_count", "ms_pmhc_pmids"):
        assert col in df.columns, f"missing aggregation column {col!r}; got {list(df.columns)}"

    assert len(df) > 0, "expected at least one MAGEA4 pMHC row"
    assert (df["gene_names"].astype(str).str.contains("MAGEA4")).any()

    # A known MAGEA4 peptide from the fixture, with PMID-backed MS support.
    gvy = df[df["peptide"] == "GVYDGREHTV"]
    assert len(gvy) > 0, "expected the GVYDGREHTV MAGEA4 peptide in the fixture"
    assert gvy["ms_pmhc_pmids"].astype(str).str.len().gt(0).any()
