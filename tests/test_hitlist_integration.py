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
import shutil
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
def test_tsarina_hits_against_real_hitlist_fixture(tmp_path: Path):
    """The real `tsarina hits` path consumes the real hitlist index schema."""
    # Run against an isolated copy so nothing a build path writes alongside the
    # artifacts can land in the committed fixture directory.
    for fixture in _FIXTURE_DIR.iterdir():
        if fixture.is_file():
            shutil.copy2(fixture, tmp_path / fixture.name)

    result = _run_tsarina_hits(tmp_path)
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


@pytest.mark.skipif(
    not (_FIXTURE_DIR / "observations.parquet").exists(),
    reason="hitlist_mini fixture index missing",
)
def test_fixture_carries_the_columns_a_current_hitlist_build_writes():
    """The fixture is the schema tripwire, so it must not lag the builder.

    Regenerate with ``python scripts/regenerate_hitlist_mini_fixture.py``; the
    same script's ``--check`` mode reports drift against a locally built index.
    """
    observations = pd.read_parquet(_FIXTURE_DIR / "observations.parquet")
    for col in (
        "mhc_restriction",
        "allele_resolution",
        "serotype",
        "serotypes",
        "mhc_allele_set",
        "mhc_allele_provenance",
        # hitlist#415: how a row establishes its named restriction.
        "restriction_evidence",
        "is_monoallelic",
        "is_binding_assay",
    ):
        assert col in observations.columns, (
            f"fixture predates the {col!r} column; regenerate it with "
            "scripts/regenerate_hitlist_mini_fixture.py"
        )

    mappings = pd.read_parquet(_FIXTURE_DIR / "peptide_mappings.parquet")
    for col in ("peptide", "gene_name", "gene_id", "protein_id", "gene_biotype"):
        assert col in mappings.columns, (
            f"fixture predates the {col!r} mapping column; regenerate it with "
            "scripts/regenerate_hitlist_mini_fixture.py"
        )


@pytest.mark.skipif(
    not (_FIXTURE_DIR / "observations.parquet").exists(),
    reason="hitlist_mini fixture index missing",
)
def test_fixture_annotations_agree_with_the_installed_hitlist():
    """Stored MHC identity must still be something the installed hitlist writes.

    The fixture is real data annotated by a real build, so it has to stay
    consistent with hitlist -- but not identical to it, and the difference
    matters because ``develop.sh`` installs sibling checkouts. The installed
    hitlist is routinely ahead of the version floor, so demanding equality
    would fail every time upstream improved.

    Scalar classifications are canonical: one restriction resolves to one
    ``allele_resolution`` and one ``mhc_species``, so a change there is real
    drift and fails here (that is how hitlist 1.55.7 reclassifying promoted
    donor sets would surface).

    Serotype membership only grows. mhcgnomes gains specificities -- ``Cw16``
    was curated in from WHO's ``hla_nom.txt`` after this slice was cut, so a
    donor set containing C*16:01 now resolves one serotype more than the
    committed rows carry. That is the library improving, not the fixture
    rotting. What must never happen is the fixture asserting a serotype the
    installed hitlist no longer assigns at all.

    Regenerate with ``python scripts/regenerate_hitlist_mini_fixture.py`` when
    the ``hitlist`` floor moves; the assertions below hold either way.
    """
    from hitlist.curation import resolve_mhc_annotation

    def _tokens(value: object) -> set[str]:
        return {token for token in str(value).split(";") if token}

    observations = pd.read_parquet(_FIXTURE_DIR / "observations.parquet")
    for _, row in observations.iterrows():
        restriction = str(row["mhc_restriction"])
        expected = resolve_mhc_annotation(restriction).as_record_fields()

        for col in ("allele_resolution", "mhc_species"):
            assert row[col] == expected[col], (
                f"{restriction!r}: stored {col}={row[col]!r} but the installed "
                f"hitlist resolves {expected[col]!r}"
            )

        stored = _tokens(row["serotypes"])
        resolved = _tokens(expected["serotypes"])
        retired = stored - resolved
        assert not retired, (
            f"{restriction!r}: stored serotypes {sorted(retired)} are no longer "
            f"assigned by the installed hitlist, which resolves {sorted(resolved)}. "
            "Regenerate with scripts/regenerate_hitlist_mini_fixture.py."
        )

        canonical = str(row["serotype"])
        assert not canonical or canonical in resolved, (
            f"{restriction!r}: stored canonical serotype {canonical!r} is not among "
            f"{sorted(resolved)} for the installed hitlist."
        )
