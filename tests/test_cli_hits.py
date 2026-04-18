import argparse
import subprocess
import sys
from unittest.mock import patch

import pandas as pd


def _run_cli(*args: str, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "tsarina.cli", *args],
        capture_output=True,
        text=True,
        check=check,
    )


def test_hits_help_exits_zero():
    r = _run_cli("hits", "--help")
    assert r.returncode == 0
    assert "--gene" in r.stdout
    assert "--allele" in r.stdout
    assert "--serotype" in r.stdout
    assert "--species" in r.stdout
    assert "--format" in r.stdout
    assert "--predict" in r.stdout
    assert "--mono-allelic-only" in r.stdout


def test_hits_help_describes_include_binding_assays_union():
    """tsarina#4 closed: cached path now routes through hitlist's
    load_all_evidence() union; help text should describe the UNION
    behavior and the evidence_kind tagging, not the old no-op limitation."""
    r = _run_cli("hits", "--help")
    assert r.returncode == 0
    assert "--include-binding-assays" in r.stdout
    assert "evidence_kind" in r.stdout


def test_hits_requires_gene_or_uniprot():
    r = _run_cli("hits", check=False)
    assert r.returncode != 0
    combined = r.stderr + r.stdout
    assert "--gene" in combined or "--uniprot" in combined


def test_hits_gene_and_uniprot_are_mutually_exclusive():
    r = _run_cli(
        "hits",
        "--gene",
        "PRAME",
        "--uniprot",
        "P08819",
        check=False,
    )
    assert r.returncode != 0


def test_include_binding_assays_routes_to_load_all_evidence(tmp_path):
    """tsarina#4 closed via hitlist 1.10.0: --include-binding-assays on the
    cached path must call load_all_evidence (union of MS + binding), not
    load_observations (MS-only).  This is the contract that makes the flag
    actually honor its name."""
    from tsarina import cli_hits

    args = argparse.Namespace(
        gene="PRAME",
        uniprot=None,
        allele=[],
        serotype=[],
        species="Homo sapiens",
        mhc_class="I",
        min_resolution=None,
        lengths=(8, 9, 10, 11),
        ensembl_release=112,
        include_binding_assays=True,
        mono_allelic_only=False,
        format="pmhc",
        predict=False,
        predictor="mhcflurry",
        iedb_path=None,
        cedar_path=None,
        skip_ms_evidence=False,
        output=str(tmp_path / "out.csv"),
    )
    empty = pd.DataFrame(
        {
            "peptide": pd.Series(dtype=str),
            "mhc_restriction": pd.Series(dtype=str),
            "evidence_kind": pd.Series(dtype=str),
        }
    )
    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_all_evidence", return_value=empty) as all_ev,
        patch("hitlist.observations.load_observations") as ms_only,
    ):
        cli_hits.handle(args)
    all_ev.assert_called_once()
    ms_only.assert_not_called()


def test_default_cached_path_uses_load_observations_not_union(tmp_path):
    """Regression guard: without --include-binding-assays, the cached path
    must still use load_observations (MS-only), not the union."""
    from tsarina import cli_hits

    args = argparse.Namespace(
        gene="PRAME",
        uniprot=None,
        allele=[],
        serotype=[],
        species="Homo sapiens",
        mhc_class="I",
        min_resolution=None,
        lengths=(8, 9, 10, 11),
        ensembl_release=112,
        include_binding_assays=False,
        mono_allelic_only=False,
        format="pmhc",
        predict=False,
        predictor="mhcflurry",
        iedb_path=None,
        cedar_path=None,
        skip_ms_evidence=False,
        output=str(tmp_path / "out.csv"),
    )
    empty = pd.DataFrame({"peptide": pd.Series(dtype=str), "mhc_restriction": pd.Series(dtype=str)})
    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_observations", return_value=empty) as ms_only,
        patch("hitlist.observations.load_all_evidence") as all_ev,
    ):
        cli_hits.handle(args)
    ms_only.assert_called_once()
    all_ev.assert_not_called()
