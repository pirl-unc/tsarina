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


def test_hits_help_warns_about_include_binding_assays_cached_noop():
    """tsarina#4 part 1: the flag is a no-op on the cached path (hitlist
    purges binding rows at build time).  Help text must say so explicitly
    — otherwise users silently get MS-only output when they asked for both."""
    r = _run_cli("hits", "--help")
    assert r.returncode == 0
    assert "--include-binding-assays" in r.stdout
    # Must mention both the limitation and the escape hatch.
    assert "no-op" in r.stdout
    assert "--iedb" in r.stdout


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


def test_include_binding_assays_on_cached_path_warns(capsys, tmp_path):
    """Calling handle() with --include-binding-assays on the default cached
    path must emit a stderr warning, not silently succeed.  Tracks tsarina#4
    until hitlist#47 ships a separate binding-assay index."""
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
            "is_binding_assay": pd.Series(dtype=bool),
        }
    )
    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_observations", return_value=empty),
    ):
        cli_hits.handle(args)
    err = capsys.readouterr().err
    assert "--include-binding-assays" in err
    assert "cached" in err
    assert "--iedb" in err
