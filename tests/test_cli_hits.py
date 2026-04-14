import subprocess
import sys


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
