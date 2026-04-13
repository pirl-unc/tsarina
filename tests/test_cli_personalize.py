import subprocess
import sys


def _run_cli(*args: str, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "tsarina.cli", *args],
        capture_output=True,
        text=True,
        check=check,
    )


def test_personalize_help_exits_zero():
    r = _run_cli("personalize", "--help")
    assert r.returncode == 0
    assert "personalize" in r.stdout.lower()
    assert "--hla" in r.stdout
    assert "--cta" in r.stdout
    assert "--predictor" in r.stdout


def test_personalize_skip_ms_evidence_no_score():
    r = _run_cli(
        "personalize",
        "--hla",
        "HLA-A*02:01",
        "--viruses",
        "",
        "--no-score",
        "--skip-ms-evidence",
    )
    assert r.returncode == 0, r.stderr


def test_personalize_requires_hla():
    r = _run_cli("personalize", check=False)
    assert r.returncode != 0
    assert "--hla" in r.stderr
