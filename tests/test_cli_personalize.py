import argparse
import subprocess
import sys

import pytest

from tsarina.cli_personalize import _parse_cta, _parse_hla


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


def test_personalize_hla_accepts_unquoted_space_separated_no_asterisk():
    """HLA-A0201 (no '*') and space-separated, unquoted alleles both work --
    the '*' is a shell glob character, so a form without it never needs
    quoting."""
    r = _run_cli(
        "personalize",
        "--hla",
        "HLA-A0201",
        "HLA-B0702",
        "--no-score",
        "--skip-ms-evidence",
    )
    assert r.returncode == 0, r.stderr


def test_personalize_cta_accepts_unquoted_space_separated():
    """Multi-token --cta parses correctly regardless of what happens next
    in the pipeline -- pin this to argument parsing, not real CTA peptide
    generation (which needs a real pyensembl reference download and
    shouldn't gate a CLI-parsing test). A TPM below --min-cta-tpm's
    default floor means personalized_targets() never reaches that step."""
    r = _run_cli(
        "personalize",
        "--hla",
        "HLA-A*02:01",
        "--cta",
        "MAGEA4=0.1",
        "PRAME=0.1",
        "--no-score",
        "--skip-ms-evidence",
        check=False,
    )
    assert r.returncode == 0, r.stderr


def test_personalize_hla_rejects_unrecognized_allele():
    r = _run_cli(
        "personalize",
        "--hla",
        "NOT-A-REAL-ALLELE",
        "--no-score",
        "--skip-ms-evidence",
        check=False,
    )
    assert r.returncode != 0
    assert "not a recognized allele" in r.stderr


def test_personalize_cta_rejects_malformed_entry():
    r = _run_cli(
        "personalize",
        "--hla",
        "HLA-A*02:01",
        "--cta",
        "NOT_GENE_EQUALS_TPM",
        "--no-score",
        "--skip-ms-evidence",
        check=False,
    )
    assert r.returncode != 0
    assert "GENE=TPM" in r.stderr


# ── _parse_hla / _parse_cta: flexible comma/space/quote handling ───────


def test_parse_hla_accepts_comma_separated_single_token():
    assert _parse_hla(["HLA-A*02:01,HLA-B*07:02"]) == ["HLA-A*02:01", "HLA-B*07:02"]


def test_parse_hla_accepts_space_separated_multi_token():
    assert _parse_hla(["HLA-A*02:01", "HLA-B*07:02"]) == ["HLA-A*02:01", "HLA-B*07:02"]


def test_parse_hla_normalizes_asterisk_free_and_bare_forms():
    assert _parse_hla(["HLA-A0201", "A0201", "hla-a*02:01"]) == [
        "HLA-A*02:01",
        "HLA-A*02:01",
        "HLA-A*02:01",
    ]


def test_parse_hla_rejects_unrecognized_token():
    with pytest.raises(argparse.ArgumentTypeError, match="not a recognized allele"):
        _parse_hla(["not-an-allele"])


def test_parse_cta_accepts_comma_and_space_mixed():
    assert _parse_cta(["MAGEA4=142.5,PRAME=87.3", "SSX1=10"]) == {
        "MAGEA4": 142.5,
        "PRAME": 87.3,
        "SSX1": 10.0,
    }


def test_parse_cta_rejects_missing_equals():
    with pytest.raises(argparse.ArgumentTypeError, match="GENE=TPM"):
        _parse_cta(["MAGEA4"])


def test_parse_cta_rejects_non_numeric_tpm():
    with pytest.raises(argparse.ArgumentTypeError, match="not a number"):
        _parse_cta(["MAGEA4=abc"])
