import argparse
import subprocess
import sys

import pytest

from tsarina.cli_personalize import _parse_cta, _parse_hla, _resolve_format


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


# A bare "--cta GENE" (no =TPM) is deliberately NOT covered by a CLI
# subprocess test: by design it bypasses the --min-cta-tpm floor and always
# runs real CTA peptide generation, which needs a warm pyensembl reference
# cache and is not something a CLI-parsing test should depend on. The
# parsing itself is covered by test_parse_cta_bare_gene_name_maps_to_nan_tpm,
# and the inclusion behavior by test_nan_tpm_included_even_below_min_cta_tpm
# in test_personalize.py.


def test_personalize_cta_rejects_non_numeric_tpm_value():
    r = _run_cli(
        "personalize",
        "--hla",
        "HLA-A*02:01",
        "--cta",
        "MAGEA4=not-a-number",
        "--no-score",
        "--skip-ms-evidence",
        check=False,
    )
    assert r.returncode != 0
    assert "not a number" in r.stderr


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


def test_parse_cta_bare_gene_name_maps_to_nan_tpm():
    import math

    out = _parse_cta(["MAGEA4"])
    assert set(out) == {"MAGEA4"}
    assert math.isnan(out["MAGEA4"])


def test_parse_cta_rejects_non_numeric_tpm():
    with pytest.raises(argparse.ArgumentTypeError, match="not a number"):
        _parse_cta(["MAGEA4=abc"])


# ── --format / --quiet ───────────────────────────────────────────────────


def test_personalize_default_format_is_table_on_stdout():
    r = _run_cli(
        "personalize", "--hla", "HLA-A*02:01", "--viruses", "", "--no-score", "--skip-ms-evidence"
    )
    assert r.returncode == 0, r.stderr
    assert "(no targets)" in r.stdout


def test_personalize_format_csv_explicit_on_stdout():
    r = _run_cli(
        "personalize",
        "--hla",
        "HLA-A*02:01",
        "--viruses",
        "",
        "--no-score",
        "--skip-ms-evidence",
        "--format",
        "csv",
    )
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip().startswith("peptide,length,category")


def test_personalize_quiet_flag_is_accepted():
    """CLI wiring smoke test -- real progress-message suppression is
    covered at the unit level (test_show_progress_false_is_silent in
    test_personalize.py), since this trivial no-candidates input never
    reaches a progress-reporting stage regardless of --quiet."""
    r = _run_cli(
        "personalize",
        "--hla",
        "HLA-A*02:01",
        "--viruses",
        "",
        "--no-score",
        "--skip-ms-evidence",
        "--quiet",
    )
    assert r.returncode == 0, r.stderr


# ── _resolve_format: --output extension sniffing ─────────────────────────


@pytest.mark.parametrize(
    "output, expected",
    [
        ("out.csv", "csv"),
        ("OUT.CSV", "csv"),
        ("out.tsv", "tsv"),
        ("out.tab", "tsv"),
        ("out.txt", "table"),
        ("/tmp/nested/path/report.txt", "table"),
        ("out.dat", "csv"),
        ("out", "csv"),
    ],
)
def test_resolve_format_sniffs_output_extension(output, expected):
    assert _resolve_format(None, output) == expected


def test_resolve_format_defaults_to_table_without_output():
    assert _resolve_format(None, None) == "table"


def test_resolve_format_explicit_flag_beats_extension():
    assert _resolve_format("table", "out.csv") == "table"
    assert _resolve_format("csv", None) == "csv"


def test_personalize_output_txt_writes_a_table(tmp_path):
    out = tmp_path / "targets.txt"
    r = _run_cli(
        "personalize",
        "--hla",
        "HLA-A*02:01",
        "--viruses",
        "",
        "--no-score",
        "--skip-ms-evidence",
        "--output",
        str(out),
    )
    assert r.returncode == 0, r.stderr
    assert "(no targets)" in out.read_text()


def test_personalize_output_csv_writes_csv(tmp_path):
    out = tmp_path / "targets.csv"
    r = _run_cli(
        "personalize",
        "--hla",
        "HLA-A*02:01",
        "--viruses",
        "",
        "--no-score",
        "--skip-ms-evidence",
        "--output",
        str(out),
    )
    assert r.returncode == 0, r.stderr
    assert out.read_text().startswith("peptide,length,category")


def test_personalize_output_tsv_writes_tabs(tmp_path):
    out = tmp_path / "targets.tsv"
    r = _run_cli(
        "personalize",
        "--hla",
        "HLA-A*02:01",
        "--viruses",
        "",
        "--no-score",
        "--skip-ms-evidence",
        "--output",
        str(out),
    )
    assert r.returncode == 0, r.stderr
    assert out.read_text().startswith("peptide\tlength\tcategory")
