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

"""CLI smoke tests for ``tsarina panel`` and deprecated ``spanning`` alias.

Mirrors the pattern in ``test_cli_hits.py`` / ``test_cli_personalize.py``:
exercise the argparse surface via subprocess so flag parsing, choices
enforcement, and help-text sanity are pinned down independently of the
library-level coverage in ``test_spanning.py``.
"""

import subprocess
import sys

HEADLINE_FLAGS = (
    "--cta-count",
    "--ctas",
    "--panel",
    "--format",
)

PANEL_ONLY_FLAGS = (
    "--alleles",
    "--monoallelic-ms-max-percentile",
    "--sample-allele-ms-max-percentile",
    "--unrestricted-ms-max-percentile",
    "--include-predicted-only",
    "--predicted-only-max-percentile",
    "--peptides-per-cell",
    "--selection-allowlist",
    "--no-vital-tissue-filter",
    "--vital-tissue-max-ntpm",
    "--allow-non-magea4-mage-family",
    "--show-empty-ctas",
    "--no-summary",
    "--progress",
    "--predictor",
)

#: Deprecated flags hidden from --help but still parsed for back-compat.
DEPRECATED_HIDDEN_FLAGS = (
    "--max-percentile",
    "--no-progress",
    "--progress-bars",
    "--no-progress-bars",
)


def _run_cli(*args: str, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "tsarina.cli", *args],
        capture_output=True,
        text=True,
        check=check,
    )


def test_panel_help_exits_zero():
    r = _run_cli("panel", "--help")
    assert r.returncode == 0
    for flag in (*HEADLINE_FLAGS, *PANEL_ONLY_FLAGS):
        assert flag in r.stdout, f"missing help text for {flag}"


def test_panel_deprecated_flags_hidden_from_help():
    r = _run_cli("panel", "--help")
    for flag in DEPRECATED_HIDDEN_FLAGS:
        assert flag not in r.stdout, f"{flag} should be hidden from --help"


def test_panel_deprecated_flags_still_parse():
    """Hidden deprecated flags must keep working (back-compat aliases)."""
    import argparse

    from tsarina import cli_spanning

    p = argparse.ArgumentParser()
    cli_spanning._configure_parser(p)
    args = p.parse_args(
        [
            "--max-percentile",
            "1.0",
            "--no-progress",
            "--progress-bars",
            "--iedb-path",
            "x.csv",
            "--cedar-path",
            "y.csv",
            "--ctas",
            "MAGEA4",
        ]
    )
    assert args.max_percentile == 1.0
    assert args.progress == "off"
    assert args.progress_bars is True
    assert args.iedb_path == "x.csv"
    assert args.cedar_path == "y.csv"


def test_panel_iedb_alias_matches_iedb_path():
    import argparse

    from tsarina import cli_spanning

    p = argparse.ArgumentParser()
    cli_spanning._configure_parser(p)
    assert p.parse_args(["--iedb", "z.csv"]).iedb_path == "z.csv"


def test_panel_unknown_panel_rejected():
    r = _run_cli("panel", "--panel", "does-not-exist", check=False)
    assert r.returncode != 0
    combined = r.stderr + r.stdout
    # argparse reports the valid choices on invalid input
    assert "iedb27_ab" in combined


def test_panel_unknown_format_rejected():
    r = _run_cli("panel", "--format", "grid", check=False)
    assert r.returncode != 0
    combined = r.stderr + r.stdout
    assert "table" in combined and "wide" in combined and "long" in combined


def test_panel_unknown_predictor_rejected():
    r = _run_cli("panel", "--predictor", "bogus", check=False)
    assert r.returncode != 0
    combined = r.stderr + r.stdout
    assert "mhcflurry" in combined


def test_panel_is_registered_in_main_help():
    """``tsarina --help`` should list panel as the visible top-level subcommand."""
    r = _run_cli("--help")
    assert r.returncode == 0
    assert "panel" in r.stdout
    assert "Build a CTA x HLA pMHC matrix for a population HLA" in r.stdout
    assert "spanning" not in r.stdout


def test_spanning_alias_help_exits_zero():
    r = _run_cli("spanning", "--help")
    assert r.returncode == 0
    for flag in HEADLINE_FLAGS:
        assert flag in r.stdout, f"missing help text for {flag}"
    assert "--include-predicted-only" in r.stdout


def test_spanning_alias_unknown_panel_rejected():
    r = _run_cli("spanning", "--panel", "does-not-exist", check=False)
    assert r.returncode != 0
    combined = r.stderr + r.stdout
    assert "iedb27_ab" in combined


def test_spanning_alias_unknown_format_rejected():
    r = _run_cli("spanning", "--format", "grid", check=False)
    assert r.returncode != 0
    combined = r.stderr + r.stdout
    assert "table" in combined and "wide" in combined and "long" in combined
