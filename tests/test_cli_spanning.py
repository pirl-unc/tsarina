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

"""CLI smoke tests for ``tsarina spanning``.

Mirrors the pattern in ``test_cli_hits.py`` / ``test_cli_personalize.py``:
exercise the argparse surface via subprocess so flag parsing, choices
enforcement, and help-text sanity are pinned down independently of the
library-level coverage in ``test_spanning.py``.
"""

from __future__ import annotations

import subprocess
import sys


def _run_cli(*args: str, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "tsarina.cli", *args],
        capture_output=True,
        text=True,
        check=check,
    )


def test_spanning_help_exits_zero():
    r = _run_cli("spanning", "--help")
    assert r.returncode == 0
    for flag in (
        "--cta-count",
        "--ctas",
        "--panel",
        "--alleles",
        "--max-percentile",
        "--format",
        "--predictor",
    ):
        assert flag in r.stdout, f"missing help text for {flag}"


def test_spanning_unknown_panel_rejected():
    r = _run_cli("spanning", "--panel", "does-not-exist", check=False)
    assert r.returncode != 0
    combined = r.stderr + r.stdout
    # argparse reports the valid choices on invalid input
    assert "iedb27_ab" in combined


def test_spanning_unknown_format_rejected():
    r = _run_cli("spanning", "--format", "grid", check=False)
    assert r.returncode != 0
    combined = r.stderr + r.stdout
    assert "wide" in combined and "long" in combined


def test_spanning_unknown_predictor_rejected():
    r = _run_cli("spanning", "--predictor", "bogus", check=False)
    assert r.returncode != 0
    combined = r.stderr + r.stdout
    assert "mhcflurry" in combined


def test_spanning_is_registered_in_main_help():
    """``tsarina --help`` should list spanning as a top-level subcommand."""
    r = _run_cli("--help")
    assert r.returncode == 0
    assert "spanning" in r.stdout
