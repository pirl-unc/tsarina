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

"""Guard: no rendered CLI line should begin with '(' — a parenthetical aside at
the start of a line reads awkwardly.

Exercises the runnable output surface (every command's --help, which renders all
argparse help/description/epilog, plus the list/reference/error paths) and
asserts no emitted line starts with an open paren after indentation. Pinned to
COLUMNS=80 because argparse wraps help to the terminal width, and a parenthetical
can land at the start of a *wrapped* continuation line at narrow widths.

Run in-process (not via subprocess) so the heavy CLI imports happen once.
"""

import contextlib
import sys

import pytest

import tsarina.cli as cli

# Commands that run without a built index or network access.
_COMMANDS = [
    ["--help"],
    ["data", "--help"],
    ["data", "list"],
    ["data", "available"],
    ["build", "--help"],
    ["build", "observations", "--help"],
    ["reference", "--help"],
    ["reference", "list"],
    ["reference", "path", "hpa_rna_consensus"],
    ["hits", "--help"],
    ["personalize", "--help"],
    ["panel", "--help"],
    ["data", "path", "definitely-not-a-dataset"],  # error path (stderr)
    ["data", "remove", "definitely-not-a-dataset"],  # warn path (stderr)
]


def _offending(text: str) -> list[str]:
    return [ln for ln in text.splitlines() if ln.lstrip().startswith("(")]


@pytest.mark.parametrize("argv", _COMMANDS, ids=lambda a: " ".join(a))
def test_no_cli_line_starts_with_paren(argv, monkeypatch, capsys, tmp_path):
    monkeypatch.setenv("TSARINA_DATA_DIR", str(tmp_path))
    monkeypatch.setenv("COLUMNS", "80")  # deterministic help wrapping
    monkeypatch.setattr(sys, "argv", ["tsarina", *argv])
    # --help and error paths exit; we only care about the emitted text.
    with contextlib.suppress(SystemExit):
        cli.main()
    captured = capsys.readouterr()
    offenders = _offending(captured.out) + _offending(captured.err)
    assert not offenders, f"`tsarina {' '.join(argv)}` emitted leading-paren line(s): {offenders}"
