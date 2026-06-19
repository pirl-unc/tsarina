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

"""tsarina reuses hitlist's shared subtle help coloring (the formatter's own
alignment/gating invariants are tested in hitlist.tests.test_cli_help)."""

import contextlib
import sys

import tsarina.cli as cli


def _help(monkeypatch, capsys, **env):
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    monkeypatch.setattr(sys, "argv", ["tsarina", "--help"])
    with contextlib.suppress(SystemExit):
        cli.main()
    return capsys.readouterr().out


def test_help_colored_when_forced(monkeypatch, capsys):
    monkeypatch.delenv("NO_COLOR", raising=False)
    out = _help(monkeypatch, capsys, FORCE_COLOR="1")
    assert "\033[36m" in out  # subcommand names in cyan
    assert "\033[1m" in out  # bold section headings


def test_help_plain_when_no_color(monkeypatch, capsys):
    monkeypatch.delenv("FORCE_COLOR", raising=False)
    out = _help(monkeypatch, capsys, NO_COLOR="1")
    assert "\033[" not in out
