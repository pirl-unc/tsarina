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

"""tsarina reuses datacache's shared subtle cyan help coloring."""

import contextlib
import sys

import tsarina.cli as cli


def _help(monkeypatch, capsys, *, tty: bool, **env):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    monkeypatch.setattr(sys.stdout, "isatty", lambda: tty)
    monkeypatch.setattr(sys, "argv", ["tsarina", "--help"])
    with contextlib.suppress(SystemExit):
        cli.main()
    return capsys.readouterr().out


def test_help_colored_on_tty(monkeypatch, capsys):
    out = _help(monkeypatch, capsys, tty=True)
    assert "\033[1;36m" in out  # bold cyan section headings
    assert "\033[36m" in out  # cyan option flags
    assert "\033[1musage:\033[0m" in out


def test_help_plain_when_no_color(monkeypatch, capsys):
    out = _help(monkeypatch, capsys, tty=True, NO_COLOR="1")
    assert "\033[" not in out


def test_help_plain_when_not_tty(monkeypatch, capsys):
    out = _help(monkeypatch, capsys, tty=False)
    assert "\033[" not in out
