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

"""Tests for datacache.cli_help integration, argcomplete, --version, and --cache-root."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import hitlist.downloads as hd
import pytest
from datacache.cli_help import ANSI_BOLD_CYAN, ANSI_CYAN, HelpArgumentParser

import tsarina.cli as cli
from tsarina.version import __version__


@pytest.fixture(autouse=True)
def _reset_hitlist_cache_state():
    hd.reset_data_dir_cache()
    hd.set_offline(None)
    yield
    hd.reset_data_dir_cache()
    hd.set_offline(None)


def _run_main(monkeypatch, capsys, argv: list[str]) -> tuple[int, str, str]:
    monkeypatch.setattr(sys, "argv", ["tsarina", *argv])
    rc = 0
    try:
        cli.main()
    except SystemExit as exc:
        rc = int(exc.code) if exc.code is not None else 0
    captured = capsys.readouterr()
    return rc, captured.out, captured.err


def test_cli_has_python_argcomplete_ok_marker():
    cli_path = Path(cli.__file__).resolve()
    header = cli_path.read_text(encoding="utf-8")[:1024]
    assert "PYTHON_ARGCOMPLETE_OK" in header


def test_parser_uses_datacache_help_argument_parser():
    parser = cli.build_parser()
    assert isinstance(parser, HelpArgumentParser)
    # format_help() must stay plain even before print_help() applies TTY styling.
    assert "\033[" not in parser.format_help()


def test_help_tty_styling_and_no_color_and_dumb_term(monkeypatch, capsys):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")
    monkeypatch.setattr(sys.stdout, "isatty", lambda: True)

    rc, out, err = _run_main(monkeypatch, capsys, ["--help"])
    assert rc == 0 and not err
    assert ANSI_BOLD_CYAN in out
    assert ANSI_CYAN in out

    monkeypatch.setenv("NO_COLOR", "1")
    rc2, out2, err2 = _run_main(monkeypatch, capsys, ["--help"])
    assert rc2 == 0 and not err2
    assert "\033[" not in out2

    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "dumb")
    rc3, out3, err3 = _run_main(monkeypatch, capsys, ["--help"])
    assert rc3 == 0 and not err3
    assert "\033[" not in out3


def test_main_invokes_enable_argcomplete(monkeypatch, capsys):
    seen: list[object] = []

    def fake_enable(parser, **kwargs):
        seen.append(parser)
        return True

    monkeypatch.setattr(cli, "enable_argcomplete", fake_enable)
    rc, out, err = _run_main(monkeypatch, capsys, ["--version"])
    assert rc == 0 and not err
    assert f"tsarina {__version__}" in out
    assert len(seen) == 1
    assert isinstance(seen[0], HelpArgumentParser)


def test_data_dirs_cache_root_overrides_env_without_leaking(monkeypatch, tmp_path: Path, capsys):
    ambient_hitlist = tmp_path / "ambient-hitlist"
    ambient_openvax = tmp_path / "ambient-openvax"
    explicit_root = tmp_path / "explicit-cache"

    monkeypatch.setenv("HITLIST_DATA_DIR", str(ambient_hitlist))
    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(ambient_openvax))
    hd.reset_data_dir_cache()

    for argv in (
        ["data", "dirs", "--cache-root", str(explicit_root), "--json"],
        ["data", "--cache-root", str(explicit_root), "dirs", "--json"],
    ):
        rc, out, err = _run_main(monkeypatch, capsys, argv)
        assert rc == 0 and not err
        payload = json.loads(out)
        assert payload["data_dir"] == str(explicit_root / "hitlist")
        assert payload["origin"] == "openvax_env"
        assert payload["exists"] is False
        assert not explicit_root.exists(), "data dirs --cache-root must remain read-only"

        # Environment and hitlist cache resolution must be fully restored.
        assert os.environ.get("HITLIST_DATA_DIR") == str(ambient_hitlist)
        assert os.environ.get("OPENVAX_DATA_CACHE") == str(ambient_openvax)
        resolved_after, origin_after = hd.resolve_data_dir()
        assert resolved_after == ambient_hitlist
        assert origin_after == "env"


def test_data_register_list_path_verify_prune_honor_cache_root(monkeypatch, tmp_path: Path, capsys):
    ambient_hitlist = tmp_path / "ambient-hitlist"
    explicit_root = tmp_path / "scoped-cache"
    monkeypatch.setenv("HITLIST_DATA_DIR", str(ambient_hitlist))
    hd.reset_data_dir_cache()

    sample = tmp_path / "hpv16.fasta"
    sample.write_text(">p1\nMAAAAA\n", encoding="utf-8")

    rc_reg, _, err_reg = _run_main(
        monkeypatch,
        capsys,
        ["data", "register", "hpv16", str(sample), "--cache-root", str(explicit_root)],
    )
    assert rc_reg == 0 and not err_reg
    assert not ambient_hitlist.exists()

    rc_path, out_path, err_path = _run_main(
        monkeypatch,
        capsys,
        ["data", "path", "hpv16", "--cache-root", str(explicit_root), "--json"],
    )
    assert rc_path == 0 and not err_path
    assert Path(json.loads(out_path)["path"]) == sample.resolve()

    rc_list, out_list, err_list = _run_main(
        monkeypatch,
        capsys,
        ["data", "list", "--cache-root", str(explicit_root), "--json"],
    )
    assert rc_list == 0 and not err_list
    assert [row["name"] for row in json.loads(out_list)] == ["hpv16"]

    rc_verify, out_verify, err_verify = _run_main(
        monkeypatch,
        capsys,
        ["data", "verify", "--cache-root", str(explicit_root), "--json"],
    )
    assert rc_verify == 0 and not err_verify
    assert json.loads(out_verify)["ok"] is True

    rc_prune, out_prune, err_prune = _run_main(
        monkeypatch,
        capsys,
        ["data", "prune", "--cache-root", str(explicit_root), "--json"],
    )
    assert rc_prune == 0 and not err_prune
    prune_payload = json.loads(out_prune)
    assert prune_payload["dry_run"] is True
    assert prune_payload["root"] == str(explicit_root / "hitlist")

    # Ambient data dir still has nothing registered.
    rc_ambient, out_ambient, err_ambient = _run_main(
        monkeypatch, capsys, ["data", "list", "--json"]
    )
    assert rc_ambient == 0 and not err_ambient
    assert json.loads(out_ambient) == []
    assert os.environ.get("HITLIST_DATA_DIR") == str(ambient_hitlist)
    assert "OPENVAX_DATA_CACHE" not in os.environ
