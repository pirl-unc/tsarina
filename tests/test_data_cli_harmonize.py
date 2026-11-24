"""Tests for harmonized tsarina data CLI (dirs, verify, prune, --offline, path --json)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import hitlist.downloads as hd
import pytest

from tsarina import cli


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


def test_data_dirs_is_read_only_and_supports_json(monkeypatch, tmp_path: Path, capsys):
    shared = tmp_path / "shared"
    monkeypatch.delenv("HITLIST_DATA_DIR", raising=False)
    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(shared))
    monkeypatch.delenv("OPENVAX_OFFLINE", raising=False)
    monkeypatch.delenv("DATACACHE_OFFLINE", raising=False)
    hd.reset_data_dir_cache()

    rc, out, err = _run_main(monkeypatch, capsys, ["data", "dirs", "--json"])
    assert rc == 0 and not err
    assert not shared.exists(), "data dirs must remain strictly read-only"
    payload = json.loads(out)
    assert payload["package"] == "tsarina"
    assert payload["data_dir"] == str(shared / "hitlist")
    assert payload["exists"] is False
    assert payload["offline"] is False

    rc2, out2, err2 = _run_main(monkeypatch, capsys, ["data", "dirs"])
    assert rc2 == 0 and not err2
    assert "data_dir" in out2
    assert "not created" in out2
    assert not shared.exists()


def test_data_verify_and_prune_dry_run_json(monkeypatch, tmp_path: Path, capsys):
    data_root = tmp_path / "hitlist-data"
    monkeypatch.setenv("HITLIST_DATA_DIR", str(data_root))
    hd.reset_data_dir_cache()

    rc, out, err = _run_main(monkeypatch, capsys, ["data", "verify", "--json"])
    assert rc == 0 and not err
    verify = json.loads(out)
    assert verify["ok"] is True
    assert verify["n_failed"] == 0

    rc2, out2, err2 = _run_main(monkeypatch, capsys, ["data", "prune", "--json"])
    assert rc2 == 0 and not err2
    prune = json.loads(out2)
    assert prune["dry_run"] is True
    assert prune["n_prunable_bundles"] == 0
    assert not data_root.exists()


def test_data_path_json(monkeypatch, tmp_path: Path, capsys):
    data_root = tmp_path / "hitlist-data"
    monkeypatch.setenv("HITLIST_DATA_DIR", str(data_root))
    hd.reset_data_dir_cache()
    sample = tmp_path / "hpv16.fasta"
    sample.write_text(">p1\nMAAAAA\n", encoding="utf-8")

    rc_reg, _, err_reg = _run_main(monkeypatch, capsys, ["data", "register", "hpv16", str(sample)])
    assert rc_reg == 0 and not err_reg

    rc, out, err = _run_main(monkeypatch, capsys, ["data", "path", "hpv16", "--json"])
    assert rc == 0 and not err
    payload = json.loads(out)
    assert payload["name"] == "hpv16"
    assert Path(payload["path"]) == sample.resolve()


def test_data_fetch_offline_blocks_network(monkeypatch, tmp_path: Path, capsys):
    data_root = tmp_path / "offline-hitlist"
    monkeypatch.setenv("HITLIST_DATA_DIR", str(data_root))
    hd.reset_data_dir_cache()

    rc, out, err = _run_main(monkeypatch, capsys, ["data", "fetch", "--offline", "hpv16"])
    assert rc == 1
    assert not out
    assert "offline" in err.lower()
    assert not hd.is_offline(), "--offline must not leak process state"
