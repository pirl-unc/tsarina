from pathlib import Path

import pytest

from tsarina.datasources import (
    DatasetNotRegisteredError,
    resolve_cedar_path,
    resolve_dataset_paths,
    resolve_iedb_path,
)


def test_explicit_iedb_path_passes_through(tmp_path: Path):
    p = tmp_path / "iedb.csv"
    p.write_text("header\n")
    assert resolve_iedb_path(p) == p


def test_explicit_cedar_path_passes_through(tmp_path: Path):
    p = tmp_path / "cedar.csv"
    p.write_text("header\n")
    assert resolve_cedar_path(p) == p


def test_missing_iedb_raises_with_register_hint(monkeypatch):
    from tsarina import datasources as ds

    def _raise(_name):
        raise KeyError("iedb")

    monkeypatch.setattr(ds, "get_path", _raise)
    with pytest.raises(DatasetNotRegisteredError) as excinfo:
        resolve_iedb_path(None)
    assert "tsarina data register iedb" in str(excinfo.value)


def test_missing_cedar_silently_returns_none(monkeypatch):
    from tsarina import datasources as ds

    def _raise(_name):
        raise KeyError("cedar")

    monkeypatch.setattr(ds, "get_path", _raise)
    assert resolve_cedar_path(None) is None


def test_resolve_dataset_paths_require_iedb_false(monkeypatch):
    from tsarina import datasources as ds

    def _raise(_name):
        raise KeyError("missing")

    monkeypatch.setattr(ds, "get_path", _raise)
    iedb, cedar = resolve_dataset_paths(None, None, require_iedb=False)
    assert iedb is None
    assert cedar is None
