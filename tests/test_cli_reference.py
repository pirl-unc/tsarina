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

"""CLI smoke tests for ``tsarina reference`` plus the bare-``data`` default.

The ``reference`` checks run end-to-end via subprocess (no network: ``list`` /
``path`` only read the manifest / compute paths), isolated to a tmp cache via
``TSARINA_DATA_DIR`` (the env var ``reference_data`` honors). The bare-``data``
routing check is a hermetic unit test that stubs the hitlist-delegated
``list_datasets`` rather than touching hitlist's own data dir.
"""

import argparse
import os
import subprocess
import sys

import tsarina.cli as cli


def _run(*args: str, cache_dir, check: bool = True) -> subprocess.CompletedProcess:
    env = {**os.environ, "TSARINA_DATA_DIR": str(cache_dir)}
    return subprocess.run(
        [sys.executable, "-m", "tsarina.cli", *args],
        capture_output=True,
        text=True,
        check=check,
        env=env,
    )


def test_reference_help_lists_subcommands(tmp_path):
    r = _run("reference", "--help", cache_dir=tmp_path)
    assert r.returncode == 0
    for sub in ("list", "fetch", "path"):
        assert sub in r.stdout


def test_reference_bare_defaults_to_list(tmp_path):
    # `tsarina reference` with no subcommand should list (not error).
    r = _run("reference", cache_dir=tmp_path)
    assert r.returncode == 0
    assert "hpa_rna_consensus" in r.stdout
    assert "Default HPA version" in r.stdout


def test_reference_list_shows_all_datasets(tmp_path):
    r = _run("reference", "list", cache_dir=tmp_path)
    assert r.returncode == 0
    for name in ("hpa_rna_consensus", "hpa_normal_tissue", "ncbi_gene_info"):
        assert name in r.stdout


def test_reference_path_prints_cache_path(tmp_path):
    r = _run("reference", "path", "hpa_rna_consensus", cache_dir=tmp_path)
    assert r.returncode == 0
    out = r.stdout.strip()
    assert out.endswith("rna_tissue_consensus.tsv")
    assert str(tmp_path) in out


def test_reference_path_unknown_dataset_errors(tmp_path):
    r = _run("reference", "path", "does_not_exist", cache_dir=tmp_path, check=False)
    assert r.returncode != 0
    assert "unknown dataset" in (r.stderr + r.stdout)


def test_reference_path_unknown_version_errors(tmp_path):
    r = _run(
        "reference",
        "path",
        "hpa_normal_tissue",
        "--hpa-version",
        "v99",
        cache_dir=tmp_path,
        check=False,
    )
    assert r.returncode != 0
    assert "no version" in (r.stderr + r.stdout)


def test_data_sources_removed(tmp_path):
    # The old nested command must be gone (replaced by top-level `reference`).
    r = _run("data", "sources", "list", cache_dir=tmp_path, check=False)
    assert r.returncode != 0
    assert "invalid choice: 'sources'" in (r.stderr + r.stdout)


def test_reference_fetch_continues_past_failures(monkeypatch, capsys):
    """`fetch all` attempts every dataset even if one fails, then exits non-zero."""
    import pytest

    from tsarina import reference_data

    attempted = []

    def _fake_download(name, version=None, force=False):
        attempted.append(name)
        if name == "hpa_normal_tissue":
            raise reference_data.ReferenceDataError("boom")
        return f"/cache/{name}"

    monkeypatch.setattr(reference_data, "download", _fake_download)
    args = argparse.Namespace(names=[], hpa_version=None, force=False)
    with pytest.raises(SystemExit) as exc:
        cli._reference_fetch(args)
    assert exc.value.code == 1
    # All three datasets were attempted despite the middle one failing.
    assert set(attempted) == set(reference_data.REFERENCE_DATASETS)
    err = capsys.readouterr().err
    assert "hpa_normal_tissue" in err and "Failed to fetch" in err


def test_is_hpa_dataset_distinguishes_kinds():
    from tsarina import reference_data

    assert reference_data.is_hpa_dataset("hpa_rna_consensus")
    assert reference_data.is_hpa_dataset("hpa_normal_tissue")
    # The rolling NCBI dataset is not version-pinned to an HPA release.
    assert not reference_data.is_hpa_dataset("ncbi_gene_info")
    # Unknown datasets are not HPA (no spec to honor --hpa-version).
    assert not reference_data.is_hpa_dataset("does_not_exist")


def test_reference_fetch_hpa_version_only_forwarded_to_hpa_datasets(monkeypatch, capsys):
    """`fetch all --hpa-version vNN` applies the version only to HPA datasets.

    Regression: the version was forwarded to *every* dataset, so the rolling
    ``ncbi_gene_info`` (which has no vNN release) raised and aborted the run
    with exit 1 even though both HPA files downloaded fine.
    """
    from tsarina import reference_data

    seen = {}

    def _fake_download(name, version=None, force=False):
        seen[name] = version
        return f"/cache/{name}"

    monkeypatch.setattr(reference_data, "download", _fake_download)
    args = argparse.Namespace(names=[], hpa_version="v23", force=False)
    cli._reference_fetch(args)  # must not raise / exit

    assert seen["hpa_rna_consensus"] == "v23"
    assert seen["hpa_normal_tissue"] == "v23"
    assert seen["ncbi_gene_info"] is None


def test_reference_path_hpa_version_not_forwarded_to_non_hpa(monkeypatch, capsys):
    """`path <ncbi> --hpa-version vNN` ignores the version for non-HPA datasets."""
    from tsarina import reference_data

    seen = {}

    def _fake_local_path(name, version=None):
        seen[name] = version
        return f"/cache/{name}"

    monkeypatch.setattr(reference_data, "local_path", _fake_local_path)
    args = argparse.Namespace(name="ncbi_gene_info", hpa_version="v23")
    cli._reference_path(args)
    assert seen["ncbi_gene_info"] is None


def test_data_bare_defaults_to_list(monkeypatch, capsys):
    """Bare `tsarina data` routes to list (not a usage error), like reference.

    Hermetic: stubs the hitlist-delegated registry calls so the test owns its
    state instead of reading hitlist's real data dir.
    """
    monkeypatch.setattr(cli, "list_datasets", lambda: {})
    monkeypatch.setattr(cli, "data_dir", lambda: "/tmp/does-not-matter")
    cli._handle_data(argparse.Namespace(data_command=None))
    out = capsys.readouterr().out
    assert "No datasets registered" in out
    # The empty-registry branch points users at the separate reference command.
    assert "tsarina reference list" in out
