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

"""CLI friendliness: --version, actionable unknown-dataset errors, and a
non-silent `data remove` on a typo."""

import json
import os
import subprocess
import sys

from tsarina.version import __version__


def _run(*args, cache_dir):
    env = {**os.environ, "HITLIST_DATA_DIR": str(cache_dir)}
    return subprocess.run(
        [sys.executable, "-m", "tsarina.cli", *args],
        capture_output=True,
        text=True,
        env=env,
    )


def test_version_flag(tmp_path):
    r = _run("--version", cache_dir=tmp_path)
    assert r.returncode == 0
    assert __version__ in r.stdout


def test_data_path_unknown_is_actionable(tmp_path):
    r = _run("data", "path", "definitely-not-a-dataset", cache_dir=tmp_path)
    assert r.returncode == 1
    # Friendly message (lists the catalog), not a leaked KeyError repr.
    assert "unknown dataset" in r.stderr.lower()
    assert "data available" in r.stderr


def test_data_remove_typo_is_not_silent_success(tmp_path):
    r = _run("data", "remove", "never-registered-xyz", cache_dir=tmp_path)
    assert r.returncode == 1
    assert "nothing to do" in r.stderr.lower()


def test_unknown_dataset_error_is_uniform_across_subcommands(tmp_path):
    for sub in ("path", "info", "fetch", "refresh"):
        r = _run("data", sub, "iedbb", cache_dir=tmp_path)
        assert r.returncode == 1, sub
        assert "unknown dataset 'iedbb'" in r.stderr, sub
        assert "Did you mean: iedb" in r.stderr, sub
        assert "Traceback" not in r.stderr, sub
        assert "[" not in r.stderr, sub  # no Python list repr of the catalog


def test_fetch_dry_run_downloads_nothing(tmp_path):
    r = _run("data", "fetch", "hpv16", "--dry-run", cache_dir=tmp_path)
    assert r.returncode == 0, r.stderr
    assert "Source:" in r.stdout and "Free space:" in r.stdout
    assert "would download" in r.stdout
    assert list(tmp_path.iterdir()) == []


def test_list_flags_missing_registered_file_and_emits_json(tmp_path):
    f = tmp_path / "ligands.csv"
    f.write_text("a,b\n1,2\n")
    assert _run("data", "register", "iedb", str(f), cache_dir=tmp_path).returncode == 0
    f.unlink()
    r = _run("data", "list", cache_dir=tmp_path)
    assert r.returncode == 0
    assert "missing" in r.stdout
    rows = json.loads(_run("data", "list", "--json", cache_dir=tmp_path).stdout)
    assert [(e["name"], e["status"]) for e in rows] == [("iedb", "missing")]
    available = json.loads(_run("data", "available", "--json", cache_dir=tmp_path).stdout)
    assert {"name": "iedb", "status": "missing"}.items() <= next(
        e for e in available if e["name"] == "iedb"
    ).items()


def test_info_reports_whether_file_exists(tmp_path):
    f = tmp_path / "ligands.csv"
    f.write_text("a,b\n1,2\n")
    _run("data", "register", "iedb", str(f), cache_dir=tmp_path)
    assert json.loads(_run("data", "info", "iedb", cache_dir=tmp_path).stdout)["exists"] is True
    f.unlink()
    assert json.loads(_run("data", "info", "iedb", cache_dir=tmp_path).stdout)["exists"] is False


def test_remove_delete_names_the_file(tmp_path):
    f = tmp_path / "ligands.csv"
    f.write_text("a,b\n1,2\n")
    _run("data", "register", "iedb", str(f), cache_dir=tmp_path)
    r = _run("data", "remove", "iedb", "--delete", cache_dir=tmp_path)
    assert r.returncode == 0
    assert str(f) in r.stderr and "other projects" in r.stderr
    assert not f.exists()
