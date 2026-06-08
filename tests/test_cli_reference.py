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

"""CLI smoke tests for ``tsarina reference`` and the bare-``data`` default.

Exercises the parser wiring + dispatch end-to-end via subprocess (no network:
``list``/``path`` only read the manifest / compute paths). The cache is
redirected to a tmp dir so the tests never touch the real cache or download.
"""

import os
import subprocess
import sys


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


def test_data_bare_defaults_to_list(tmp_path):
    # `tsarina data` with no subcommand should list (not error), like reference.
    r = _run("data", cache_dir=tmp_path)
    assert r.returncode == 0
    # Either the registered table or the empty-registry hint, but never a usage error.
    assert "Usage: tsarina data" not in (r.stderr + r.stdout)
