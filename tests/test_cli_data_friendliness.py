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

import os
import subprocess
import sys

from tsarina.version import __version__


def _run(*args, cache_dir):
    env = {**os.environ, "TSARINA_DATA_DIR": str(cache_dir)}
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
