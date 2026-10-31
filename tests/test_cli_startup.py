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

"""``--version`` must not import oncoref; the lazy public API must still resolve."""

import subprocess
import sys

import pytest

import tsarina


def test_cli_trivial_commands_skip_heavy_imports():
    code = (
        "import sys; sys.argv = ['tsarina', '--version']\n"
        "import tsarina.cli as c\n"
        "try:\n    c.main()\nexcept SystemExit:\n    pass\n"
        "assert 'oncoref' not in sys.modules\n"  # pandas is hitlist's to import, not ours
    )
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_every_public_name_resolves_lazily():
    for name in tsarina.__all__:
        assert getattr(tsarina, name) is not None, name
    assert tsarina.personalize.__name__ == "tsarina.personalize"  # submodule, not the function
    assert callable(tsarina.personalized_targets)


def test_unknown_attribute_raises_attribute_error():
    with pytest.raises(AttributeError):
        tsarina.definitely_not_a_name  # noqa: B018
