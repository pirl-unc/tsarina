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

"""The headline analysis verbs are importable from the top level.

(`personalize` is the exception — it stays at `tsarina.personalize.personalize`
to avoid shadowing the same-named submodule.)"""

import tsarina


def test_headline_verbs_importable():
    from tsarina import (  # noqa: F401
        HOTSPOT_MUTATIONS,
        mutant_peptides,
        score_presentation,
        target_peptides,
        viral_peptides,
    )

    assert callable(target_peptides)
    assert callable(mutant_peptides)


def test_headline_verbs_in_all():
    for name in (
        "target_peptides",
        "mutant_peptides",
        "score_presentation",
        "viral_peptides",
    ):
        assert name in tsarina.__all__


def test_personalize_importable_from_submodule():
    from tsarina.personalize import personalize

    assert callable(personalize)
