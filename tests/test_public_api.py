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

The personalize workflow is exported as `personalized_targets` (a distinct name
that doesn't shadow the `tsarina.personalize` submodule); `personalize` remains a
back-compat alias inside that module.
"""

import tsarina


def test_headline_verbs_importable():
    from tsarina import (  # noqa: F401
        HOTSPOT_MUTATIONS,
        mutant_peptides,
        personalized_targets,
        score_presentation,
        target_peptides,
        viral_peptides,
    )

    assert callable(personalized_targets)
    assert callable(target_peptides)
    assert callable(mutant_peptides)


def test_headline_verbs_in_all():
    for name in (
        "personalized_targets",
        "target_peptides",
        "mutant_peptides",
        "score_presentation",
        "viral_peptides",
    ):
        assert name in tsarina.__all__


def test_personalize_alias_still_works():
    # Back-compat: the old submodule import resolves to the renamed function.
    from tsarina import personalized_targets
    from tsarina.personalize import personalize

    assert personalize is personalized_targets
