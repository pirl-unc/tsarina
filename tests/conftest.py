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

"""Shared pytest fixtures.

The autouse fixture in this file installs a sentinel object into the
cached predictor singleton at
``tsarina.scoring._MHCFLURRY_PRESENTATION_PREDICTOR``. The sentinel
prevents the real ``Class1PresentationPredictor`` — which loads ~140
Keras models and consumes 3-6 GB per process — from ever being
instantiated during the test suite, and refuses to fabricate
predictions if a test path reaches it.

See issue #68: without this stub, any test path that reaches
``score_presentation(..., predictor="mhcflurry")`` without per-test
monkeypatching loaded real MHCflurry. Multiplied across pytest-xdist
workers this produced an OOM (~55 GB observed) on developer machines.

Per-test stubs of ``tsarina.scoring.score_presentation`` or
``tsarina.scoring._load_mhcflurry_presentation_predictor`` continue to
override this sentinel at the call site, so tests with their own scoring
fakes are unaffected. Tests that need to exercise the real loader path
(e.g. asserting the ImportError raised when mhcflurry isn't installed)
reset the singleton to ``None`` themselves.
"""

from __future__ import annotations

import pytest


class _UnstubbedMHCflurryPredictor:
    """Sentinel ``Class1PresentationPredictor`` that refuses to predict.

    Installed into ``tsarina.scoring._MHCFLURRY_PRESENTATION_PREDICTOR`` by
    the autouse fixture below so that no test path can accidentally load
    real MHCflurry. If a test reaches this object's ``predict`` method,
    that's a bug in the test (it forgot to stub
    ``tsarina.scoring.score_presentation`` or
    ``tsarina.scoring._load_mhcflurry_presentation_predictor``) — fail
    loudly rather than return dummy values that could be silently
    asserted on.
    """

    def predict(self, *, peptides, alleles, include_affinity_percentile, verbose):
        raise AssertionError(
            "Real MHCflurry predictor was reached during the test suite. "
            "Stub `tsarina.scoring.score_presentation` (preferred) or "
            "`tsarina.scoring._load_mhcflurry_presentation_predictor` in "
            "the test that triggered this call. To exercise the real "
            "loader path, reset the singleton to None first: "
            "`monkeypatch.setattr(scoring, '_MHCFLURRY_PRESENTATION_PREDICTOR', None)`."
        )


@pytest.fixture(autouse=True)
def _stub_mhcflurry_predictor(monkeypatch):
    import tsarina.scoring as scoring

    monkeypatch.setattr(
        scoring,
        "_MHCFLURRY_PRESENTATION_PREDICTOR",
        _UnstubbedMHCflurryPredictor(),
        raising=True,
    )
