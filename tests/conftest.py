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

The autouse fixture in this file installs a tiny in-memory stand-in for
``mhcflurry.Class1PresentationPredictor`` into the cached singleton at
``tsarina.scoring._MHCFLURRY_PRESENTATION_PREDICTOR``. This prevents the
real predictor — which loads ~140 Keras models and consumes 3-6 GB per
process — from ever being instantiated during the test suite.

See issue #68: without this stub, any test path that reaches
``score_presentation(..., predictor="mhcflurry")`` without per-test
monkeypatching can trigger an OOM (~55 GB observed) when multiplied
across pytest-cov subprocesses.

Tests that need to exercise the real loader path (e.g. asserting the
ImportError raised when mhcflurry isn't installed) reset the singleton
to ``None`` themselves; the loader function is unchanged and will
attempt the real import when the singleton is cleared.
"""

from __future__ import annotations

import pandas as pd
import pytest


class _StubMHCflurryPresentationPredictor:
    """Deterministic stand-in for ``Class1PresentationPredictor``.

    Mirrors the keyword-only ``predict`` surface that
    ``tsarina.scoring._score_mhcflurry`` invokes and returns a DataFrame
    with the columns the downstream code reads.
    """

    def predict(self, *, peptides, alleles, include_affinity_percentile, verbose):
        rows = []
        for sample_name, sample_alleles in alleles.items():
            assert sample_alleles == [sample_name]
            for peptide in peptides:
                rows.append(
                    {
                        "peptide": peptide,
                        "sample_name": sample_name,
                        "best_allele": sample_name,
                        "affinity": 50.0,
                        "presentation_score": 0.95,
                        "presentation_percentile": 0.5,
                    }
                )
        return pd.DataFrame(rows)


@pytest.fixture(autouse=True)
def _stub_mhcflurry_predictor(monkeypatch):
    import tsarina.scoring as scoring

    monkeypatch.setattr(
        scoring,
        "_MHCFLURRY_PRESENTATION_PREDICTOR",
        _StubMHCflurryPresentationPredictor(),
        raising=True,
    )
