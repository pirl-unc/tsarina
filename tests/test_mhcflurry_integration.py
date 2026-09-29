"""Opt-in application contracts against the installed, downloaded model bundle.

CI installs the declared extra and downloads models before enabling these tests.
Missing dependencies or models must fail when enabled, never silently skip.
"""

import numpy as np
import pytest

from tsarina import scoring
from tsarina.alleles import get_panel

pytestmark = pytest.mark.mhcflurry


@pytest.fixture(scope="module")
def presentation_predictor():
    from mhcflurry import Class1PresentationPredictor

    return Class1PresentationPredictor.load()


def test_panel_affinity_calibration(presentation_predictor):
    predictor = presentation_predictor.affinity_predictor
    # Inspect the affinity model bundled with the presentation predictor that
    # tsarina actually uses, rather than an independently installed pan model.
    alleles = get_panel("global54_abc")
    for affinity in [50.0, 500.0, 5000.0]:
        ranks = predictor.percentile_ranks([affinity] * len(alleles), alleles=alleles)
        assert np.isfinite(ranks).all()
        assert ((ranks >= 0) & (ranks <= 100)).all()
    # No assertion about equivalence between C*14:02 and C*14:03: that depends
    # on the bundle's pseudosequences and calibration, not tsarina's contract.


def test_presentation_scores_cover_requested_pairs(presentation_predictor, monkeypatch):
    monkeypatch.setattr(scoring, "_MHCFLURRY_PRESENTATION_PREDICTOR", presentation_predictor)
    peptides = ["SIINFEKL", "SLLMWITQC"]
    # C*15:05 also exercises scoring without requiring affinity calibration.
    alleles = [*get_panel("global54_abc"), "HLA-C*15:05"]
    result = scoring.score_presentation(peptides, alleles)
    assert len(result) == len(peptides) * len(alleles)
    assert set(zip(result.peptide, result.allele)) == {
        (peptide, allele) for peptide in peptides for allele in alleles
    }
    for column in ["affinity_nm", "presentation_score", "presentation_percentile"]:
        assert np.isfinite(result[column]).all(), column
    assert (result.affinity_nm > 0).all()
    assert result.presentation_score.between(0, 1).all()
    assert result.presentation_percentile.between(0, 100).all()
    assert result.affinity_percentile.isna().all()
