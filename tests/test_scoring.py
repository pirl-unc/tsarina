import pytest

from perseus.scoring import (
    AFFINITY_THRESHOLDS_NM,
    PRESENTATION_PERCENTILE_THRESHOLDS,
    PRESENTATION_SCORE_THRESHOLDS,
)


def test_affinity_thresholds_are_sorted():
    assert list(AFFINITY_THRESHOLDS_NM) == sorted(AFFINITY_THRESHOLDS_NM)


def test_presentation_score_thresholds_are_sorted():
    assert list(PRESENTATION_SCORE_THRESHOLDS) == sorted(PRESENTATION_SCORE_THRESHOLDS)


def test_presentation_percentile_thresholds_are_sorted():
    assert list(PRESENTATION_PERCENTILE_THRESHOLDS) == sorted(PRESENTATION_PERCENTILE_THRESHOLDS)


def test_score_presentation_import_error():
    """Verify that score_presentation raises ImportError with helpful message
    if mhcflurry is not installed."""
    try:
        import mhcflurry  # noqa: F401

        pytest.skip("mhcflurry is installed")
    except ImportError:
        pass

    from perseus.scoring import score_presentation

    with pytest.raises(ImportError, match="mhcflurry"):
        score_presentation(peptides=["SLYNTVATL"], alleles=["HLA-A*02:01"])
