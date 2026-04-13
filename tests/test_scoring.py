import pytest

from tsarina.scoring import (
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


def test_score_presentation_import_error(monkeypatch):
    """Verify that score_presentation raises ImportError with helpful
    message if topiary (the prediction backend) is not installed."""
    import builtins

    real_import = builtins.__import__

    def _block_topiary(name, *args, **kwargs):
        if name == "topiary" or name.startswith("topiary."):
            raise ImportError("No module named 'topiary'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _block_topiary)
    from tsarina.scoring import score_presentation

    with pytest.raises(ImportError, match="topiary"):
        score_presentation(peptides=["SLYNTVATL"], alleles=["HLA-A*02:01"])
