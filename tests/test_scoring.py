import sys

import pandas as pd
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


def test_mhcflurry_import_error(monkeypatch):
    """Verify that score_presentation raises ImportError with helpful
    message if MHCflurry is not installed."""
    import builtins

    import tsarina.scoring as scoring

    monkeypatch.setattr(scoring, "_MHCFLURRY_PRESENTATION_PREDICTOR", None)
    real_import = builtins.__import__

    def _block_mhcflurry(name, *args, **kwargs):
        if name == "mhcflurry" or name.startswith("mhcflurry."):
            raise ImportError("No module named 'mhcflurry'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _block_mhcflurry)
    from tsarina.scoring import score_presentation

    with pytest.raises(ImportError, match="MHCflurry is required") as excinfo:
        score_presentation(peptides=["SLYNTVATL"], alleles=["HLA-A*02:01"])
    assert sys.executable in str(excinfo.value)


def test_non_mhcflurry_predictor_import_error_mentions_topiary(monkeypatch):
    """Non-MHCflurry predictors still use topiary/mhctools."""
    import builtins

    real_import = builtins.__import__

    def _block_topiary(name, *args, **kwargs):
        if name == "topiary" or name.startswith("topiary."):
            raise ImportError("No module named 'topiary'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _block_topiary)
    from tsarina.scoring import score_presentation

    with pytest.raises(ImportError, match="Topiary is required") as excinfo:
        score_presentation(
            peptides=["SLYNTVATL"],
            alleles=["HLA-A*02:01"],
            predictor="netmhcpan",
        )
    assert sys.executable in str(excinfo.value)


def test_mhcflurry_direct_path_batches_and_skips_affinity_percentiles(monkeypatch):
    """Regression for alleles with affinity support but no affinity-percentile calibration."""
    import tsarina.scoring as scoring

    calls = []

    class FakeMHCflurryPresentationPredictor:
        def predict(
            self,
            *,
            peptides,
            alleles,
            include_affinity_percentile,
            verbose,
        ):
            calls.append(
                {
                    "peptides": list(peptides),
                    "alleles": dict(alleles),
                    "include_affinity_percentile": include_affinity_percentile,
                    "verbose": verbose,
                }
            )
            if include_affinity_percentile:
                raise ValueError("missing affinity percentile calibration")
            rows = []
            for sample_name, sample_alleles in alleles.items():
                assert sample_alleles == [sample_name]
                for peptide in peptides:
                    rows.append(
                        {
                            "peptide": peptide,
                            "sample_name": sample_name,
                            "best_allele": sample_name,
                            "affinity": 86.7 if sample_name == "HLA-C*15:05" else 35.4,
                            "presentation_score": 0.78,
                            "presentation_percentile": 0.31,
                        }
                    )
            return pd.DataFrame(rows)

    monkeypatch.setattr(
        scoring,
        "_load_mhcflurry_presentation_predictor",
        lambda: FakeMHCflurryPresentationPredictor(),
    )

    out = scoring.score_presentation(
        peptides=["SLYNTVATL", "SLYNTVATL", "GILGFVFTL"],
        alleles=["HLA-A*02:01", "HLA-C*15:05", "HLA-C*15:05"],
    )

    assert len(calls) == 1
    assert calls[0] == {
        "peptides": ["SLYNTVATL", "GILGFVFTL"],
        "alleles": {
            "HLA-A*02:01": ["HLA-A*02:01"],
            "HLA-C*15:05": ["HLA-C*15:05"],
        },
        "include_affinity_percentile": False,
        "verbose": 0,
    }
    assert list(out.columns) == [
        "peptide",
        "allele",
        "presentation_score",
        "presentation_percentile",
        "affinity_nm",
    ]
    assert set(out["allele"]) == {"HLA-A*02:01", "HLA-C*15:05"}
    assert len(out) == 4
