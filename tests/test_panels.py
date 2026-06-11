import pandas as pd
import pytest


def _panel_targets() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "peptide": ["PEPTIDEA", "PEPTIDEB", "PEPTIDEC"],
            "source": ["MAGEA4", "MAGEA4", "PRAME"],
            "category": ["cta", "cta", "cta"],
            "has_ms_evidence": [True, True, True],
            "ms_alleles": [
                "HLA-A*02:01;HLA-B*07:02",
                "HLA-A*24:02",
                "HLA-B*07:02",
            ],
        }
    )


def test_ms_peptide_count_uses_literal_normalized_alleles(monkeypatch):
    from tsarina.panels import build_panel_matrix

    calls = {}

    def _fake_target_peptides(**kwargs):
        calls.update(kwargs)
        return _panel_targets()

    monkeypatch.setattr("tsarina.targets.target_peptides", _fake_target_peptides, raising=True)

    matrix = build_panel_matrix(
        category="cta",
        alleles=["A*02:01", "HLA-A*24:02", "HLA-B*07:02"],
        metric="ms_peptide_count",
    )

    assert calls["attach_ms_evidence"] is True
    magea4 = matrix[matrix["source"] == "MAGEA4"].iloc[0]
    assert magea4["A*02:01"] == 1
    assert magea4["HLA-A*24:02"] == 1
    assert magea4["HLA-B*07:02"] == 1


def test_ms_confirmed_only_still_filters_target_peptides(monkeypatch):
    from tsarina.panels import build_panel_matrix

    calls = {}

    def _fake_target_peptides(**kwargs):
        calls.update(kwargs)
        return _panel_targets()

    monkeypatch.setattr("tsarina.targets.target_peptides", _fake_target_peptides, raising=True)

    build_panel_matrix(
        category="cta",
        alleles=["HLA-A*02:01"],
        metric="ms_peptide_count",
        ms_confirmed_only=True,
    )

    assert calls["require_ms_evidence"] is True
    assert calls["attach_ms_evidence"] is True


def test_best_percentile_import_error_propagates(monkeypatch):
    from tsarina.panels import build_panel_matrix

    monkeypatch.setattr("tsarina.targets.target_peptides", lambda **kw: _panel_targets())

    def _raise_import_error(**kwargs):
        raise ImportError("missing predictor")

    monkeypatch.setattr("tsarina.scoring.score_presentation", _raise_import_error, raising=True)

    with pytest.raises(ImportError, match="missing predictor"):
        build_panel_matrix(
            category="cta",
            alleles=["HLA-A*02:01"],
            metric="best_percentile",
        )


def test_best_percentile_uses_lowest_percentile_per_cell(monkeypatch):
    from tsarina.panels import build_panel_matrix

    monkeypatch.setattr("tsarina.targets.target_peptides", lambda **kw: _panel_targets())

    def _fake_score_presentation(peptides, alleles):
        return pd.DataFrame(
            {
                "peptide": ["PEPTIDEA", "PEPTIDEB", "PEPTIDEC"],
                "allele": ["HLA-A*02:01", "HLA-A*02:01", "HLA-A*02:01"],
                "presentation_percentile": [2.0, 0.4, 1.5],
            }
        )

    monkeypatch.setattr(
        "tsarina.scoring.score_presentation",
        _fake_score_presentation,
        raising=True,
    )

    matrix = build_panel_matrix(
        category="cta",
        alleles=["HLA-A*02:01"],
        metric="best_percentile",
    )

    assert set(matrix["source"]) == {"MAGEA4", "PRAME"}
    magea4 = matrix[matrix["source"] == "MAGEA4"].iloc[0]
    prame = matrix[matrix["source"] == "PRAME"].iloc[0]
    assert magea4["HLA-A*02:01"] == 0.4
    assert prame["HLA-A*02:01"] == 1.5


def test_viral_category_defaults_to_all_viruses(monkeypatch):
    from tsarina.panels import build_panel_matrix

    calls = {}

    def _fake_target_peptides(**kwargs):
        import pandas as pd

        calls.update(kwargs)
        return pd.DataFrame(columns=["source", "category", "peptide"])

    monkeypatch.setattr("tsarina.targets.target_peptides", _fake_target_peptides, raising=True)

    build_panel_matrix(category="viral", alleles=["HLA-A*02:01"])

    assert calls["viruses"] is True


def test_unknown_metric_raises():
    from tsarina.panels import build_panel_matrix

    with pytest.raises(ValueError, match="Unknown metric"):
        build_panel_matrix(metric="not_a_metric")
