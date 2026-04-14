from pathlib import Path
from unittest.mock import patch

import pandas as pd

from tsarina.indexing import ensure_index_built, load_ms_evidence


def test_ensure_index_built_skips_when_already_built(tmp_path: Path, capsys):
    fake_path = tmp_path / "obs.parquet"
    fake_path.write_text("dummy")
    with (
        patch("hitlist.observations.is_built", return_value=True),
        patch("hitlist.observations.observations_path", return_value=fake_path),
        patch("hitlist.builder.build_observations") as build,
    ):
        p = ensure_index_built()
    assert p == fake_path
    build.assert_not_called()
    # No progress message when cache hit.
    assert "Building" not in capsys.readouterr().err


def test_ensure_index_built_triggers_build_when_missing(tmp_path: Path, capsys):
    fake_path = tmp_path / "obs.parquet"
    with (
        patch("hitlist.observations.is_built", return_value=False),
        patch("hitlist.observations.observations_path", return_value=fake_path),
        patch("hitlist.builder.build_observations") as build,
    ):
        ensure_index_built()
    build.assert_called_once()
    assert "Building hitlist observations index" in capsys.readouterr().err


def test_ensure_index_built_force_triggers_rebuild(tmp_path: Path):
    fake_path = tmp_path / "obs.parquet"
    fake_path.write_text("dummy")
    with (
        patch("hitlist.observations.is_built", return_value=True),
        patch("hitlist.observations.observations_path", return_value=fake_path),
        patch("hitlist.builder.build_observations") as build,
    ):
        ensure_index_built(force=True, verbose=False)
    build.assert_called_once_with(force=True)


def _fake_loader(full_frame: pd.DataFrame):
    """Return a mock load_observations that honors the peptide= pushdown."""

    def _load(**kwargs):
        df = full_frame
        peptide = kwargs.get("peptide")
        if peptide is not None:
            wanted = set(peptide) if not isinstance(peptide, str) else {peptide}
            df = df[df["peptide"].isin(wanted)]
        return df.reset_index(drop=True)

    return _load


def test_load_ms_evidence_pushes_peptide_filter_to_loader():
    fake = pd.DataFrame(
        {
            "peptide": ["AAA", "BBB", "CCC"],
            "mhc_restriction": ["HLA-A*02:01"] * 3,
            "is_binding_assay": [False, False, True],
        }
    )
    with (
        patch("hitlist.observations.is_built", return_value=True),
        patch("hitlist.observations.load_observations", side_effect=_fake_loader(fake)) as loader,
    ):
        out = load_ms_evidence(peptides={"AAA", "CCC"})
    loader.assert_called_once()
    kwargs = loader.call_args.kwargs
    assert kwargs["mhc_class"] == "I"
    assert kwargs["species"] == "Homo sapiens"
    # Pushdown sent AAA+CCC to the loader; binding-assay drop removed CCC.
    assert sorted(kwargs["peptide"]) == ["AAA", "CCC"]
    assert out["peptide"].tolist() == ["AAA"]


def test_load_ms_evidence_falls_back_when_peptide_kwarg_unsupported():
    fake = pd.DataFrame(
        {
            "peptide": ["AAA", "BBB"],
            "mhc_restriction": ["HLA-A*02:01"] * 2,
            "is_binding_assay": [False, False],
        }
    )

    def _load(**kwargs):
        if "peptide" in kwargs:
            raise TypeError("unexpected kwarg peptide (simulating hitlist < 1.6.0)")
        return fake

    with (
        patch("hitlist.observations.is_built", return_value=True),
        patch("hitlist.observations.load_observations", side_effect=_load),
    ):
        out = load_ms_evidence(peptides={"AAA"})
    # Fallback path still produces the right answer via in-memory isin.
    assert out["peptide"].tolist() == ["AAA"]


def test_load_ms_evidence_can_skip_binding_assay_drop():
    fake = pd.DataFrame(
        {
            "peptide": ["AAA"],
            "mhc_restriction": ["HLA-A*02:01"],
            "is_binding_assay": [True],
        }
    )
    with (
        patch("hitlist.observations.is_built", return_value=True),
        patch("hitlist.observations.load_observations", side_effect=_fake_loader(fake)),
    ):
        out = load_ms_evidence(peptides={"AAA"}, drop_binding_assays=False)
    assert len(out) == 1
