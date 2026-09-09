"""Tests for tsarina.indexing.

``ensure_index_built`` owns one policy: ask hitlist whether either artifact
needs rebuilding.  These tests pin that delegation, the stderr reporting, and
the fact that hitlist's own stdout never reaches tsarina's data stream.
"""

from pathlib import Path
from unittest.mock import patch

import pandas as pd

from tsarina.indexing import ensure_index_built, load_ms_evidence


def _artifact_patches(tmp_path: Path, observations: Path, *, is_built: bool):
    """Patch set common to every ensure_index_built test."""
    return (
        patch("hitlist.observations.is_built", return_value=is_built),
        patch("hitlist.observations.observations_path", return_value=observations),
        patch("hitlist.mappings.mappings_path", return_value=tmp_path / "mappings.parquet"),
    )


def test_ensure_index_built_delegates_freshness_to_hitlist(tmp_path: Path, capsys):
    """An existing artifact is validated, not assumed current.

    Existence is not freshness: hitlist compares its curation fingerprints and
    artifact contracts, so tsarina must hand the decision over on every call.
    """
    observations = tmp_path / "obs.parquet"
    observations.write_text("dummy")
    is_built_p, obs_path_p, map_path_p = _artifact_patches(tmp_path, observations, is_built=True)
    with (
        is_built_p,
        obs_path_p,
        map_path_p,
        patch("tsarina.indexing._sources_registered", return_value=True),
        patch("hitlist.builder.build_observations") as build,
    ):
        assert ensure_index_built() == observations
    build.assert_called_once_with(force=False)
    captured = capsys.readouterr()
    assert captured.err == ""
    assert captured.out == ""


def test_ensure_index_built_can_report_a_current_index(tmp_path: Path, capsys):
    """``tsarina build observations`` exists to describe the index."""
    observations = tmp_path / "obs.parquet"
    observations.write_text("dummy")
    is_built_p, obs_path_p, map_path_p = _artifact_patches(tmp_path, observations, is_built=True)
    with (
        is_built_p,
        obs_path_p,
        map_path_p,
        patch("tsarina.indexing._sources_registered", return_value=True),
        patch(
            "hitlist.builder.build_observations",
            side_effect=lambda force=False: print("Observations already built (40 rows)"),
        ),
    ):
        ensure_index_built(report_current=True)
    captured = capsys.readouterr()
    assert "Observations already built (40 rows)" in captured.err
    assert "stale" not in captured.err
    assert captured.out == ""


def test_ensure_index_built_reports_a_stale_rebuild(tmp_path: Path, capsys):
    observations = tmp_path / "obs.parquet"
    observations.write_text("legacy")

    def _rebuild(force: bool = False):
        print("Observations rebuilt (4,439,321 rows)")
        observations.write_text("rebuilt with the current curation")

    is_built_p, obs_path_p, map_path_p = _artifact_patches(tmp_path, observations, is_built=True)
    with (
        is_built_p,
        obs_path_p,
        map_path_p,
        patch("tsarina.indexing._sources_registered", return_value=True),
        patch("hitlist.builder.build_observations", side_effect=_rebuild),
    ):
        ensure_index_built()
    captured = capsys.readouterr()
    assert "stale" in captured.err
    assert "Observations rebuilt (4,439,321 rows)" in captured.err
    # hitlist reports on stdout; tsarina's stdout carries query results.
    assert captured.out == ""


def test_ensure_index_built_keeps_hitlist_output_off_stdout(tmp_path: Path, capsys):
    observations = tmp_path / "obs.parquet"
    is_built_p, obs_path_p, map_path_p = _artifact_patches(tmp_path, observations, is_built=False)
    with (
        is_built_p,
        obs_path_p,
        map_path_p,
        patch(
            "hitlist.builder.build_observations",
            side_effect=lambda force=False: print("scanning IEDB rows..."),
        ),
    ):
        ensure_index_built()
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "scanning IEDB rows..." in captured.err


def test_ensure_index_built_triggers_build_when_missing(tmp_path: Path, capsys):
    observations = tmp_path / "obs.parquet"
    is_built_p, obs_path_p, map_path_p = _artifact_patches(tmp_path, observations, is_built=False)
    with (
        is_built_p,
        obs_path_p,
        map_path_p,
        patch("hitlist.builder.build_observations") as build,
    ):
        ensure_index_built()
    build.assert_called_once_with(force=False)
    assert "Building hitlist observations index" in capsys.readouterr().err


def test_ensure_index_built_missing_index_does_not_need_registered_sources(tmp_path: Path):
    """A missing index always builds; hitlist raises its own registration error."""
    observations = tmp_path / "obs.parquet"
    is_built_p, obs_path_p, map_path_p = _artifact_patches(tmp_path, observations, is_built=False)
    with (
        is_built_p,
        obs_path_p,
        map_path_p,
        patch("tsarina.indexing._sources_registered") as registered,
        patch("hitlist.builder.build_observations") as build,
    ):
        ensure_index_built(verbose=False)
    build.assert_called_once_with(force=False)
    registered.assert_not_called()


def test_ensure_index_built_force_triggers_rebuild(tmp_path: Path):
    observations = tmp_path / "obs.parquet"
    observations.write_text("dummy")
    is_built_p, obs_path_p, map_path_p = _artifact_patches(tmp_path, observations, is_built=True)
    with (
        is_built_p,
        obs_path_p,
        map_path_p,
        patch("hitlist.builder.build_observations") as build,
    ):
        ensure_index_built(force=True, verbose=False)
    build.assert_called_once_with(force=True)


def test_ensure_index_built_uses_prebuilt_index_without_registered_sources(tmp_path: Path):
    """A prebuilt index copied in without the raw exports stays usable.

    ``build_observations`` raises when no IEDB/CEDAR export is registered, so
    validation is impossible rather than merely inconvenient here.
    """
    observations = tmp_path / "obs.parquet"
    observations.write_text("dummy")
    is_built_p, obs_path_p, map_path_p = _artifact_patches(tmp_path, observations, is_built=True)
    with (
        is_built_p,
        obs_path_p,
        map_path_p,
        # tsarina.datasources binds get_path at import time.
        patch("tsarina.datasources.get_path", side_effect=KeyError("iedb")),
        patch("hitlist.builder.build_observations") as build,
    ):
        assert ensure_index_built() == observations
    build.assert_not_called()


def test_sources_registered_accepts_cedar_alone():
    from tsarina import indexing

    def _get_path(name):
        if name == "cedar":
            return "/tmp/cedar.csv"
        raise FileNotFoundError(name)

    with patch("tsarina.datasources.get_path", side_effect=_get_path):
        assert indexing._sources_registered()


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
        patch("tsarina.indexing.ensure_index_built") as ensure,
        patch("hitlist.observations.load_observations", side_effect=_fake_loader(fake)) as loader,
    ):
        out = load_ms_evidence(peptides={"AAA", "CCC"})
    ensure.assert_called_once_with()
    loader.assert_called_once()
    kwargs = loader.call_args.kwargs
    assert kwargs["mhc_class"] == "I"
    assert kwargs["species"] == "Homo sapiens"
    # Pushdown sent AAA+CCC to the loader; binding-assay drop removed CCC.
    assert sorted(kwargs["peptide"]) == ["AAA", "CCC"]
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
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_observations", side_effect=_fake_loader(fake)),
    ):
        out = load_ms_evidence(peptides={"AAA"}, drop_binding_assays=False)
    assert len(out) == 1


def test_load_ms_evidence_validates_the_index_for_a_gene_filter():
    fake = pd.DataFrame(
        {
            "peptide": ["AAAAAAAAAAAAAAA"],
            "mhc_restriction": ["HLA-DRB1*04:01"],
            "is_binding_assay": [False],
        }
    )
    with (
        patch("tsarina.indexing.ensure_index_built") as ensure,
        patch("hitlist.observations.load_observations", side_effect=_fake_loader(fake)),
    ):
        out = load_ms_evidence(gene_name="MAGEA4", mhc_class="II")
    ensure.assert_called_once_with()
    assert len(out) == 1


def test_load_ms_evidence_can_opt_out_of_validation():
    """``auto_build=False`` is the one path that reads whatever is on disk."""
    fake = pd.DataFrame(
        {
            "peptide": ["AAA"],
            "mhc_restriction": ["HLA-A*02:01"],
            "is_binding_assay": [False],
        }
    )
    with (
        patch("tsarina.indexing.ensure_index_built") as ensure,
        patch("hitlist.observations.load_observations", side_effect=_fake_loader(fake)),
    ):
        out = load_ms_evidence(peptides={"AAA"}, auto_build=False)
    ensure.assert_not_called()
    assert len(out) == 1
