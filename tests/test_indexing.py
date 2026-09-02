from pathlib import Path
from unittest.mock import patch

import pandas as pd

from tsarina.indexing import ensure_index_built, load_ms_evidence


def test_ensure_index_built_skips_when_already_built(tmp_path: Path, capsys):
    fake_path = tmp_path / "obs.parquet"
    fake_path.write_text("dummy")
    with (
        patch("hitlist.observations.is_built", return_value=True),
        patch("hitlist.mappings.is_mappings_built", return_value=True),
        patch("hitlist.observations.observations_path", return_value=fake_path),
        patch("hitlist.builder.build_observations") as build,
        patch("hitlist.mappings.build_peptide_mappings") as build_mappings,
        patch("tsarina.indexing._mapping_cache_is_compatible", return_value=True),
    ):
        p = ensure_index_built()
    assert p == fake_path
    build.assert_not_called()
    build_mappings.assert_not_called()
    # No progress message when cache hit.
    assert "Building" not in capsys.readouterr().err


def test_ensure_index_built_triggers_build_when_missing(tmp_path: Path, capsys):
    fake_path = tmp_path / "obs.parquet"
    with (
        patch("hitlist.observations.is_built", return_value=False),
        patch("hitlist.mappings.is_mappings_built", return_value=False),
        patch("hitlist.observations.observations_path", return_value=fake_path),
        patch("hitlist.builder.build_observations") as build,
        patch("hitlist.mappings.mappings_path", return_value=tmp_path / "mappings.parquet"),
    ):
        ensure_index_built()
    build.assert_called_once()
    assert "Building hitlist observations index" in capsys.readouterr().err


def test_ensure_index_built_rebuilds_only_sidecar_when_it_is_missing(tmp_path: Path, capsys):
    """Observations parquet present but peptide_mappings.parquet missing
    (e.g. partial build via ``build_mappings=False``, or sidecar deleted)
    must trigger a rebuild — otherwise downstream gene-attach paths
    raise ``FileNotFoundError`` mid-query."""
    fake_path = tmp_path / "obs.parquet"
    fake_path.write_text("dummy")
    with (
        patch("hitlist.observations.is_built", return_value=True),
        patch("hitlist.mappings.is_mappings_built", return_value=False),
        patch("hitlist.observations.observations_path", return_value=fake_path),
        patch("hitlist.builder.build_observations") as build_observations,
        patch("hitlist.mappings.build_peptide_mappings") as build_mappings,
        patch("hitlist.mappings.mappings_path", return_value=tmp_path / "mappings.parquet"),
    ):
        ensure_index_built()
    build_observations.assert_not_called()
    build_mappings.assert_called_once_with(force=True)
    assert "Building hitlist peptide mappings (missing" in capsys.readouterr().err


def test_ensure_index_built_rebuilds_stale_sidecar(tmp_path: Path, capsys):
    fake_path = tmp_path / "obs.parquet"
    fake_path.write_text("dummy")
    with (
        patch("hitlist.observations.is_built", return_value=True),
        patch("hitlist.mappings.is_mappings_built", return_value=True),
        patch("hitlist.observations.observations_path", return_value=fake_path),
        patch("hitlist.builder.build_observations") as build_observations,
        patch("hitlist.mappings.build_peptide_mappings") as build_mappings,
        patch("hitlist.mappings.mappings_path", return_value=tmp_path / "mappings.parquet"),
        patch("tsarina.indexing._mapping_cache_is_compatible", return_value=False),
    ):
        ensure_index_built()
    build_observations.assert_not_called()
    build_mappings.assert_called_once_with(force=True)
    assert "predates hitlist 1.55.2" in capsys.readouterr().err


def test_ensure_index_built_force_triggers_rebuild(tmp_path: Path):
    fake_path = tmp_path / "obs.parquet"
    fake_path.write_text("dummy")
    with (
        patch("hitlist.observations.is_built", return_value=True),
        patch("hitlist.mappings.is_mappings_built", return_value=True),
        patch("hitlist.observations.observations_path", return_value=fake_path),
        patch("hitlist.builder.build_observations") as build,
        patch("hitlist.mappings.mappings_path", return_value=tmp_path / "mappings.parquet"),
    ):
        ensure_index_built(force=True, verbose=False)
    build.assert_called_once_with(force=True)


def test_mapping_compatibility_probe_accepts_current_sidecar(tmp_path: Path):
    from tsarina import indexing

    mapping_path = tmp_path / "peptide_mappings.parquet"
    mapping_path.write_text("current")
    class_ii = pd.DataFrame({"peptide": ["A" * 15, "B" * 15]})
    length_7 = pd.DataFrame({"peptide": ["C" * 7]})

    def _mapped(*, peptide, columns):
        assert columns == ["peptide"]
        return pd.DataFrame({"peptide": peptide})

    with (
        patch("hitlist.observations.load_observations", side_effect=[class_ii, length_7]),
        patch("hitlist.mappings.load_peptide_mappings", side_effect=_mapped),
    ):
        assert indexing._probe_mapping_compatibility(mapping_path)

    marker = indexing._mapping_compatibility_path(mapping_path)
    assert marker.exists()
    assert indexing._mapping_compatibility_is_recorded(mapping_path)


def test_mapping_compatibility_probe_rejects_pre_1_55_sidecar(tmp_path: Path):
    from tsarina import indexing

    mapping_path = tmp_path / "peptide_mappings.parquet"
    mapping_path.write_text("stale")
    class_ii = pd.DataFrame({"peptide": ["A" * 15]})
    with (
        patch("hitlist.observations.load_observations", return_value=class_ii),
        patch(
            "hitlist.mappings.load_peptide_mappings",
            return_value=pd.DataFrame({"peptide": pd.Series(dtype=str)}),
        ),
    ):
        assert not indexing._probe_mapping_compatibility(mapping_path)

    assert not indexing._mapping_compatibility_path(mapping_path).exists()


def test_recorded_mapping_compatibility_is_invalidated_by_file_change(tmp_path: Path):
    from tsarina import indexing

    mapping_path = tmp_path / "peptide_mappings.parquet"
    mapping_path.write_text("first")
    indexing._record_mapping_compatibility(mapping_path)
    assert indexing._mapping_compatibility_is_recorded(mapping_path)

    mapping_path.write_text("replacement-with-a-different-size")
    assert not indexing._mapping_compatibility_is_recorded(mapping_path)


def test_recorded_mapping_compatibility_skips_semantic_probe(tmp_path: Path):
    from tsarina import indexing

    mapping_path = tmp_path / "peptide_mappings.parquet"
    mapping_path.write_text("verified")
    indexing._record_mapping_compatibility(mapping_path)

    with (
        patch("hitlist.mappings.mappings_path", return_value=mapping_path),
        patch("tsarina.indexing._probe_mapping_compatibility") as probe,
    ):
        assert indexing._mapping_cache_is_compatible()
    probe.assert_not_called()


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


def test_load_ms_evidence_validates_mapping_cache_for_gene_filter():
    fake = pd.DataFrame(
        {
            "peptide": ["AAAAAAAAAAAAAAA"],
            "mhc_restriction": ["HLA-DRB1*04:01"],
            "is_binding_assay": [False],
        }
    )
    with (
        patch("hitlist.observations.is_built", return_value=True),
        patch("tsarina.indexing.ensure_index_built") as ensure,
        patch("hitlist.observations.load_observations", side_effect=_fake_loader(fake)),
    ):
        out = load_ms_evidence(gene_name="MAGEA4", mhc_class="II")
    ensure.assert_called_once_with()
    assert len(out) == 1
