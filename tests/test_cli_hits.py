import argparse
import subprocess
import sys
from unittest.mock import patch

import pandas as pd


def _run_cli(*args: str, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "tsarina.cli", *args],
        capture_output=True,
        text=True,
        check=check,
    )


def test_hits_help_exits_zero():
    r = _run_cli("hits", "--help")
    assert r.returncode == 0
    assert "--gene" in r.stdout
    assert "--allele" in r.stdout
    assert "--serotype" in r.stdout
    assert "--species" in r.stdout
    assert "--format" in r.stdout
    assert "--predict" in r.stdout
    assert "--mono-allelic-only" in r.stdout


def test_hits_help_describes_include_binding_assays_union():
    """tsarina#4 closed: cached path now routes through hitlist's
    load_all_evidence() union; help text should describe the UNION
    behavior and the evidence_kind tagging, not the old no-op limitation."""
    r = _run_cli("hits", "--help")
    assert r.returncode == 0
    assert "--include-binding-assays" in r.stdout
    assert "evidence_kind" in r.stdout


def test_hits_format_refs_is_listed_as_choice():
    """The argparse ``--format`` help text should advertise ``refs`` as
    a valid choice alongside the older peptides / pmhc / raw options."""
    r = _run_cli("hits", "--help")
    assert r.returncode == 0
    # argparse renders choices inline with the flag, e.g. "{peptides,pmhc,refs,raw}"
    assert "refs" in r.stdout
    # Sibling options should still be advertised
    assert "peptides" in r.stdout
    assert "pmhc" in r.stdout


def test_hits_format_refs_rejects_unknown_value():
    r = _run_cli("hits", "--gene", "PRAME", "--format", "not-a-format", check=False)
    assert r.returncode != 0
    combined = r.stderr + r.stdout
    # argparse error message surfaces the valid choices on invalid input
    assert "refs" in combined and "peptides" in combined


def test_hits_requires_gene_or_uniprot():
    r = _run_cli("hits", check=False)
    assert r.returncode != 0
    combined = r.stderr + r.stdout
    assert "--gene" in combined or "--uniprot" in combined


def test_hits_gene_and_uniprot_are_mutually_exclusive():
    r = _run_cli(
        "hits",
        "--gene",
        "PRAME",
        "--uniprot",
        "P08819",
        check=False,
    )
    assert r.returncode != 0


def test_include_binding_assays_routes_to_load_all_evidence(tmp_path):
    """tsarina#4 closed via hitlist 1.10.0: --include-binding-assays on the
    cached path must call load_all_evidence (union of MS + binding), not
    load_observations (MS-only).  This is the contract that makes the flag
    actually honor its name."""
    from tsarina import cli_hits

    args = argparse.Namespace(
        gene="PRAME",
        uniprot=None,
        allele=[],
        serotype=[],
        species="Homo sapiens",
        mhc_class="I",
        min_resolution=None,
        lengths=(8, 9, 10, 11),
        ensembl_release=112,
        include_binding_assays=True,
        mono_allelic_only=False,
        format="pmhc",
        predict=False,
        predictor="mhcflurry",
        iedb_path=None,
        cedar_path=None,
        skip_ms_evidence=False,
        output=str(tmp_path / "out.csv"),
    )
    empty = pd.DataFrame(
        {
            "peptide": pd.Series(dtype=str),
            "mhc_restriction": pd.Series(dtype=str),
            "evidence_kind": pd.Series(dtype=str),
        }
    )
    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_all_evidence", return_value=empty) as all_ev,
        patch("hitlist.observations.load_observations") as ms_only,
    ):
        cli_hits.handle(args)
    all_ev.assert_called_once()
    ms_only.assert_not_called()


def test_default_cached_path_uses_load_observations_not_union(tmp_path):
    """Regression guard: without --include-binding-assays, the cached path
    must still use load_observations (MS-only), not the union."""
    from tsarina import cli_hits

    args = argparse.Namespace(
        gene="PRAME",
        uniprot=None,
        allele=[],
        serotype=[],
        species="Homo sapiens",
        mhc_class="I",
        min_resolution=None,
        lengths=(8, 9, 10, 11),
        ensembl_release=112,
        include_binding_assays=False,
        mono_allelic_only=False,
        format="pmhc",
        predict=False,
        predictor="mhcflurry",
        iedb_path=None,
        cedar_path=None,
        skip_ms_evidence=False,
        output=str(tmp_path / "out.csv"),
    )
    empty = pd.DataFrame({"peptide": pd.Series(dtype=str), "mhc_restriction": pd.Series(dtype=str)})
    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_observations", return_value=empty) as ms_only,
        patch("hitlist.observations.load_all_evidence") as all_ev,
    ):
        cli_hits.handle(args)
    ms_only.assert_called_once()
    all_ev.assert_not_called()


def _base_cached_args(tmp_path, lengths, *, include_binding_assays=False, mhc_class="I"):
    return argparse.Namespace(
        gene="PRAME",
        uniprot=None,
        allele=[],
        serotype=[],
        species="Homo sapiens",
        mhc_class=mhc_class,
        min_resolution=None,
        lengths=lengths,
        ensembl_release=112,
        include_binding_assays=include_binding_assays,
        mono_allelic_only=False,
        format="pmhc",
        predict=False,
        predictor="mhcflurry",
        iedb_path=None,
        cedar_path=None,
        skip_ms_evidence=False,
        output=str(tmp_path / "out.csv"),
    )


def test_cached_path_pushes_length_bounds_to_load_observations(tmp_path):
    """hitlist 1.15.1+ exposes length_min/length_max on load_observations.
    tsarina's --lengths flag controls enumeration; on the cached path those
    bounds must also push down so we don't silently pull e.g. 13-mer MHC-II
    rows when the user asked for 8-11."""
    from tsarina import cli_hits

    args = _base_cached_args(tmp_path, lengths=(8, 9, 10, 11))
    empty = pd.DataFrame({"peptide": pd.Series(dtype=str), "mhc_restriction": pd.Series(dtype=str)})
    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_observations", return_value=empty) as ms_only,
    ):
        cli_hits.handle(args)
    kwargs = ms_only.call_args.kwargs
    assert kwargs.get("length_min") == 8
    assert kwargs.get("length_max") == 11


def test_cached_path_pushes_length_bounds_to_load_all_evidence(tmp_path):
    """Parity with load_observations: the union path must also forward
    length_min/length_max derived from --lengths.  Use an MHC-II window
    so the class / length pairing matches realistic usage."""
    from tsarina import cli_hits

    args = _base_cached_args(
        tmp_path,
        lengths=(12, 13, 14, 15),
        include_binding_assays=True,
        mhc_class="II",
    )
    empty = pd.DataFrame(
        {
            "peptide": pd.Series(dtype=str),
            "mhc_restriction": pd.Series(dtype=str),
            "evidence_kind": pd.Series(dtype=str),
        }
    )
    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_all_evidence", return_value=empty) as all_ev,
    ):
        cli_hits.handle(args)
    kwargs = all_ev.call_args.kwargs
    assert kwargs.get("length_min") == 12
    assert kwargs.get("length_max") == 15


def test_cached_path_non_contiguous_lengths_exact_set_filter(tmp_path):
    """--lengths 9,11 (non-contiguous) pushes a [9, 11] bound to hitlist
    but must post-filter to the exact set so 10-mers don't leak through.
    Verifies both the bound-forwarding and the post-filter independently
    so a refactor that drops one and keeps the other still fails a test."""
    from tsarina import cli_hits

    args = _base_cached_args(tmp_path, lengths=(9, 11))
    hits = pd.DataFrame(
        {
            "peptide": ["AAAAAAAAA", "AAAAAAAAAA", "AAAAAAAAAAA"],  # 9, 10, 11
            "mhc_restriction": ["HLA-A*02:01"] * 3,
        }
    )
    # Post-1.30.46 the fast path re-attaches gene columns via
    # ``load_peptide_mappings``; the test focus is length filtering, so
    # provide an empty-mappings stub to keep the loader real-import-free.
    empty_mappings = pd.DataFrame(
        {
            "peptide": pd.Series(dtype=str),
            "gene_name": pd.Series(dtype=str),
            "gene_id": pd.Series(dtype=str),
            "protein_id": pd.Series(dtype=str),
        }
    )
    captured: dict[str, pd.DataFrame] = {}

    def _capture_write(path, df):
        captured["df"] = df

    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_observations", return_value=hits) as ms_only,
        patch("hitlist.mappings.load_peptide_mappings", return_value=empty_mappings),
        patch("tsarina.cli_hits._write", side_effect=_capture_write),
        patch(
            "hitlist.aggregate.aggregate_per_pmhc",
            side_effect=lambda df: df[["peptide", "mhc_restriction"]].copy(),
        ),
    ):
        cli_hits.handle(args)

    kwargs = ms_only.call_args.kwargs
    assert kwargs.get("length_min") == 9
    assert kwargs.get("length_max") == 11

    out = captured["df"]
    got_lengths = {len(p) for p in out["peptide"]}
    assert got_lengths == {9, 11}


def test_cached_path_omitted_lengths_and_class_skips_bounds_pushdown(tmp_path):
    """When neither --lengths nor --mhc-class is set, the cached path
    must NOT push length_min / length_max.  Preserves "load everything"
    semantics for unqualified queries so off-window rows aren't silently
    dropped by an opinionated default."""
    from tsarina import cli_hits

    args = _base_cached_args(tmp_path, lengths=None, mhc_class=None)
    empty = pd.DataFrame({"peptide": pd.Series(dtype=str), "mhc_restriction": pd.Series(dtype=str)})
    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_observations", return_value=empty) as ms_only,
    ):
        cli_hits.handle(args)
    kwargs = ms_only.call_args.kwargs
    assert kwargs.get("length_min") is None
    assert kwargs.get("length_max") is None


def test_cached_path_class_I_derives_8_to_15_bounds(tmp_path):
    """--mhc-class I without --lengths populates length bounds as 8-15
    (wider than textbook 8-11 to catch phospho-extended class I ligands
    hitlist actually curates).  Explicit --lengths would override."""
    from tsarina import cli_hits

    args = _base_cached_args(tmp_path, lengths=None, mhc_class="I")
    empty = pd.DataFrame({"peptide": pd.Series(dtype=str), "mhc_restriction": pd.Series(dtype=str)})
    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_observations", return_value=empty) as ms_only,
    ):
        cli_hits.handle(args)
    kwargs = ms_only.call_args.kwargs
    assert kwargs.get("length_min") == 8
    assert kwargs.get("length_max") == 15


def test_cached_path_class_II_derives_12_to_45_bounds(tmp_path):
    """--mhc-class II without --lengths populates length bounds as 12-45
    (wider than textbook 12-25 to accommodate long class II tails)."""
    from tsarina import cli_hits

    args = _base_cached_args(tmp_path, lengths=None, mhc_class="II")
    empty = pd.DataFrame({"peptide": pd.Series(dtype=str), "mhc_restriction": pd.Series(dtype=str)})
    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_observations", return_value=empty) as ms_only,
    ):
        cli_hits.handle(args)
    kwargs = ms_only.call_args.kwargs
    assert kwargs.get("length_min") == 12
    assert kwargs.get("length_max") == 45


def test_cached_path_explicit_lengths_override_class_default(tmp_path):
    """Explicit --lengths must win over the --mhc-class-derived default:
    requesting --lengths 9,10 with --mhc-class I should push [9, 10],
    not the class I default [8, 15]."""
    from tsarina import cli_hits

    args = _base_cached_args(tmp_path, lengths=(9, 10), mhc_class="I")
    empty = pd.DataFrame({"peptide": pd.Series(dtype=str), "mhc_restriction": pd.Series(dtype=str)})
    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_observations", return_value=empty) as ms_only,
    ):
        cli_hits.handle(args)
    kwargs = ms_only.call_args.kwargs
    assert kwargs.get("length_min") == 9
    assert kwargs.get("length_max") == 10


def test_cached_path_reattaches_gene_columns_on_post_1_30_46_parquets(tmp_path):
    """hitlist 1.30.46 (#238) split gene_names / gene_ids / protein_ids out of
    observations.parquet — no-projection loads return them only if the caller
    already had them in columns=.  tsarina's fast path doesn't project, so it
    must re-attach those columns via peptide_mappings; otherwise per-peptide
    aggregation drops gene identifiers on freshly-built parquets."""
    from tsarina import cli_hits

    args = _base_cached_args(tmp_path, lengths=(9,))
    # Post-1.30.46 shape: no gene_names / gene_ids / protein_ids on the obs frame.
    hits = pd.DataFrame(
        {
            "peptide": ["AAAAAAAAA", "BBBBBBBBB"],
            "mhc_restriction": ["HLA-A*02:01"] * 2,
        }
    )
    mappings = pd.DataFrame(
        {
            "peptide": ["AAAAAAAAA", "AAAAAAAAAA", "BBBBBBBBB"],
            "gene_name": ["PRAME", "PRAME", "MAGEA4"],
            "gene_id": ["ENSG_PRAME", "ENSG_PRAME", "ENSG_MAGEA4"],
            "protein_id": ["ENSP_PRAME_1", "ENSP_PRAME_2", "ENSP_MAGEA4"],
        }
    )
    captured: dict[str, pd.DataFrame] = {}

    def _capture_write(path, df):
        captured["df"] = df

    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_observations", return_value=hits),
        patch("hitlist.mappings.load_peptide_mappings", return_value=mappings) as mappings_loader,
        patch("tsarina.cli_hits._write", side_effect=_capture_write),
        patch(
            "hitlist.aggregate.aggregate_per_pmhc",
            side_effect=lambda df: df[["peptide", "mhc_restriction"]].copy(),
        ),
    ):
        cli_hits.handle(args)

    mappings_loader.assert_called_once()
    out = captured["df"]
    assert "gene_names" in out.columns
    by_peptide = out.set_index("peptide")["gene_names"].to_dict()
    assert by_peptide["AAAAAAAAA"] == "PRAME"
    assert by_peptide["BBBBBBBBB"] == "MAGEA4"


def test_cached_path_reattaches_gene_columns_on_load_all_evidence_path(tmp_path):
    """Parity with ``test_cached_path_reattaches_gene_columns_on_post_1_30_46_parquets``
    for the ``include_binding_assays=True`` branch — both branches share
    the re-attach guard, so the binding-assays + observations union path
    must also pick up gene columns from peptide_mappings on freshly-built
    parquets."""
    from tsarina import cli_hits

    args = _base_cached_args(tmp_path, lengths=(9,), include_binding_assays=True)
    hits = pd.DataFrame(
        {
            "peptide": ["AAAAAAAAA"],
            "mhc_restriction": ["HLA-A*02:01"],
            "evidence_kind": ["ms"],
        }
    )
    mappings = pd.DataFrame(
        {
            "peptide": ["AAAAAAAAA"],
            "gene_name": ["PRAME"],
            "gene_id": ["ENSG_PRAME"],
            "protein_id": ["ENSP_PRAME"],
        }
    )
    captured: dict[str, pd.DataFrame] = {}

    def _capture_write(path, df):
        captured["df"] = df

    with (
        patch("tsarina.indexing.ensure_index_built"),
        patch("hitlist.observations.load_all_evidence", return_value=hits) as all_ev,
        patch("hitlist.mappings.load_peptide_mappings", return_value=mappings) as mappings_loader,
        patch("tsarina.cli_hits._write", side_effect=_capture_write),
        patch(
            "hitlist.aggregate.aggregate_per_pmhc",
            side_effect=lambda df: df[["peptide", "mhc_restriction"]].copy(),
        ),
    ):
        cli_hits.handle(args)

    all_ev.assert_called_once()
    mappings_loader.assert_called_once()
    out = captured["df"]
    assert "gene_names" in out.columns
    assert out.set_index("peptide")["gene_names"]["AAAAAAAAA"] == "PRAME"


def test_enumerate_gene_peptides_caches_proteome_index_across_calls():
    """Second invocation with the same (release, lengths) must reuse the
    cached ~8-15 GB Ensembl ProteomeIndex instead of rebuilding."""
    from tsarina import cli_hits

    class _StubProteomeIndex:
        proteins = {"P1": "MAGEAPEPTIDESEQX"}
        protein_meta = {"P1": {"gene_name": "TARGET", "gene_id": "ENSG"}}

    cli_hits._cached_proteome_index.cache_clear()
    try:
        with patch(
            "hitlist.proteome.ProteomeIndex.from_ensembl",
            return_value=_StubProteomeIndex(),
        ) as from_ensembl:
            df = cli_hits._enumerate_gene_peptides("TARGET", 112, (9,))
            cli_hits._enumerate_gene_peptides("TARGET", 112, (9,))
        from_ensembl.assert_called_once()
        assert "peptide" in df.columns
    finally:
        cli_hits._cached_proteome_index.cache_clear()
