from types import SimpleNamespace

import pandas as pd

from tsarina import CTA_evidence, CTA_excluded_gene_ids, CTA_gene_ids, CTA_gene_names, CTA_selection
from tsarina.vaccine_construct import VaccineConfig


def test_neither_and_unassigned_are_not_relabelled_unsafe_or_selected():
    rows = CTA_evidence()
    neither = rows.loc[rows.selection_category.eq("Neither")]
    unassigned = rows.loc[rows.selection_category.isna()]
    assert not neither.empty and not unassigned.empty
    for subset in (neither, unassigned):
        assert not subset.is_current_cta.any()
        assert set(subset.Ensembl_Gene_ID).isdisjoint(CTA_excluded_gene_ids())
        assert set(subset.Ensembl_Gene_ID).isdisjoint(CTA_gene_ids())


def test_current_partition_retains_every_nonselected_category_in_background(monkeypatch):
    import tsarina.partition as partition

    frame = pd.DataFrame(
        {
            "Ensembl_Gene_ID": ["selected", "unsafe", "neither", "unknown"],
            "Symbol": ["CTA", "UNSAFE", "NEITHER", "UNKNOWN"],
            "selection_category": ["CTA", "Unsafe", "Neither", pd.NA],
        }
    )
    genes = [
        SimpleNamespace(gene_id=g, gene_name=g, biotype="protein_coding")
        for g in [*frame.Ensembl_Gene_ID, "outside"]
    ]
    monkeypatch.setattr(partition, "CTA_evidence", lambda: frame)
    monkeypatch.setattr(partition, "CTA_gene_ids", lambda: {"selected"})
    monkeypatch.setattr(
        "pyensembl.EnsemblRelease", lambda release: SimpleNamespace(genes=lambda: genes)
    )
    result = partition.CTA_partition_gene_ids()
    assert result.cta == {"selected"}
    assert result.cta_never_expressed == set()
    assert result.non_cta == {"unsafe", "neither", "unknown", "outside"}
    names = partition.CTA_partition_gene_names()
    assert "selected" not in names.non_cta  # Sequence and selection symbols can differ.
    assert names.cta == {"CTA"}


def test_human_vaccine_and_cli_default_to_current_policy():
    import argparse

    from tsarina.cli_vaccine import build_parser

    assert VaccineConfig().definition == "current"
    assert VaccineConfig(definition="strict").definition == "strict"
    assert VaccineConfig(definition="loose").definition == "loose"
    assert VaccineConfig.for_canine("test").definition == "strict"
    parser = argparse.ArgumentParser()
    build_parser(parser.add_subparsers(dest="command"))
    args = parser.parse_args(["vaccine", "--output-dir", "/tmp/test-current-policy"])
    assert args.cta_definition is None  # Species-aware VaccineConfig chooses current for human.


def test_current_selection_returns_detached_upstream_decisions():
    selected = CTA_selection(category="CTA")
    selected.loc[:, "category"] = "Unsafe"
    assert CTA_selection(category="CTA").category.eq("CTA").all()


def test_default_panel_source_candidates_preserve_every_current_member():
    from tsarina.spanning import _cta_member_symbols, _resolve_ctas

    labels = _resolve_ctas(ctas=None, cta_count=None, cta_rank_by="Symbol", selection_allowlist=[])
    assert {gene for label in labels for gene in _cta_member_symbols(label)} == CTA_gene_names()


def test_unassigned_personalized_target_reports_current_reason_and_provenance():
    import pytest

    from tsarina import personalized_targets

    with pytest.warns(UserWarning, match="GAGE10.*unassigned"):
        result = personalized_targets(
            hla_alleles=[],
            cta_expression={"GAGE10": None},
            score_presentation=False,
            skip_ms_evidence=True,
            ensembl_release=112,
        )
    assert result.empty
    assert result.attrs["ensembl_release"] == 112
    assert result.attrs["cta_selection_reference_ids"]
