from types import SimpleNamespace

import pandas as pd

from tsarina import (
    CTA_TIERS,
    CTA_by_tier,
    CTA_cancer_antigen_gene_ids,
    CTA_cancer_antigen_gene_names,
    CTA_evidence,
    CTA_excluded_gene_ids,
    CTA_gene_ids,
    CTA_gene_names,
    CTA_reproductive_antigen_gene_ids,
    CTA_reproductive_antigen_gene_names,
    CTA_selection,
    CTA_selection_coverage,
    CTA_target_candidate_gene_ids,
    CTA_target_candidate_gene_names,
)
from tsarina.vaccine_construct import VaccineConfig


def test_non_cta_and_unassigned_categories_are_not_relabelled_unsafe_or_selected():
    rows = CTA_evidence()
    other = rows.loc[rows.selection_category.eq("Other")]
    ca_only = rows.loc[rows.selection_category.eq("Cancer Antigen only")]
    ra_only = rows.loc[rows.selection_category.eq("Reproductive Antigen only")]
    # The rebuilt target_selection policy classifies every protein, so real data
    # may hold no unassigned rows. A null category must still never read as Unsafe or CTA.
    unassigned = rows.loc[rows.selection_category.isna()]
    assert not other.empty
    assert not ca_only.empty
    assert not ra_only.empty
    for subset in (other, ca_only, ra_only, unassigned):
        assert not subset.is_current_cta.any()
        assert set(subset.Ensembl_Gene_ID).isdisjoint(CTA_excluded_gene_ids())
        assert set(subset.Ensembl_Gene_ID).isdisjoint(CTA_gene_ids())


def test_rebuilt_target_selection_sets_and_cta_tiers_are_consistent():
    cta_ids = CTA_gene_ids()
    ca_ids = CTA_cancer_antigen_gene_ids()
    ra_ids = CTA_reproductive_antigen_gene_ids()
    cand_ids = CTA_target_candidate_gene_ids()

    assert cta_ids == ca_ids & ra_ids
    assert cand_ids == ca_ids | ra_ids
    assert (
        CTA_gene_names() == CTA_cancer_antigen_gene_names() & CTA_reproductive_antigen_gene_names()
    )
    assert CTA_target_candidate_gene_names() == (
        CTA_cancer_antigen_gene_names() | CTA_reproductive_antigen_gene_names()
    )

    tier_sets = [CTA_by_tier(tier, column="Ensembl_Gene_ID") for tier in CTA_TIERS]
    assert all(tier_sets)
    assert set().union(*tier_sets) == cta_ids
    assert sum(len(s) for s in tier_sets) == len(cta_ids)
    assert "PRAME" in CTA_by_tier("core_with_extended_exposure")
    assert "MAGEA4" in CTA_by_tier("core")
    assert "MAGEC2" in CTA_by_tier("core")


def test_nonhuman_cta_selection_delegates_to_species_catalog():
    required_cols = {
        "species",
        "entity_id",
        "entity_kind",
        "label",
        "gene_ids",
        "category",
        "subtype",
        "policy_id",
        "reference_id",
        "source_snapshot_id",
    }
    for sp, expected_species in (("dog", "dog"), ("mouse", "mouse")):
        df = CTA_selection(species=sp)
        assert not df.empty, sp
        assert required_cols <= set(df.columns), sp
        assert set(df["species"].dropna().unique()) == {expected_species}
        assert df["policy_id"].notna().all(), sp
        assert df["reference_id"].notna().all(), sp

    mouse_cov = CTA_selection_coverage(species="mouse")
    assert mouse_cov["selection_status"] == "evaluated"
    assert "ensembl115" in mouse_cov.get("registered_references", [])

    dog_cov = CTA_selection_coverage(species="dog")
    assert dog_cov["selection_status"] == "pending_reconciled_panels"

    human_cov = CTA_selection_coverage(species="human")
    assert isinstance(human_cov, dict) and human_cov


def test_oncoref_146_curation_reclassifications():
    ev = CTA_evidence().set_index("Symbol")
    assert ev.loc["SP7", "selection_category"] == "Other"
    assert ev.loc["SP7", "selection_subtype"] == "Unscreened tissue"
    assert "SP7" not in CTA_cancer_antigen_gene_names()
    for sym in ("RGSL1", "OOSP4B"):
        assert ev.loc[sym, "selection_category"] == "Unsafe"
        assert sym not in CTA_gene_names()


def test_current_partition_retains_every_nonselected_category_in_background(monkeypatch):
    import tsarina.partition as partition

    frame = pd.DataFrame(
        {
            "Ensembl_Gene_ID": ["selected", "unsafe", "other", "unknown"],
            "Symbol": ["CTA", "UNSAFE", "OTHER", "UNKNOWN"],
            "selection_category": ["CTA", "Unsafe", "Other", pd.NA],
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
    assert result.non_cta == {"unsafe", "other", "unknown", "outside"}
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


def test_other_personalized_target_reports_current_reason_and_provenance():
    import pytest

    from tsarina import personalized_targets

    other = CTA_evidence().query(
        "selection_category == 'Other' and Symbol.str.fullmatch('[A-Z0-9]+')"
    )
    symbol = other.Symbol.sort_values().iloc[0]
    with pytest.warns(UserWarning, match=f"{symbol}.*current OncoRef category: Other"):
        result = personalized_targets(
            hla_alleles=[],
            cta_expression={symbol: None},
            score_presentation=False,
            skip_ms_evidence=True,
            ensembl_release=112,
        )
    assert result.empty
    assert result.attrs["ensembl_release"] == 112
    assert result.attrs["cta_selection_reference_ids"]


def test_resolve_cta_alternate_haplotype_ids_excludes_primary_unnamed_and_cross_chrom():
    """tsarina#187: alternate-haplotype resolution matches only HSCHR* contigs
    backed by OncoRef alternate IDs or same-chromosome non-empty symbols."""
    from tsarina.partition import resolve_cta_alternate_haplotype_ids

    genes = [
        SimpleNamespace(gene_id="ENSG_PRAME", gene_name="PRAME", contig="22"),
        SimpleNamespace(gene_id="ENSG_PRAME_ALT", gene_name="PRAME", contig="HSCHR22_1_CTG3"),
        SimpleNamespace(gene_id="ENSG_PRAME_CHR1_ALT", gene_name="PRAME", contig="HSCHR1_1_CTG1"),
        SimpleNamespace(gene_id="ENSG_PRAME_PRIMARY", gene_name="PRAME", contig="1"),
        SimpleNamespace(gene_id="ENSG_UNNAMED_CTA", gene_name="", contig="X"),
        SimpleNamespace(gene_id="ENSG_UNNAMED_ALT", gene_name="", contig="HSCHRX_1_CTG1"),
        SimpleNamespace(gene_id="ENSG_CAT_ALT", gene_name="DIFFERENT", contig="HSCHR22_2_CTG1"),
    ]
    resolved_ids, aliases = resolve_cta_alternate_haplotype_ids(
        genes,
        {"ENSG_PRAME", "ENSG_UNNAMED_CTA"},
        catalog_alternates={"ENSG_CAT_ALT": {"ENSG_PRAME"}},
    )
    assert resolved_ids == {"ENSG_PRAME_ALT", "ENSG_CAT_ALT"}
    assert "ENSG_PRAME_PRIMARY" not in resolved_ids
    assert "ENSG_PRAME_CHR1_ALT" not in resolved_ids
    assert "ENSG_UNNAMED_ALT" not in resolved_ids
    by_id = aliases.set_index("gene_id")
    assert by_id.loc["ENSG_PRAME_ALT", "resolution_basis"] == "hschr_same_symbol"
    assert by_id.loc["ENSG_CAT_ALT", "resolution_basis"] == "oncoref_alternate_id"


def test_partition_excludes_hschr_alternate_haplotypes_from_non_cta(monkeypatch):
    """tsarina#187: CTA_partition_gene_ids excludes audited HSCHR* alternate
    haplotypes from non_cta while keeping partition.cta equal to CTA_gene_ids()."""
    import tsarina.partition as partition

    frame = pd.DataFrame(
        {
            "Ensembl_Gene_ID": ["ENSG_PRAME", "ENSG_OTHER"],
            "Symbol": ["PRAME", "OTHER"],
            "selection_category": ["CTA", "Other"],
        }
    )
    genes = [
        SimpleNamespace(
            gene_id="ENSG_PRAME", gene_name="PRAME", contig="22", biotype="protein_coding"
        ),
        SimpleNamespace(
            gene_id="ENSG_PRAME_ALT",
            gene_name="PRAME",
            contig="HSCHR22_1_CTG3",
            biotype="protein_coding",
        ),
        SimpleNamespace(
            gene_id="ENSG_PRAME_PRIMARY_BG",
            gene_name="PRAME",
            contig="1",
            biotype="protein_coding",
        ),
        SimpleNamespace(
            gene_id="ENSG_OTHER", gene_name="OTHER", contig="2", biotype="protein_coding"
        ),
    ]
    monkeypatch.setattr(partition, "CTA_evidence", lambda: frame)
    monkeypatch.setattr(partition, "CTA_gene_ids", lambda: {"ENSG_PRAME"})
    monkeypatch.setattr(partition, "_oncoref_alternate_id_to_canonical", lambda: {})
    monkeypatch.setattr(
        "pyensembl.EnsemblRelease", lambda release: SimpleNamespace(genes=lambda: genes)
    )

    result = partition.CTA_partition_gene_ids(ensembl_release=112)
    aliases = partition.CTA_partition_alternate_haplotype_aliases(ensembl_release=112)
    assert result.cta == {"ENSG_PRAME"}
    assert result.cta_alternate_haplotypes == {"ENSG_PRAME_ALT"}
    assert result.non_cta == {"ENSG_PRAME_PRIMARY_BG", "ENSG_OTHER"}
    assert list(aliases["gene_id"]) == ["ENSG_PRAME_ALT"]


def test_ensembl_112_partition_resolves_prame_and_magea6_alternate_haplotypes():
    """tsarina#187: real Ensembl 112 partition resolves PRAME and MAGEA6 HSCHR*
    alternate haplotypes so they do not self-veto their primary CTAs."""
    import pytest

    try:
        from pyensembl import EnsemblRelease

        EnsemblRelease(112).gene_by_id("ENSG00000185686")
    except Exception as exc:
        pytest.skip(f"pyensembl release 112 data unavailable: {exc}")

    from tsarina.partition import (
        CTA_partition_alternate_haplotype_aliases,
        CTA_partition_gene_ids,
    )
    from tsarina.peptides import cta_exclusive_peptides

    partition = CTA_partition_gene_ids(ensembl_release=112)
    aliases = CTA_partition_alternate_haplotype_aliases(ensembl_release=112)

    assert partition.cta == CTA_gene_ids()
    assert partition.cta.isdisjoint(partition.non_cta)
    assert partition.cta.isdisjoint(partition.cta_alternate_haplotypes)
    assert partition.cta_alternate_haplotypes.isdisjoint(partition.non_cta)

    # PRAME (ENSG00000185686) and MAGEA6 (ENSG00000197172) are current CTAs whose
    # identical HSCHR* annotations must not appear in non_cta.
    assert "ENSG00000185686" in partition.cta
    assert "ENSG00000275013" in partition.cta_alternate_haplotypes
    assert "ENSG00000275013" not in partition.non_cta
    assert "ENSG00000197172" in partition.cta
    assert "ENSG00000292195" in partition.cta_alternate_haplotypes
    assert "ENSG00000292195" not in partition.non_cta

    # Unnamed HSCHR loci must never be swept into the alias table.
    assert not aliases.empty
    assert (aliases["gene_name"].astype(str).str.strip() != "").all()
    assert aliases["contig"].str.startswith("HSCHR").all()

    # PRAME now retains CTA-exclusive 9-mers rather than being self-vetoed to zero.
    prame_peptides = cta_exclusive_peptides(
        gene_names=["PRAME"],
        lengths=(9,),
        ensembl_release=112,
    )
    assert not prame_peptides.empty
    assert "SLLQHLIGL" in set(prame_peptides["peptide"])
