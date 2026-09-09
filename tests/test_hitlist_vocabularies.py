"""Tripwires for the hitlist vocabularies tsarina keys on.

tsarina restates two upstream vocabularies as literals — the CLI's
``--min-resolution`` choices, so ``--help`` does not pay for importing
``hitlist.curation``, and the provenance values that count as sample-narrowed
evidence.  Restating them is only safe while something fails when hitlist's own
vocabulary moves; that is what these tests are.
"""

from hitlist import MHC_ALLELE_PROVENANCE_VALUES
from hitlist.curation import ALLELE_RESOLUTION_ORDER

from tsarina.cli_hits import _SUPPORTED_RESOLUTIONS
from tsarina.spanning import _SAMPLE_NARROWED_PROVENANCES

#: Provenance values that name no donor, so they carry no sample-narrowed
#: allele bag: ``exact`` needs none, and the other two are study-wide pools.
_NOT_SAMPLE_NARROWED = {"exact", "pmid_class_pool", "unmatched"}


def test_supported_resolutions_are_hitlist_order_most_specific_first():
    """``--min-resolution`` ranks by hitlist's order, so it must mirror it."""
    most_specific = tuple(ALLELE_RESOLUTION_ORDER[: len(_SUPPORTED_RESOLUTIONS)])
    assert most_specific == _SUPPORTED_RESOLUTIONS


def test_supported_resolutions_omit_only_the_catch_all_level():
    """Every level a user could sensibly ask for is offered.

    ``unresolved`` is the one deliberate omission: it ranks last, so asking
    for it keeps every row and filters nothing.
    """
    assert tuple(ALLELE_RESOLUTION_ORDER[len(_SUPPORTED_RESOLUTIONS) :]) == ("unresolved",)


def test_sample_narrowed_provenances_partition_hitlist_vocabulary():
    """A new upstream provenance value has to be classified, not ignored."""
    assert _SAMPLE_NARROWED_PROVENANCES.issubset(MHC_ALLELE_PROVENANCE_VALUES)
    assert set(MHC_ALLELE_PROVENANCE_VALUES) - _SAMPLE_NARROWED_PROVENANCES == (
        _NOT_SAMPLE_NARROWED
    )


def test_coding_gene_biotypes_track_hitlist():
    """tsarina's coding universe is hitlist's, not a copy of it.

    Four sites used to say ``protein_coding`` independently while hitlist's
    ``proteome_kmer_set`` had widened to include germline IG/TR segments, so
    ``human_exclusive_viral_peptides`` subtracted a larger self set than the
    CTA / non-CTA partition covered.
    """
    from hitlist.proteome import ENSEMBL_CODING_GENE_BIOTYPES

    from tsarina.gene_sets import CODING_GENE_BIOTYPES

    assert set(ENSEMBL_CODING_GENE_BIOTYPES) == CODING_GENE_BIOTYPES


def test_coding_universe_includes_germline_ig_and_tr_segments():
    """Named explicitly: these are translated and presented, so they are self.

    Ensembl gives an IG_V gene's transcripts the biotype ``IG_V_gene`` rather
    than ``protein_coding``, so a transcript filter keyed on ``protein_coding``
    silently excludes them even when the gene is in the universe.
    """
    from tsarina.gene_sets import CODING_GENE_BIOTYPES, is_coding_gene, is_coding_transcript

    for biotype in (
        "IG_V_gene",
        "IG_D_gene",
        "IG_J_gene",
        "IG_C_gene",
        "TR_V_gene",
        "TR_D_gene",
        "TR_J_gene",
        "TR_C_gene",
        "protein_coding",
    ):
        assert biotype in CODING_GENE_BIOTYPES

    class _Feature:
        def __init__(self, biotype: str) -> None:
            self.biotype = biotype

    assert is_coding_gene(_Feature("IG_V_gene"))
    assert is_coding_transcript(_Feature("TR_V_gene"))
    # Non-coding transcripts of a coding gene stay out.
    assert not is_coding_transcript(_Feature("retained_intron"))
    assert not is_coding_transcript(_Feature("processed_transcript"))
    assert not is_coding_gene(_Feature("lncRNA"))
    # A feature with no biotype at all must not slip through.
    assert not is_coding_transcript(object())
