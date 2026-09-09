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
