import warnings

import pandas as pd

from tsarina.iedb import scan_public_ms


def test_scan_no_sources_returns_empty():
    df = scan_public_ms(peptides={"SLYNTVATL"}, iedb_path=None, cedar_path=None)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 0


def test_scan_nonexistent_path_returns_empty():
    df = scan_public_ms(
        peptides={"SLYNTVATL"},
        iedb_path="/nonexistent/path.csv",
    )
    assert len(df) == 0


def test_scan_empty_peptides_returns_empty():
    df = scan_public_ms(peptides=set(), iedb_path=None)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 0


def test_scan_does_not_emit_deprecation_warnings():
    """After the hitlist#72 / tsarina follow-up, scan_public_ms passes
    mhc_species= directly to hitlist.  If we ever regress back to the old
    human_only= kwarg, hitlist 1.9.0+ emits a DeprecationWarning — fail here."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        scan_public_ms(peptides={"X"}, iedb_path=None, cedar_path=None)
        scan_public_ms(peptides={"X"}, iedb_path=None, cedar_path=None, mhc_species=None)
        scan_public_ms(peptides={"X"}, iedb_path=None, cedar_path=None, mhc_species="Mus musculus")
