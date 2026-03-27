import pandas as pd

from ctabase.iedb import scan_public_ms


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
