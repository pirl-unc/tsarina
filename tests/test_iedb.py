import pandas as pd

from tsarina.iedb import _species_kwargs, scan_public_ms


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


# ── _species_kwargs ──────────────────────────────────────────────────────
# Regression guard: the Homo-sapiens branch must stay on human_only=True so
# the host/species fallback kicks in for unparseable MHC restrictions (e.g.
# bare "HLA class I").  If mhcgnomes ever changes its canonical species
# form, these tests catch the silent drop to the strict path.


def test_species_kwargs_none_disables_filter():
    assert _species_kwargs(None) == {"human_only": False}


def test_species_kwargs_homo_sapiens_uses_human_only_fallback():
    assert _species_kwargs("Homo sapiens") == {"human_only": True}


def test_species_kwargs_human_alias_normalizes_to_human_only():
    # mhcgnomes should normalize "human" → "Homo sapiens".
    assert _species_kwargs("human") == {"human_only": True}


def test_species_kwargs_non_human_uses_strict_path():
    assert _species_kwargs("Mus musculus") == {
        "mhc_species": "Mus musculus",
        "human_only": False,
    }
