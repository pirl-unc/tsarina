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


def test_scan_public_ms_filters_to_positive_ms_by_default(monkeypatch):
    """tsarina#199: scan_public_ms applies Hitlist positive-MS admission by default."""
    raw = pd.DataFrame(
        [
            {
                "peptide": "SLYNTVATL",
                "assay_method": "cellular MHC/mass spectrometry",
                "qualitative_measurement": "Positive",
                "source": "iedb",
            },
            {
                "peptide": "GVYDGREHTV",
                "assay_method": "x-ray crystallography",
                "qualitative_measurement": "Positive",
                "source": "iedb",
            },
            {
                "peptide": "ACDEFGHIK",
                "assay_method": "mass spectrometry",
                "qualitative_measurement": "unknown",
                "source": "iedb",
            },
            {
                "peptide": "KLMNPQRST",
                "assay_method": "purified MHC/direct/fluorescence",
                "response_measured": "dissociation constant KD",
                "qualitative_measurement": "Positive",
                "source": "iedb",
            },
        ]
    )
    monkeypatch.setattr("tsarina.iedb._hitlist_scan", lambda **kwargs: raw.copy())

    filtered = scan_public_ms(peptides=set(raw["peptide"]), iedb_path="/tmp/fake.csv")
    assert list(filtered["peptide"]) == ["SLYNTVATL"]

    unfiltered = scan_public_ms(
        peptides=set(raw["peptide"]),
        iedb_path="/tmp/fake.csv",
        positive_ms_only=False,
    )
    assert set(unfiltered["peptide"]) == set(raw["peptide"])
