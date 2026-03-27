import pandas as pd

from hitlist.targets import target_summary


def test_target_summary_on_empty():
    df = pd.DataFrame(columns=["peptide", "length", "category", "source", "source_detail"])
    summary = target_summary(df)
    assert isinstance(summary, pd.DataFrame)
    assert len(summary) == 0


def test_target_summary_categories():
    df = pd.DataFrame(
        {
            "peptide": ["AAAA", "BBBB", "CCCC", "DDDD"],
            "length": [4, 4, 4, 4],
            "category": ["cta", "cta", "viral", "mutant"],
            "source": ["MAGEA4", "MAGEA4", "HPV-16", "KRAS G12D"],
            "source_detail": ["ENSG1", "ENSG1", "P1", "G12D"],
            "has_ms_evidence": [True, False, True, False],
            "ms_in_cancer": [True, False, False, False],
        }
    )
    summary = target_summary(df)
    assert len(summary) == 3
    assert set(summary["category"]) == {"cta", "viral", "mutant"}
