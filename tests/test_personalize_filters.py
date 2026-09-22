"""Selection diagnostics and clinical-target overrides must compose with gates."""

import warnings

import pandas as pd
import pytest

from tsarina.personalize import personalize


@pytest.fixture
def cta_selection(monkeypatch):
    generated = {"strict": [], "flagged": []}
    monkeypatch.setattr("tsarina.gene_sets.CTA_gene_names", lambda: {"STRICT"})
    monkeypatch.setattr(
        "tsarina.gene_sets.CTA_clinical_target_gene_names", lambda: {"STRICT", "FLAGGED"}
    )
    monkeypatch.setattr("tsarina.gene_sets.CTA_by_axes", lambda **kw: {"STRICT"})
    monkeypatch.setattr("tsarina.gene_sets.CTA_excluded_gene_names", lambda: {"EXCLUDED"})
    monkeypatch.setattr("tsarina.gene_sets.CTA_never_expressed_gene_names", lambda: {"SILENT"})
    monkeypatch.setattr(
        "tsarina.gene_sets.CTA_unfiltered_gene_names",
        lambda: {"STRICT", "FLAGGED", "EXCLUDED", "SILENT"},
    )
    monkeypatch.setattr(
        "tsarina.personalize._cta_flag_rationale", lambda genes: dict.fromkeys(genes, "kidney IHC")
    )

    def peptides(genes, category):
        generated[category].extend(genes)
        return pd.DataFrame(
            [(gene, "ENSG_" + gene, "PEPTIDEAA", 9) for gene in genes],
            columns=["gene_name", "gene_id", "peptide", "length"],
        )

    monkeypatch.setattr(
        "tsarina.peptides.cta_exclusive_peptides", lambda **kw: peptides(["STRICT"], "strict")
    )
    monkeypatch.setattr(
        "tsarina.personalize._cta_flagged_gene_peptides",
        lambda genes, **kw: peptides(genes, "flagged"),
    )
    return generated


def select(expression, **kwargs):
    return personalize(
        hla_alleles=["HLA-A*02:01"],
        cta_expression=expression,
        score_presentation=False,
        skip_ms_evidence=True,
        drop_weak_tier=False,
        show_progress=False,
        proteoform_rollup=False,
        **kwargs,
    )


def test_confidence_rejection_cannot_turn_strict_cta_into_flagged(cta_selection, monkeypatch):
    monkeypatch.setattr("tsarina.gene_sets.CTA_by_axes", lambda **kw: set())
    with pytest.warns(UserWarning, match="STRICT.*restriction.confidence"):
        result = select({"STRICT": 20.0})
    assert result.empty
    assert cta_selection == {"strict": [], "flagged": []}


@pytest.mark.parametrize("gene", ["STRICT", "FLAGGED"])
@pytest.mark.parametrize("measured", [True, False])
def test_mtec_gate_applies_to_both_categories(cta_selection, tmp_path, gene, measured):
    path = tmp_path / "mtec.tsv"
    path.write_text("gene_symbol\tsample\n" + (f"{gene}\t5\n" if measured else "OTHER\t0\n"))
    with pytest.warns(UserWarning, match=gene + ".*mTEC"):
        result = select({gene: 20.0}, mtec_matrix_path=path)
    assert result.empty
    assert cta_selection == {"strict": [], "flagged": []}


def test_mtec_pass_preserves_both_categories_and_flag_caveat(cta_selection, tmp_path):
    path = tmp_path / "mtec.tsv"
    path.write_text("gene_symbol\tsample\nSTRICT\t0\nFLAGGED\t1\n")
    result = select({"STRICT": 20.0, "FLAGGED": None}, mtec_matrix_path=path)
    assert set(zip(result.source, result.category)) == {
        ("STRICT", "cta"),
        ("FLAGGED", "cta_flagged"),
    }
    flagged = result[result.source == "FLAGGED"].iloc[0]
    assert "kidney IHC" in flagged.flag_reason
    assert "not screened for peptide overlap" in flagged.flag_reason
    assert pd.isna(flagged.source_tpm)


def test_upstream_exclusion_low_expression_and_unknown_have_distinct_warnings(cta_selection):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = select({"EXCLUDED": 20.0, "SILENT": 20.0, "TYPO": 20.0})
    messages = [str(w.message) for w in caught]
    assert any("EXCLUDED" in m and "kidney IHC" in m for m in messages)
    assert any("SILENT" in m and "expression floor" in m for m in messages)
    assert any("TYPO" in m and "not a recognized CTA symbol" in m for m in messages)
    assert result.empty
    assert cta_selection == {"strict": [], "flagged": []}


@pytest.mark.parametrize("gene", ["STRICT", "FLAGGED"])
def test_low_measured_tpm_warns_for_both_categories(cta_selection, gene):
    with pytest.warns(UserWarning, match=gene + ".*TPM.*2"):
        result = select({gene: 0.0})
    assert result.empty
    assert cta_selection == {"strict": [], "flagged": []}


def test_disabling_confidence_retains_strict_exclusivity_path(cta_selection, monkeypatch):
    monkeypatch.setattr("tsarina.gene_sets.CTA_by_axes", lambda **kw: set())
    result = select({"STRICT": 20.0}, min_restriction_confidence=None)
    assert list(result.category) == ["cta"]
    assert cta_selection == {"strict": ["STRICT"], "flagged": []}
