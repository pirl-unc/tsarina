"""Tests for OncoRef 1.52+ CTA safety caveats and 1-aa near-self screening (#259)."""

from __future__ import annotations

import sys

import pandas as pd
import pytest

import tsarina
from tsarina.cli import main
from tsarina.cli_spanning import format_panel_table
from tsarina.near_self import (
    annotate_cta_safety_caveats,
    annotate_near_self_1mm,
    screen_near_self_1mm,
)
from tsarina.personalize import format_table, personalized_targets
from tsarina.spanning import spanning_pmhc_set


def test_top_level_exports_near_self_and_caveat_helpers():
    assert tsarina.annotate_cta_safety_caveats is annotate_cta_safety_caveats
    assert tsarina.annotate_near_self_1mm is annotate_near_self_1mm
    assert tsarina.screen_near_self_1mm is screen_near_self_1mm


def test_annotate_cta_safety_caveats_covers_ctag2_nuclear_and_high_near_self():
    df = pd.DataFrame(
        [
            {"source": "CTAG2", "peptide": "SLLMWITQC", "category": "cta"},
            {"source": "CTAG2", "peptide": "GARGPESRL", "category": "cta"},
            {"source": "PLAC1", "peptide": "VLLLVLLLV", "category": "cta"},
            {"source": "DUX4", "peptide": "ALQALQALQ", "category": "cta"},
            {"source": "NUTM1", "peptide": "KLLKLLKLL", "category": "cta"},
            {"source": "MAGEA6", "peptide": "KVADLVGFL", "category": "cta"},
            {"source": "MAGEA4", "peptide": "GVYDGREHTV", "category": "cta"},
            {"source": "PRAME", "peptide": "SLLQHLIGL", "category": "cta"},
            {"source": "KRAS G12D", "peptide": "VVVGADGVGK", "category": "mutant"},
        ]
    )

    out = annotate_cta_safety_caveats(df)

    # CTAG2 SLLMWITQC is clinically cleared with no unverified-epitope caveat.
    assert bool(out.loc[0, "clinical_cleared_epitope"]) is True
    assert out.loc[0, "cta_safety_note"] == ""
    assert out.loc[0, "unsafe_near_self_1mm_pct"] == pytest.approx(10.47, abs=0.5)

    # Non-SLLMWITQC CTAG2 peptide is flagged as outside epitope-scoped clearance.
    assert bool(out.loc[1, "clinical_cleared_epitope"]) is False
    assert "SLLMWITQC" in out.loc[1, "cta_safety_note"]
    assert "unverified epitope" in out.loc[1, "cta_safety_note"]

    # PLAC1, DUX4, NUTM1 surface OncoRef nuclear_rna_caution.
    assert "snrnaseq_intronic_spillover_cleared" in out.loc[2, "cta_safety_note"]
    assert "snrnaseq_intronic_spillover_cleared" in out.loc[3, "cta_safety_note"]
    assert "bidirectional_promoter_prompt_nop10_cleared" in out.loc[4, "cta_safety_note"]

    # MAGEA6 has >80% 1-mm near-self overlap.
    assert float(out.loc[5, "unsafe_near_self_1mm_pct"]) >= 80.0
    assert "high 1-mm near-self overlap" in out.loc[5, "cta_safety_note"]

    # MAGEA4 GVYDGREHTV and PRAME SLLQHLIGL are clinically cleared epitopes.
    assert bool(out.loc[6, "clinical_cleared_epitope"]) is True
    assert out.loc[6, "cta_safety_note"] == ""
    assert bool(out.loc[7, "clinical_cleared_epitope"]) is True
    assert out.loc[7, "cta_safety_note"] == ""

    # Non-CTA rows get clean defaults.
    assert bool(out.loc[8, "clinical_cleared_epitope"]) is False
    assert out.loc[8, "cta_safety_note"] == ""
    assert pd.isna(out.loc[8, "unsafe_near_self_1mm_pct"])


def test_screen_near_self_1mm_exact_hamming_semantics():
    bg = {
        "MAGEA12": "XXXXKMAELVHFLXXXX",  # 1-mm from KVAELVHFL (V->M)
        "TITIN": "XXXXESDPIVAQYXXXX",  # 2-mm from EVDPIGHLY (not 1-mm)
        "INDEL_BG": "XXXXSLLMWITQXXXX",  # 1-aa shorter deletion from SLLMWITQC (not Hamming 1-mm)
        "EXACT_BG": "XXXXSLLMWITQCXXXX",  # exact 0-mm match (excluded from 1-mm count)
        "ALT_BG": "XXXXKVAELVHFIXXXX",  # second 1-mm from KVAELVHFL (L->I)
    }

    res = screen_near_self_1mm(
        ["KVAELVHFL", "SLLMWITQC", "EVDPIGHLY"],
        background_proteins=bg,
    ).set_index("peptide")

    assert res.loc["KVAELVHFL", "near_self_1mm_count"] == 2
    assert res.loc["KVAELVHFL", "near_self_1mm_genes"] == "ALT_BG;MAGEA12"
    assert res.loc["KVAELVHFL", "near_self_1mm_peptides"] == "KMAELVHFL;KVAELVHFI"

    # Exact match and indel are NOT counted as 1-mm substitutions.
    assert res.loc["SLLMWITQC", "near_self_1mm_count"] == 0
    assert res.loc["SLLMWITQC", "near_self_1mm_genes"] == ""
    assert res.loc["EVDPIGHLY", "near_self_1mm_count"] == 0
    assert "near_self_contract_id" in res.attrs


def test_annotate_near_self_1mm_and_max_matches_gate():
    df = pd.DataFrame(
        [
            {"peptide": "KVAELVHFL", "source": "MAGEA3"},
            {"peptide": "SLLMWITQC", "source": "CTAG2"},
        ]
    )
    bg = {"MAGEA12": "AAAKMAELVHFLAAA"}

    annotated = annotate_near_self_1mm(df, background_proteins=bg)
    assert annotated["near_self_1mm_count"].tolist() == [1, 0]
    assert annotated["near_self_1mm_genes"].tolist() == ["MAGEA12", ""]

    gated = annotate_near_self_1mm(df, background_proteins=bg, max_matches=0)
    assert gated["peptide"].tolist() == ["SLLMWITQC"]

    with pytest.raises(ValueError, match="max_matches"):
        annotate_near_self_1mm(df, background_proteins=bg, max_matches=-1)
    with pytest.raises(ValueError, match="near_self_background"):
        screen_near_self_1mm(["SLLMWITQC"], background="Invalid", background_proteins=bg)


def test_personalized_targets_surfaces_cta_caveats_and_optional_1mm_screen():
    custom = pd.DataFrame(
        [
            {"peptide": "SLLMWITQC", "category": "cta", "source": "CTAG2"},
            {"peptide": "GARGPESRL", "category": "cta", "source": "CTAG2"},
            {"peptide": "KVADLVGFL", "category": "cta", "source": "MAGEA6"},
        ]
    )
    bg = {"MAGEA12": "ZZZKVADLVGFMZZZ"}  # 1-mm from KVADLVGFL

    out = personalized_targets(
        hla_alleles=["HLA-A*02:01"],
        custom_targets=custom,
        score_presentation=False,
        skip_ms_evidence=True,
        drop_weak_tier=False,
        show_progress=False,
        screen_near_self_1mm=True,
        near_self_background_proteins=bg,
    )
    by_pep = out.set_index("peptide")
    assert bool(by_pep.loc["SLLMWITQC", "clinical_cleared_epitope"]) is True
    assert by_pep.loc["SLLMWITQC", "cta_safety_note"] == ""
    assert "SLLMWITQC" in by_pep.loc["GARGPESRL", "cta_safety_note"]
    assert "high 1-mm near-self overlap" in by_pep.loc["KVADLVGFL", "cta_safety_note"]
    assert by_pep.loc["KVADLVGFL", "near_self_1mm_count"] == 1
    assert by_pep.loc["KVADLVGFL", "near_self_1mm_genes"] == "MAGEA12"
    assert by_pep.loc["SLLMWITQC", "near_self_1mm_count"] == 0
    assert out.attrs["near_self_1mm_screened"] is True

    # format_table shows compact 1mm_self column + concise footnotes.
    rendered = format_table(out)
    assert "1mm_self" in rendered
    assert "OncoRef CTA safety caveats" in rendered
    assert "1-aa near-self matches" in rendered

    # Strict max_near_self_1mm_matches=0 drops KVADLVGFL.
    strict = personalized_targets(
        hla_alleles=["HLA-A*02:01"],
        custom_targets=custom,
        score_presentation=False,
        skip_ms_evidence=True,
        drop_weak_tier=False,
        show_progress=False,
        max_near_self_1mm_matches=0,
        near_self_background_proteins=bg,
    )
    assert "KVADLVGFL" not in set(strict["peptide"])
    assert set(strict["peptide"]) == {"SLLMWITQC", "GARGPESRL"}


def test_spanning_long_and_panel_table_surface_cta_caveats_and_1mm_screen(monkeypatch):
    def fake_peptides(**kwargs):
        return pd.DataFrame(
            [
                {
                    "gene_name": "CTAG2",
                    "gene_id": "ENSG00000126890",
                    "peptide": "SLLMWITQC",
                    "length": 9,
                },
                {
                    "gene_name": "CTAG2",
                    "gene_id": "ENSG00000126890",
                    "peptide": "GARGPESRL",
                    "length": 9,
                },
                {
                    "gene_name": "MAGEA6",
                    "gene_id": "ENSG00000197172",
                    "peptide": "KVADLVGFL",
                    "length": 9,
                },
            ]
        )

    def fake_score(peptides, alleles, predictor="mhcflurry", **kwargs):
        return pd.DataFrame(
            [
                {
                    "peptide": pep,
                    "allele": allele,
                    "presentation_percentile": 0.1,
                    "presentation_score": 0.9,
                    "affinity_nm": 15.0,
                    "affinity_percentile": 0.2,
                }
                for pep in peptides
                for allele in alleles
            ]
        )

    def fake_hits(peptides, **kwargs):
        return pd.DataFrame(
            [
                {
                    "peptide": pep,
                    "mhc_restriction": "HLA-A*02:01",
                    "is_monoallelic": True,
                    "pmid": "12345",
                    "src_cancer": True,
                    "src_healthy_tissue": False,
                }
                for pep in peptides
            ]
        )

    monkeypatch.setattr("tsarina.peptides.cta_exclusive_peptides", fake_peptides)
    monkeypatch.setattr("tsarina.scoring.score_presentation", fake_score)
    monkeypatch.setattr("tsarina.ms_evidence.load_public_ms_hits", fake_hits)

    bg = {"MAGEA12": "ZZZKVADLVGFMZZZ"}
    long_df = spanning_pmhc_set(
        ctas=["CTAG2", "MAGEA6"],
        alleles=["HLA-A*02:01"],
        output_format="long",
        peptides_per_cell=3,
        screen_near_self_1mm=True,
        near_self_background_proteins=bg,
    )
    by_pep = long_df.set_index("peptide")
    assert bool(by_pep.loc["SLLMWITQC", "clinical_cleared_epitope"]) is True
    assert by_pep.loc["SLLMWITQC", "cta_safety_note"] == ""
    assert "SLLMWITQC" in by_pep.loc["GARGPESRL", "cta_safety_note"]
    assert "high 1-mm near-self overlap" in by_pep.loc["KVADLVGFL", "cta_safety_note"]
    assert by_pep.loc["KVADLVGFL", "near_self_1mm_count"] == 1
    assert by_pep.loc["KVADLVGFL", "near_self_1mm_genes"] == "MAGEA12"

    table = format_panel_table(long_df)
    assert "OncoRef CTA safety caveats" in table
    assert "1-aa near-self matches" in table


def test_cli_personalize_and_panel_forward_near_self_flags(monkeypatch, capsys):
    captured: dict[str, object] = {}

    def fake_personalize(**kwargs):
        captured["personalize"] = kwargs
        return pd.DataFrame(columns=["peptide", "category", "source"])

    def fake_panel(**kwargs):
        captured["panel"] = kwargs
        out = pd.DataFrame(columns=["cta", "allele", "peptide"])
        out.attrs["panel_summary"] = None
        return out

    monkeypatch.setattr("tsarina.personalize.personalized_targets", fake_personalize)
    monkeypatch.setattr("tsarina.spanning.spanning_pmhc_set", fake_panel)

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "tsarina",
            "personalize",
            "--hla",
            "HLA-A*02:01",
            "--cta",
            "MAGEA4=10",
            "--skip-ms-evidence",
            "--no-score",
            "--screen-near-self",
            "--max-near-self-1mm",
            "0",
            "--near-self-background",
            "Unsafe|Other",
        ],
    )
    main()
    assert captured["personalize"]["screen_near_self_1mm"] is True
    assert captured["personalize"]["max_near_self_1mm_matches"] == 0
    assert captured["personalize"]["near_self_background"] == "Unsafe|Other"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "tsarina",
            "panel",
            "--ctas",
            "MAGEA4",
            "--alleles",
            "HLA-A*02:01",
            "--format",
            "long",
            "--no-summary",
            "--screen-near-self",
            "--max-near-self-1mm",
            "1",
            "--near-self-background",
            "Unsafe",
        ],
    )
    main()
    capsys.readouterr()
    assert captured["panel"]["screen_near_self_1mm"] is True
    assert captured["panel"]["max_near_self_1mm_matches"] == 1
    assert captured["panel"]["near_self_background"] == "Unsafe"
