"""Tests for composable OncoRef HostContext and 4-way sex_restriction_code filtering (#260)."""

from __future__ import annotations

import json
import sys
import warnings

import pandas as pd
import pytest

import tsarina
from tsarina.cli import main
from tsarina.cta_selection import (
    SEX_RESTRICTION_CODES,
    CTA_selection,
    CTA_sex_eligible_gene_ids,
    CTA_sex_eligible_gene_names,
    CTA_sex_restricted_antigens,
    HostContext,
    build_host_context,
    host_compartments,
    selection_representative_transcripts,
)
from tsarina.evidence import CTA_evidence
from tsarina.personalize import personalized_targets


def test_top_level_exports_host_context_and_sex_codes():
    assert tsarina.HostContext is HostContext
    assert tsarina.SEX_RESTRICTION_CODES == SEX_RESTRICTION_CODES
    assert tsarina.build_host_context is build_host_context
    assert tsarina.host_compartments is host_compartments


def test_four_way_sex_code_filtering_across_human_mouse_and_dog():
    # Human CTA-F includes LALBA; CTA-P includes DUX4; RA-M-only includes PRM2.
    cta_f = CTA_selection(sex_code="cta-f")
    assert not cta_f.empty
    assert set(cta_f["sex_restriction_code"]) == {"CTA-F"}
    assert "LALBA" in set(cta_f["label"])

    cta_p = CTA_selection(sex_code=["CTA-P"])
    assert not cta_p.empty
    assert set(cta_p["sex_restriction_code"]) == {"CTA-P"}
    assert "DUX4" in set(cta_p["label"])

    ra_m = CTA_selection(sex_code="RA-M-only")
    assert not ra_m.empty
    assert set(ra_m["sex_restriction_code"]) == {"RA-M-only"}
    assert "PRM2" in set(ra_m["label"])

    # Comma-separated multi-code filter works on CTA_sex_restricted_antigens and gene helpers.
    multi = CTA_sex_restricted_antigens(sex="female", sex_code="CTA-M,CTA-P")
    assert set(multi["sex_restriction_code"]) <= {"CTA-M", "CTA-P"}
    assert "MAGEA4" in set(multi["label"])

    names = CTA_sex_eligible_gene_names("female", sex_code="CTA-P")
    ids = CTA_sex_eligible_gene_ids("female", sex_code="CTA-P")
    assert "DUX4" in names
    assert len(ids) >= len(names)

    # Mouse and dog support the exact same 4-way codes.
    mouse_m = CTA_selection(species="mouse", sex_code="CTA-M")
    assert not mouse_m.empty
    assert set(mouse_m["sex_restriction_code"]) == {"CTA-M"}

    dog_p = CTA_selection(species="dog", sex_code="CTA-P")
    assert not dog_p.empty
    assert set(dog_p["sex_restriction_code"]) == {"CTA-P"}

    with pytest.raises(ValueError, match="Unsupported sex_restriction_code"):
        CTA_selection(sex_code="CTA-INVALID")


def test_host_context_filtering_across_human_mouse_and_dog():
    # Human lactating female excludes lactation CTAs (LALBA, CSN2, BTN1A1) while keeping MAGEA4.
    lact_ctx = build_host_context(patient_sex="female", lactating=True, species="human")
    assert isinstance(lact_ctx, HostContext)
    assert "lactating_breast" in lact_ctx.protected

    assessed = CTA_selection(category="CTA", host_context=lact_ctx)
    assert {"host_status", "host_reasons", "pass_host"} <= set(assessed.columns)
    assert assessed.attrs["host_context"]["species"] == "human"
    excluded = set(assessed.loc[~assessed["pass_host"], "label"])
    assert {"LALBA", "CSN2", "BTN1A1"} <= excluded
    assert "MAGEA4" not in excluded

    filtered = CTA_selection(category="CTA", host_context=lact_ctx, filter_host=True)
    assert "LALBA" not in set(filtered["label"])
    assert "MAGEA4" in set(filtered["label"])

    # Human protect_tissues=("testis",) excludes testis CTAs while retaining non-testis CTAs.
    testis_ctx = build_host_context(protect_tissues=("testis",), species="human")
    testis_kept = CTA_selection(category="CTA", host_context=testis_ctx, filter_host=True)
    assert 0 < len(testis_kept) < len(assessed)
    assert "LALBA" in set(testis_kept["label"])

    # Mouse protect testis filters mouse CTAs.
    mouse_ctx = build_host_context(protect_tissues=("testis",), species="mouse")
    mouse_kept = CTA_selection(
        category="CTA", species="mouse", host_context=mouse_ctx, filter_host=True
    )
    assert 0 < len(mouse_kept) < len(CTA_selection(category="CTA", species="mouse"))

    # Dog lactating protects lactating_mammary and excludes canine casein CTAs.
    dog_ctx = build_host_context(patient_sex="female", lactating=True, species="dog")
    assert "lactating_mammary" in dog_ctx.protected
    dog_assessed = CTA_selection(category="CTA", species="dog", host_context=dog_ctx)
    dog_excluded = set(dog_assessed.loc[~dog_assessed["pass_host"], "label"])
    assert "CSN1S1" in dog_excluded or "CSN3" in dog_excluded

    # Species mismatch fails closed.
    with pytest.raises(ValueError, match="does not match"):
        CTA_selection(category="CTA", species="human", host_context=dog_ctx)


def test_personalized_targets_applies_host_context_for_human_and_dog(monkeypatch):
    def fake_cta_exclusive(**kwargs):
        return pd.DataFrame(
            [
                {
                    "gene_name": "MAGEA4",
                    "gene_id": "ENSG00000147381",
                    "peptide": "GVYDGREHTV",
                    "length": 10,
                },
                {
                    "gene_name": "LALBA",
                    "gene_id": "ENSG00000167531",
                    "peptide": "KQFTKCELS",
                    "length": 9,
                },
            ]
        )

    monkeypatch.setattr("tsarina.peptides.cta_exclusive_peptides", fake_cta_exclusive)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = personalized_targets(
            hla_alleles=["HLA-A*02:01"],
            cta_expression={"MAGEA4": 50.0, "LALBA": 100.0},
            lactating=True,
            score_presentation=False,
            skip_ms_evidence=True,
            drop_weak_tier=False,
            show_progress=False,
        )

    assert set(out["source"]) == {"MAGEA4"}
    assert "LALBA" not in set(out["source"])
    assert "host_context" in out.attrs
    assert "lactating_breast" in out.attrs["host_context"]["protected"]
    assert any("LALBA" in str(w.message) and "HostContext" in str(w.message) for w in caught)

    # Pregnant protects embryo/placenta/lactating_breast and drops PLAC1/MAGEA4 conservatively.
    with warnings.catch_warnings(record=True) as preg_warns:
        warnings.simplefilter("always")
        preg_out = personalized_targets(
            hla_alleles=["HLA-A*02:01"],
            cta_expression={"MAGEA4": 50.0},
            pregnant=True,
            score_presentation=False,
            skip_ms_evidence=True,
            drop_weak_tier=False,
            show_progress=False,
        )
    assert preg_out.empty
    assert "host_context" in preg_out.attrs
    assert "embryo" in preg_out.attrs["host_context"]["protected"]
    assert any("MAGEA4" in str(w.message) and "HostContext" in str(w.message) for w in preg_warns)

    # Passing both explicit host_context and individual flags raises ValueError.
    ctx = build_host_context(lactating=True, species="human")
    with pytest.raises(ValueError, match="not both"):
        personalized_targets(
            hla_alleles=["HLA-A*02:01"],
            cta_expression={"MAGEA4": 50.0},
            host_context=ctx,
            lactating=True,
            score_presentation=False,
            skip_ms_evidence=True,
            show_progress=False,
        )

    # Dog HostContext drops canine lactating-mammary CTA CSN1S1.
    monkeypatch.setattr(
        "tsarina.peptides.gene_peptides",
        lambda gene_ids, **kwargs: pd.DataFrame(
            [{"gene_name": "CSN1S1", "gene_id": gene_ids[0], "peptide": "AAAAAAAAA", "length": 9}]
        ),
    )
    with warnings.catch_warnings(record=True) as dog_warns:
        warnings.simplefilter("always")
        dog_out = personalized_targets(
            hla_alleles=["DLA-88*508:01"],
            cta_expression={"CSN1S1": 80.0},
            host_species="dog",
            lactating=True,
            host_proteins={"BG": "MMMMMMMMMMMM"},
            score_presentation=False,
            skip_ms_evidence=True,
            drop_weak_tier=False,
            show_progress=False,
        )
    assert dog_out.empty
    assert "host_context" in dog_out.attrs
    assert dog_out.attrs["host_context"]["species"] == "dog"
    assert any("CSN1S1" in str(w.message) and "HostContext" in str(w.message) for w in dog_warns)


def test_cli_catalog_sex_code_host_context_and_personalize_forwarding(monkeypatch, capsys):
    # catalog ctas --sex-code CTA-F --json
    monkeypatch.setattr(
        sys, "argv", ["tsarina", "catalog", "ctas", "--sex-code", "CTA-F", "--top", "50", "--json"]
    )
    main()
    payload = json.loads(capsys.readouterr().out)
    assert payload["sex_codes"] == ["CTA-F"]
    assert payload["total"] > 0
    assert all(row["sex_restriction_code"] == "CTA-F" for row in payload["ctas"])
    assert any(row["label"] == "LALBA" for row in payload["ctas"])

    # catalog ctas --sex-code RA-M-only --json (does not get zeroed out by default category=CTA)
    monkeypatch.setattr(
        sys,
        "argv",
        ["tsarina", "catalog", "ctas", "--sex-code", "RA-M-only", "--top", "200", "--json"],
    )
    main()
    ra_payload = json.loads(capsys.readouterr().out)
    assert ra_payload["total"] > 0
    assert any(row["label"] == "PRM2" for row in ra_payload["ctas"])

    # catalog ctas --lactating --sex-code CTA-F --json excludes LALBA/CSN2/BTN1A1
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "tsarina",
            "catalog",
            "ctas",
            "--sex-code",
            "CTA-F",
            "--lactating",
            "--top",
            "50",
            "--json",
        ],
    )
    main()
    lact_payload = json.loads(capsys.readouterr().out)
    assert "host_context" in lact_payload
    assert not any(row["label"] in {"LALBA", "CSN2", "BTN1A1"} for row in lact_payload["ctas"])

    # personalize CLI forwards all HostContext flags cleanly
    captured: dict[str, object] = {}

    def fake_personalize(**kwargs):
        captured.update(kwargs)
        return pd.DataFrame(columns=["peptide", "category", "source"])

    monkeypatch.setattr("tsarina.personalize.personalized_targets", fake_personalize)
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
            "--sex",
            "female",
            "--lactating",
            "--protect-tissue",
            "ovary,uterus",
            "--allow-tissue",
            "cervix",
            "--absent-tissue",
            "testis",
            "--skip-ms-evidence",
            "--no-score",
        ],
    )
    main()
    capsys.readouterr()
    assert captured["patient_sex"] == "female"
    assert captured["lactating"] is True
    assert captured["pregnant"] is False
    assert captured["protect_tissues"] == ["ovary", "uterus"]
    assert captured["accept_loss_tissues"] == ["cervix"]
    assert captured["absent_tissues"] == ["testis"]
    assert captured["allow_incomplete_scope"] is False


def test_allow_incomplete_scope_ectopic_replicated_and_cross_species_transcripts(
    monkeypatch, capsys
):
    # Pregnant fail-closed (default) vs allow_incomplete_scope=True
    preg_closed = build_host_context(pregnant=True, species="human")
    preg_open = build_host_context(pregnant=True, allow_incomplete_scope=True, species="human")
    assert preg_closed is not None and preg_closed.allow_incomplete_scope is False
    assert preg_open is not None and preg_open.allow_incomplete_scope is True

    closed_ctas = CTA_selection(category="CTA", host_context=preg_closed, filter_host=True)
    open_ctas = CTA_selection(category="CTA", host_context=preg_open, filter_host=True)
    assert len(closed_ctas) == 0
    assert len(open_ctas) >= 100
    assert "NY-ESO-1" in set(open_ctas["label"])
    assert "PLAC1" not in set(open_ctas["label"])

    # ectopic_only and replicated_only filters across human, mouse, and dog
    all_human = CTA_selection(category="CTA")
    ectopic_human = CTA_selection(category="CTA", ectopic_only=True)
    replicated_human = CTA_selection(category="CTA", replicated_only=True)
    assert 0 < len(ectopic_human) < len(all_human)
    assert set(ectopic_human["tumor_activation_mode"]) <= {
        "ectopic",
        "ectopic_and_lineage_retained",
    }
    assert 0 < len(replicated_human) < len(all_human)
    assert set(replicated_human["tumor_support"]) == {"replicated"}

    dog_ectopic = CTA_selection(category="CTA", species="dog", ectopic_only=True)
    dog_replicated = CTA_selection(category="CTA", species="dog", replicated_only=True)
    assert 0 < len(dog_ectopic) < len(CTA_selection(category="CTA", species="dog"))
    assert 0 < len(dog_replicated) < len(CTA_selection(category="CTA", species="dog"))
    assert "PRAME" in set(dog_ectopic["label"])

    # Evidence table projects tumor_activation_mode / tumor_support / qualifying_cancer_codes
    ev = CTA_evidence()
    for col in (
        "selection_tumor_activation_mode",
        "selection_tumor_support",
        "selection_qualifying_cancer_codes",
    ):
        assert col in ev.columns
    cta_ev = ev.loc[ev["is_current_cta"]]
    assert cta_ev["selection_tumor_activation_mode"].notna().all()
    assert cta_ev["selection_tumor_support"].notna().all()

    # Cross-species representative_transcripts are populated for human, mouse, and dog
    for sp, prefix in (("human", "ENST"), ("mouse", "ENSMUST"), ("dog", "ENSCAFT")):
        rep = selection_representative_transcripts(species=sp)
        assert len(rep) > 0
        sample_tx = next(iter(rep.values()))
        assert all(tid.startswith(prefix) and "." not in tid for tid in sample_tx)

    # CLI catalog ctas supports --ectopic-only, --replicated, and --pregnant --allow-incomplete-scope
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "tsarina",
            "catalog",
            "ctas",
            "--pregnant",
            "--allow-incomplete-scope",
            "--ectopic-only",
            "--replicated",
            "--top",
            "200",
            "--json",
        ],
    )
    main()
    payload = json.loads(capsys.readouterr().out)
    assert payload["host_context"]["allow_incomplete_scope"] is True
    assert payload["total"] > 50
    assert all(
        row["tumor_activation_mode"] in {"ectopic", "ectopic_and_lineage_retained"}
        for row in payload["ctas"]
    )
    assert all(row["tumor_support"] == "replicated" for row in payload["ctas"])

    # personalize CLI forwards --allow-incomplete-scope
    captured: dict[str, object] = {}

    def fake_personalize(**kwargs):
        captured.update(kwargs)
        return pd.DataFrame(columns=["peptide", "category", "source"])

    monkeypatch.setattr("tsarina.personalize.personalized_targets", fake_personalize)
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
            "--pregnant",
            "--allow-incomplete-scope",
            "--skip-ms-evidence",
            "--no-score",
        ],
    )
    main()
    capsys.readouterr()
    assert captured["pregnant"] is True
    assert captured["allow_incomplete_scope"] is True
