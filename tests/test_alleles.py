import pytest

from tsarina.alleles import (
    GLOBAL44_ADDON,
    GLOBAL48_ADDON,
    GLOBAL51_AB_BACKBONE,
    GLOBAL51_COMMON_AB_COMPLEMENT,
    GLOBAL51_HLA_C,
    GLOBAL51_SSA_ADDON,
    GLOBAL53_CTA_MS_ADDON,
    GLOBAL53_HLA_C,
    HLA_C_ADDON,
    IEDB27_AB,
    PANEL_DEFINITIONS,
    PANEL_SOURCE_CATEGORIES,
    get_panel,
    panel_names,
)


def test_iedb27_has_27_alleles():
    assert len(IEDB27_AB) == 27


def test_hla_c_addon_has_9_alleles():
    assert len(HLA_C_ADDON) == 9


def test_panel_names_order():
    names = panel_names()
    assert names == [
        "iedb27_ab",
        "iedb36_abc",
        "global44_abc",
        "global48_abc",
        "global51_abc_ssa",
        "global51_abc",
        "global53_abc",
    ]


def test_panels_are_nested_supersets():
    names = ["iedb27_ab", "iedb36_abc", "global44_abc", "global48_abc", "global51_abc_ssa"]
    for i in range(1, len(names)):
        smaller = set(get_panel(names[i - 1]))
        larger = set(get_panel(names[i]))
        assert smaller < larger, f"{names[i - 1]} is not a strict subset of {names[i]}"


def test_global51_has_51_alleles():
    panel = get_panel("global51_abc_ssa")
    assert len(panel) == 51


def test_global51_has_51_alleles_for_reference_global_panel():
    panel = get_panel("global51_abc")
    assert len(panel) == 51


def test_global53_has_53_alleles_for_default_global_panel():
    panel = get_panel("global53_abc")
    assert len(panel) == 53


def test_global51_keeps_reference_backbone():
    panel = set(get_panel("global51_abc"))
    assert set(IEDB27_AB) <= panel
    assert "HLA-A*24:02" in panel
    assert "HLA-A*24:03" not in panel


def test_global51_includes_frequent_hla_c_panel():
    panel = set(get_panel("global51_abc"))
    assert set(GLOBAL51_HLA_C) <= panel
    assert "HLA-C*14:02" in panel
    assert "HLA-C*14:03" in panel


def test_global51_ab_complement_is_reference_backed():
    panel = set(get_panel("global51_abc"))
    assert set(GLOBAL51_COMMON_AB_COMPLEMENT) <= panel
    assert GLOBAL51_COMMON_AB_COMPLEMENT == [
        "HLA-B*18:01",
        "HLA-B*40:02",
        "HLA-B*46:01",
    ]
    assert not (
        panel
        & {
            "HLA-A*02:11",  # India-specific local add-on, not in the IEDB 38 A/B set.
            "HLA-A*74:01",  # SSA-focused local add-on, not in the IEDB 38 A/B set.
            "HLA-B*14:02",  # IEDB 38 A/B candidate below the selected top-three cap.
            "HLA-B*15:03",  # SSA-focused local add-on, not in the IEDB 38 A/B set.
            "HLA-B*50:01",  # Kuwait-only local add-on, not in the IEDB 38 A/B set.
            "HLA-B*58:02",  # SSA-focused local add-on, not in the IEDB 38 A/B set.
        }
    )


def test_global53_default_adds_cta_ms_supported_alleles():
    panel = set(get_panel("global53_abc"))
    assert set(get_panel("global51_abc")) - {"HLA-C*14:03"} < panel
    assert set(GLOBAL53_CTA_MS_ADDON) <= panel
    assert GLOBAL53_CTA_MS_ADDON == [
        "HLA-A*29:02",
        "HLA-B*15:02",
        "HLA-B*27:05",
    ]


def test_global53_default_keeps_one_c14_representative():
    panel = set(get_panel("global53_abc"))
    assert "HLA-C*14:02" in panel
    assert "HLA-C*14:03" not in panel
    assert set(GLOBAL53_HLA_C) <= panel


def test_global53_default_uses_mhcflurry_runtime_calibration_when_available():
    try:
        mhcflurry = pytest.importorskip("mhcflurry")
    except Exception as e:
        pytest.skip(f"mhcflurry is installed but not importable: {e}")

    predictor = mhcflurry.Class1AffinityPredictor.load()
    missing = sorted(
        allele
        for allele in set(get_panel("global53_abc"))
        if predictor.percent_rank_calibrated_allele(allele) is None
    )
    assert missing == []
    assert predictor.percent_rank_calibrated_allele("HLA-A*24:02") is not None
    assert (
        predictor.allele_to_sequence["HLA-C*14:02"] == predictor.allele_to_sequence["HLA-C*14:03"]
    )
    assert predictor.percent_rank_calibrated_allele("HLA-C*14:03") == "HLA-C*14:02"
    assert predictor.percent_rank_calibrated_allele("HLA-C*15:05") is None


def test_all_alleles_have_source_category():
    all_alleles = set()
    for defn in PANEL_DEFINITIONS.values():
        all_alleles.update(defn["alleles"])
    for allele in all_alleles:
        assert allele in PANEL_SOURCE_CATEGORIES, f"Missing source category for {allele}"


def test_alleles_are_sorted():
    for name, defn in PANEL_DEFINITIONS.items():
        alleles = defn["alleles"]
        expected = sorted(alleles, key=lambda a: (a.split("*", 1)[0], a))
        assert alleles == expected, f"Panel {name} alleles are not sorted"


def test_all_alleles_start_with_hla():
    for name, defn in PANEL_DEFINITIONS.items():
        for allele in defn["alleles"]:
            assert allele.startswith("HLA-"), f"{allele} in {name} does not start with HLA-"


def test_addon_lists_no_overlap():
    addons = [HLA_C_ADDON, GLOBAL44_ADDON, GLOBAL48_ADDON, GLOBAL51_SSA_ADDON]
    for i, a in enumerate(addons):
        for j, b in enumerate(addons):
            if i != j:
                overlap = set(a) & set(b)
                assert not overlap, f"Addon lists {i} and {j} overlap: {overlap}"


def test_global51_default_components_have_expected_sizes():
    assert len(GLOBAL51_AB_BACKBONE) == 27
    assert len(GLOBAL51_HLA_C) == 21
    assert len(GLOBAL51_COMMON_AB_COMPLEMENT) == 3
    assert len(GLOBAL53_HLA_C) == 20
    assert len(GLOBAL53_CTA_MS_ADDON) == 3
