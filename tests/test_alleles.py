from ctabase.alleles import (
    GLOBAL44_ADDON,
    GLOBAL48_ADDON,
    GLOBAL51_SSA_ADDON,
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
    assert names == ["iedb27_ab", "iedb36_abc", "global44_abc", "global48_abc", "global51_abc_ssa"]


def test_panels_are_nested_supersets():
    names = panel_names()
    for i in range(1, len(names)):
        smaller = set(get_panel(names[i - 1]))
        larger = set(get_panel(names[i]))
        assert smaller < larger, f"{names[i - 1]} is not a strict subset of {names[i]}"


def test_global51_has_51_alleles():
    panel = get_panel("global51_abc_ssa")
    assert len(panel) == 51


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
