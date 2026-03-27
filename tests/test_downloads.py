from ctabase.downloads import (
    FETCHABLE_DATASETS,
    MANUAL_DATASETS,
    available_datasets,
)


def test_all_fetchable_have_required_keys():
    required = {"url", "filename", "description"}
    for name, info in FETCHABLE_DATASETS.items():
        assert required.issubset(info.keys()), f"{name} missing: {required - info.keys()}"


def test_all_manual_have_required_keys():
    required = {"download_url", "description", "expected_filename"}
    for name, info in MANUAL_DATASETS.items():
        assert required.issubset(info.keys()), f"{name} missing: {required - info.keys()}"


def test_available_datasets_includes_all():
    available = available_datasets()
    for name in FETCHABLE_DATASETS:
        assert name in available
    for name in MANUAL_DATASETS:
        assert name in available


def test_viral_proteomes_in_fetchable():
    viruses = {"hpv16", "hpv18", "ebv", "htlv1", "hbv", "hcv", "kshv", "mcpyv", "hiv1"}
    assert viruses.issubset(FETCHABLE_DATASETS.keys())


def test_iedb_cedar_in_manual():
    assert "iedb" in MANUAL_DATASETS
    assert "cedar" in MANUAL_DATASETS
