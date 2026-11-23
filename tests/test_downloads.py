from tsarina.downloads import (
    DATASET_CATEGORIES,
    FETCHABLE_DATASETS,
    MANUAL_DATASETS,
    available_dataset_catalog,
    available_datasets,
    dataset_categories,
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


def test_iedb_cedar_now_fetchable():
    # IEDB/CEDAR moved MANUAL -> FETCHABLE in hitlist 1.43.0: the downloader.php
    # endpoints serve the zip directly, so `tsarina data fetch iedb` works.
    assert "iedb" in FETCHABLE_DATASETS
    assert "cedar" in FETCHABLE_DATASETS
    assert "iedb" not in MANUAL_DATASETS
    assert "cedar" not in MANUAL_DATASETS


def test_structured_dataset_catalog_and_categories_reexported():
    assert dataset_categories() == DATASET_CATEGORIES
    assert "core_immunopeptidomics" in DATASET_CATEGORIES
    core = available_dataset_catalog(category="core_immunopeptidomics")
    assert {"iedb", "cedar"} <= set(core)
    assert core["iedb"]["category"] == "core_immunopeptidomics"
