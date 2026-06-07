from pathlib import Path

import pandas as pd
import pytest

from tsarina.ms_evidence import (
    aggregate_ms_hits_by_peptide,
    aggregate_ms_hits_for_iedb_columns,
    load_public_ms_hits,
)


def test_load_public_ms_hits_raw_path_drops_binding_assays(monkeypatch):
    calls = {}

    def _fake_resolve(iedb_path, cedar_path, require_iedb=True):
        calls["resolve"] = (iedb_path, cedar_path, require_iedb)
        return Path("iedb.csv"), None

    def _fake_scan_public_ms(**kwargs):
        calls["scan"] = kwargs
        return pd.DataFrame(
            {
                "peptide": ["MSPEPTIDE", "BINDPEPTD"],
                "mhc_restriction": ["HLA-A*02:01", "HLA-A*02:01"],
                "is_binding_assay": [False, True],
            }
        )

    monkeypatch.setattr("tsarina.datasources.resolve_dataset_paths", _fake_resolve, raising=True)
    monkeypatch.setattr("tsarina.iedb.scan_public_ms", _fake_scan_public_ms, raising=True)

    out = load_public_ms_hits(
        peptides={"MSPEPTIDE", "BINDPEPTD"},
        iedb_path="iedb.csv",
        mhc_class="I",
    )

    assert calls["resolve"] == ("iedb.csv", None, True)
    assert calls["scan"]["mhc_species"] == "Homo sapiens"
    assert out["peptide"].tolist() == ["MSPEPTIDE"]


def test_aggregate_ms_hits_by_peptide_carries_source_flags_and_cell_lines():
    hits = pd.DataFrame(
        {
            "peptide": ["MSPEPTIDE", "MSPEPTIDE"],
            "mhc_restriction": ["HLA-A*02:01", "HLA-B*07:02"],
            "src_cancer": [True, False],
            "src_healthy_tissue": [False, True],
            "cell_line_name": ["A375", ""],
        }
    )

    out = aggregate_ms_hits_by_peptide(hits)

    row = out.iloc[0]
    assert row["ms_hit_count"] == 2
    assert row["ms_alleles"] == "HLA-A*02:01;HLA-B*07:02"
    assert row["ms_allele_count"] == 2
    assert bool(row["ms_cancer"]) is True
    assert bool(row["ms_healthy_tissue"]) is True
    assert row["ms_cell_lines"] == "A375"


def test_aggregate_ms_hits_by_peptide_ignores_blank_and_missing_strings():
    hits = pd.DataFrame(
        {
            "peptide": ["MSPEPTIDE", "MSPEPTIDE", "MSPEPTIDE"],
            "mhc_restriction": ["HLA-A*02:01", None, ""],
            "cell_line_name": ["A375", None, ""],
        }
    )

    out = aggregate_ms_hits_by_peptide(hits)

    row = out.iloc[0]
    assert row["ms_hit_count"] == 3
    assert row["ms_alleles"] == "HLA-A*02:01"
    assert row["ms_allele_count"] == 1
    assert row["ms_cell_lines"] == "A375"


def test_aggregate_ms_hits_for_iedb_columns_keeps_legacy_shape():
    hits = pd.DataFrame(
        {
            "peptide": ["MSPEPTIDE"],
            "mhc_restriction": ["HLA-A*02:01"],
        }
    )

    out = aggregate_ms_hits_for_iedb_columns(hits)

    assert out.to_dict(orient="records") == [
        {
            "peptide": "MSPEPTIDE",
            "iedb_hit_count": 1,
            "iedb_alleles": "HLA-A*02:01",
        }
    ]


def test_cta_healthy_tissue_ms_hits_surfaces_tissue_allele_pmid_vital(monkeypatch):
    """Per-peptide healthy-somatic MS hits with vital-organ flag (tsarina#76)."""
    import tsarina.indexing

    rows = pd.DataFrame(
        {
            "peptide": ["HEARTPEP", "BLOODPEP", "CANCERPEP", "TESTISPEP"],
            "mhc_restriction": ["HLA-B*49:01", "HLA-B*44:03", "HLA-A*02:01", "HLA-A*01:01"],
            "mhc_allele_set": ["HLA-B*49:01", "HLA-B*44:03", "HLA-A*02:01", "HLA-A*01:01"],
            "mhc_allele_provenance": ["exact", "sample_allele_match", "exact", "exact"],
            "pmid": ["33858848", "38920720", "111", "222"],
            "source_tissue": ["Heart", "Blood", "Melanoma", "Testis"],
            "cell_name": ["", "", "", ""],
            # only the first two are healthy somatic
            "src_healthy_tissue": [True, True, False, False],
        }
    )
    monkeypatch.setattr(tsarina.indexing, "load_ms_evidence", lambda **kw: rows, raising=True)

    from tsarina import cta_healthy_tissue_ms_hits

    out = cta_healthy_tissue_ms_hits("MAGEA4")
    # Only the two healthy-somatic rows survive (cancer/testis dropped).
    assert list(out["peptide"]) == ["HEARTPEP", "BLOODPEP"]
    assert list(out["tissue"]) == ["Heart", "Blood"]
    assert list(out["pmid"]) == ["33858848", "38920720"]
    # Heart is a vital organ; blood is not.
    assert list(out["vital_organ"]) == [True, False]
    assert set(out.columns) == {
        "peptide",
        "tissue",
        "allele",
        "allele_set",
        "provenance",
        "pmid",
        "vital_organ",
    }


def test_cta_healthy_tissue_ms_hits_empty_when_no_healthy_hits(monkeypatch):
    import tsarina.indexing

    rows = pd.DataFrame(
        {
            "peptide": ["CANCERPEP"],
            "mhc_restriction": ["HLA-A*02:01"],
            "source_tissue": ["Melanoma"],
            "src_healthy_tissue": [False],
        }
    )
    monkeypatch.setattr(tsarina.indexing, "load_ms_evidence", lambda **kw: rows, raising=True)

    from tsarina import cta_healthy_tissue_ms_hits

    out = cta_healthy_tissue_ms_hits("PRAME")
    assert out.empty


# ── Vital-organ matcher (_is_vital_organ_tissue) ──────────────────────────────


def test_vital_organ_matcher_covers_safety_source_of_truth():
    """The MS screen must flag every tissue in the RNA-side vital vocabulary.

    Guards against drift: if a tissue is added to SAFETY_TISSUE_GROUPS or
    VITAL_TISSUE_MS_NAMES, this fails until the MS screen covers it too.
    """
    from tsarina.ms_evidence import _is_vital_organ_tissue
    from tsarina.tiers import SAFETY_TISSUE_GROUPS, VITAL_TISSUE_MS_NAMES

    source_of_truth = set().union(*SAFETY_TISSUE_GROUPS.values()) | set(VITAL_TISSUE_MS_NAMES)
    missed = sorted(name for name in source_of_truth if not _is_vital_organ_tissue(name))
    assert missed == []


@pytest.mark.parametrize(
    "tissue",
    [
        # heart
        "Heart",
        "heart muscle",
        "Cardiac muscle",
        "myocardium",
        "left ventricle myocardium",
        # lung
        "Lung",
        "pulmonary",
        # liver
        "Liver",
        "livers",
        "hepatic tissue",
        "hepatocyte",
        # pancreas
        "Pancreas",
        "pancreatic islets",
        "islets of Langerhans",
        # brain / CNS (incl. cortex spellings and subregions)
        "Brain",
        "midbrain",
        "brainstem",
        "Cerebellum",
        "cerebellar cortex",
        "Cerebral cortex",
        "frontal cortex",
        "visual cortex",
        "Central nervous system (CNS)",
        "CNS",
        "spinal cord",
        "Hippocampus",
        "hippocampal formation",
        "amygdala",
        "thalamus",
        "hypothalamus",
        "pons",
        "basal ganglia",
        "choroid plexus",
        "medulla oblongata",
        "white matter",
        "retina",
    ],
)
def test_vital_organ_matcher_true(tissue):
    from tsarina.ms_evidence import _is_vital_organ_tissue

    assert _is_vital_organ_tissue(tissue) is True


@pytest.mark.parametrize(
    "tissue",
    [
        "",
        "Blood",
        "Bone Marrow",
        "Skin",
        "Breast",
        "Colon",
        "Kidney",
        "Stomach",
        "Lymph Node",
        # non-CNS cortices must NOT masquerade as brain
        "adrenal cortex",
        "renal cortex",
        "kidney cortex",
        "ovarian cortex",
        "lymph node cortex",
        "thymic cortex",
        # bare "medulla" names adrenal/renal medulla, not brain
        "adrenal medulla",
        "renal medulla",
        # expected-tissue / unrelated
        "thymus",
        "testis",
        "ovary",
        "placenta",
        # substring traps
        "delivery",
        "striated muscle",
        "skeletal muscle",
        "gastric cardia",
        "melanoma",
    ],
)
def test_vital_organ_matcher_false(tissue):
    from tsarina.ms_evidence import _is_vital_organ_tissue

    assert _is_vital_organ_tissue(tissue) is False
