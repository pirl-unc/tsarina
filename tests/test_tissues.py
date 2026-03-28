from perseo.tissues import (
    CORE_REPRODUCTIVE_TISSUES,
    EXTENDED_REPRODUCTIVE_TISSUES,
    HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS,
    PERMISSIVE_REPRODUCTIVE_TISSUES,
    adaptive_rna_threshold,
    is_tissue_restricted,
    nonreproductive_tissues,
)


def test_core_reproductive_has_three_tissues():
    assert frozenset({"testis", "ovary", "placenta"}) == CORE_REPRODUCTIVE_TISSUES


def test_extended_is_superset_of_core():
    assert CORE_REPRODUCTIVE_TISSUES < EXTENDED_REPRODUCTIVE_TISSUES


def test_permissive_is_superset_of_extended():
    assert EXTENDED_REPRODUCTIVE_TISSUES < PERMISSIVE_REPRODUCTIVE_TISSUES


def test_is_tissue_restricted_core():
    assert is_tissue_restricted({"testis", "ovary"})
    assert not is_tissue_restricted({"testis", "liver"})


def test_thymus_excluded_by_default():
    assert is_tissue_restricted({"testis", "thymus"})


def test_thymus_not_excluded_when_requested():
    assert not is_tissue_restricted({"testis", "thymus"}, exclude_thymus=False)


def test_nonreproductive_core():
    result = nonreproductive_tissues({"testis", "liver", "brain"})
    assert result == frozenset({"liver", "brain"})


def test_nonreproductive_extended():
    result = nonreproductive_tissues(
        {"testis", "prostate", "liver"},
        definition=EXTENDED_REPRODUCTIVE_TISSUES,
    )
    assert result == frozenset({"liver"})


def test_adaptive_thresholds():
    assert adaptive_rna_threshold("Enhanced") == 0.80
    assert adaptive_rna_threshold("Supported") == 0.90
    assert adaptive_rna_threshold("Approved") == 0.95
    assert adaptive_rna_threshold("Uncertain") == 0.99
    assert adaptive_rna_threshold("no data") == 0.99


def test_threshold_order_is_monotonic():
    tiers = ["Enhanced", "Supported", "Approved", "Uncertain"]
    thresholds = [HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS[t] for t in tiers]
    assert thresholds == sorted(thresholds)
