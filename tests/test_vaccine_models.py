"""Real final-sequence affinity and cleavage contracts (explicit model opt-in)."""

import pytest

from tsarina.vaccine_construct import PredictionAudit, VaccineConfig
from tsarina.vaccine_sequences import Segment, assemble_layers

pytestmark = pytest.mark.mhcflurry


def test_real_final_junction_binding_and_local_cleavage_context(monkeypatch):
    from mhctools import Pepsickle

    import tsarina.scoring as scoring

    monkeypatch.setattr(scoring, "_MHCFLURRY_PRESENTATION_PREDICTOR", None)

    protein = "MALWGQDPYVKVLEHVVRVATKVRRQ"
    a = Segment("a", "a", "a", 1, 1.0, protein, 0, 12, 0, 12, (), ())
    b = Segment("b", "b", "b", 2, 0.5, protein, 12, len(protein), 12, len(protein), (), ())
    sequence, layers = assemble_layers([(a, 0, 0, ""), (b, 0, 0, "AAY")])
    audit = PredictionAudit(["HLA-A*02:01", "HLA-A*24:02"], VaccineConfig())
    audit.prepare([(sequence, layers)])
    key, junctions, cleavage = audit.assess(sequence, layers)
    assert junctions and len(cleavage) == 2
    assert all(row["affinity_nm"] > 0 for row in junctions)
    assert key[0] == sum(row["below_threshold"] for row in junctions)
    assert audit.cleavage_model["name"]
    # Local 17-aa contexts must agree with the actual complete protein profile;
    # using 16 residues would incorrectly pad one required upstream residue.
    full = Pepsickle(human_only=True, isolate_subprocess=True).cleavage_probs(sequence)
    for row in cleavage:
        assert row["cleavage_probability"] == pytest.approx(full[row["bond"] - 1], abs=1e-7)


def test_real_codon_encoding_with_utr_and_polya():
    from tsarina.vaccine_construct import encode_construct, nucleotide_elements

    cfg = VaccineConfig(include_utrs=True, poly_a_length=120)
    utr5, utr3, polya = nucleotide_elements(cfg)
    cds, full = encode_construct(
        {"protein": "MACDEFGHIK", "utr5": utr5, "utr3": utr3, "poly_a": polya}, cfg
    )
    assert len(utr5) == 50 and len(utr3) == 268
    assert len(cds) == 33 and len(full) == 471
    assert cds.startswith("AUG") and cds.endswith("UAA") and "T" not in full
