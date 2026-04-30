import pandas as pd

from tsarina.mutations import HOTSPOT_MUTATIONS


def test_hotspot_mutations_nonempty():
    assert len(HOTSPOT_MUTATIONS) >= 15


def test_all_mutations_have_required_fields():
    required = {
        "gene",
        "gene_id",
        "transcript_id",
        "protein_position",
        "ref_aa",
        "alt_aa",
        "label",
        "cancer_types",
        "frequency_note",
    }
    for i, mut in enumerate(HOTSPOT_MUTATIONS):
        assert required.issubset(mut.keys()), f"Mutation {i} missing: {required - mut.keys()}"


def test_amino_acids_are_single_letter():
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    for mut in HOTSPOT_MUTATIONS:
        assert mut["ref_aa"] in valid_aa, f"{mut['label']}: invalid ref_aa {mut['ref_aa']}"
        assert mut["alt_aa"] in valid_aa, f"{mut['label']}: invalid alt_aa {mut['alt_aa']}"
        assert mut["ref_aa"] != mut["alt_aa"], f"{mut['label']}: ref == alt"


def test_kras_g12d_present():
    labels = {m["label"] for m in HOTSPOT_MUTATIONS}
    assert "KRAS G12D" in labels


def test_braf_v600e_present():
    labels = {m["label"] for m in HOTSPOT_MUTATIONS}
    assert "BRAF V600E" in labels


def test_idh1_r132h_present():
    labels = {m["label"] for m in HOTSPOT_MUTATIONS}
    assert "IDH1 R132H" in labels


def test_all_have_cancers():
    for mut in HOTSPOT_MUTATIONS:
        assert len(mut["cancer_types"]) > 0, f"{mut['label']} has no cancer types"


def test_gene_ids_are_ensembl():
    for mut in HOTSPOT_MUTATIONS:
        assert mut["gene_id"].startswith("ENSG"), f"{mut['label']}: bad gene_id"


def test_transcript_ids_are_ensembl():
    for mut in HOTSPOT_MUTATIONS:
        assert mut["transcript_id"].startswith("ENST"), f"{mut['label']}: bad transcript_id"


def test_mutant_iedb_overlap_uses_public_ms_loader(monkeypatch):
    from tsarina.mutations import mutant_iedb_overlap

    calls = {}
    monkeypatch.setattr(
        "tsarina.mutations.mutant_peptides",
        lambda **kw: pd.DataFrame(
            {
                "peptide": ["MSPEPTIDE", "NOEVIDENC"],
                "mutation": ["G12D", "G12D"],
            }
        ),
        raising=True,
    )

    def _fake_load_public_ms_hits(peptides, **kwargs):
        calls["peptides"] = peptides
        calls["kwargs"] = kwargs
        return pd.DataFrame(
            {
                "peptide": ["MSPEPTIDE"],
                "mhc_restriction": ["HLA-A*02:01"],
            }
        )

    monkeypatch.setattr(
        "tsarina.ms_evidence.load_public_ms_hits",
        _fake_load_public_ms_hits,
        raising=True,
    )

    out = mutant_iedb_overlap()

    assert calls["peptides"] == {"MSPEPTIDE", "NOEVIDENC"}
    assert calls["kwargs"]["drop_binding_assays"] is True
    hit_row = out[out["peptide"] == "MSPEPTIDE"].iloc[0]
    assert bool(hit_row["has_iedb_hit"]) is True
    assert hit_row["iedb_alleles"] == "HLA-A*02:01"
