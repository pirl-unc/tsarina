from perseo.mutations import HOTSPOT_MUTATIONS


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
