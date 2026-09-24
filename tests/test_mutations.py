from types import SimpleNamespace

import pandas as pd
import pytest

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


@pytest.fixture
def mutation_reference(monkeypatch):
    """Ensembl 112 sequence windows reproducing the PLK1/RHOT2 overlaps."""
    from tsarina.mutations import HOTSPOT_MUTATIONS

    sequences = {
        "ENST00000646891": "A" * 589 + "VKIGDFGLATVKSRWSGSHQFEQLSG",
        "ENST00000256078": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVD",
    }

    def transcript(sequence, biotype="protein_coding"):
        return SimpleNamespace(protein_sequence=sequence, biotype=biotype)

    genes = [
        SimpleNamespace(biotype="protein_coding", transcripts=[transcript(seq)])
        for seq in sequences.values()
    ]
    # The overlapping protein is not the first isoform: screen every coding
    # transcript, rather than just one representative translation per gene.
    genes.append(
        SimpleNamespace(
            biotype="protein_coding",
            transcripts=[
                transcript("AAAAAAA"),
                transcript("NLFLNEDLEVKIGDFGLATKVEYDGERKKTLCGTPNYIAP"),
                transcript("GQTQRSVLLCKVVGARGVGKSAFLQAFLGRGLGHQDTREQPPGYA"),
            ],
        )
    )

    class Reference:
        def __init__(self, release):
            assert release == 112

        def transcript_by_id(self, transcript_id):
            return transcript(sequences[transcript_id])

        def genes(self):
            return iter(genes)

    monkeypatch.setattr("pyensembl.EnsemblRelease", Reference)
    return SimpleNamespace(
        genes=genes,
        transcript=transcript,
        mutations=[m for m in HOTSPOT_MUTATIONS if m["label"] in {"BRAF V600K", "KRAS G12R"}],
        self_peptides={"KIGDFGLATK", "VKIGDFGLATK", "VVGARGVGK", "VVGARGVGKSA"},
    )


def test_mutant_peptides_exclude_other_human_proteins(mutation_reference):
    from tsarina.mutations import mutant_peptides

    out = mutant_peptides(mutations=mutation_reference.mutations)
    assert not out.empty
    assert mutation_reference.self_peptides.isdisjoint(out.peptide)
    assert "KIGDFGLATKK" in set(out.peptide)
    assert (out.peptide != out.wildtype_peptide).all()


def test_mutant_raw_enumeration_is_explicit(mutation_reference):
    from tsarina.mutations import mutant_peptides

    raw = mutant_peptides(mutations=mutation_reference.mutations, require_human_exclusive=False)
    assert mutation_reference.self_peptides <= set(raw.peptide)


def test_mutant_self_screen_includes_ig_transcripts(mutation_reference):
    from tsarina.mutations import mutant_peptides

    mutation_reference.genes.append(
        SimpleNamespace(
            biotype="IG_V_gene",
            transcripts=[mutation_reference.transcript("KIGDFGLATKK", biotype="IG_V_gene")],
        )
    )
    out = mutant_peptides(mutations=mutation_reference.mutations)
    assert "KIGDFGLATKK" not in set(out.peptide)


def test_mutant_background_errors_propagate(mutation_reference, monkeypatch):
    def broken_genes(self):
        raise RuntimeError("reference unavailable")

    monkeypatch.setattr("pyensembl.EnsemblRelease.genes", broken_genes)
    from tsarina.mutations import mutant_peptides

    with pytest.raises(RuntimeError, match="reference unavailable"):
        mutant_peptides(mutations=mutation_reference.mutations)


def test_mutant_self_peptides_never_reach_personalized_ranking(mutation_reference, monkeypatch):
    from tsarina.personalize import personalized_targets

    def score(peptides, alleles, **kwargs):
        assert mutation_reference.self_peptides.isdisjoint(peptides)
        return pd.DataFrame(
            {"peptide": peptides, "allele": alleles[0], "presentation_percentile": 0.01}
        )

    monkeypatch.setattr("tsarina.scoring.score_presentation", score)
    out = personalized_targets(
        hla_alleles=["HLA-A*03:01"],
        mutations=["BRAF V600K", "KRAS G12R"],
        skip_ms_evidence=True,
        show_progress=False,
    )
    assert not out.empty
    assert out.tier.eq(1).all()
    assert mutation_reference.self_peptides.isdisjoint(out.peptide)


def test_unified_targets_screen_mutant_self_peptides(mutation_reference):
    from tsarina.targets import target_peptides

    out = target_peptides(cta=False, mutations=True)
    assert not out.empty
    assert mutation_reference.self_peptides.isdisjoint(out.peptide)


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
