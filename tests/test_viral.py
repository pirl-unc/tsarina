import pandas as pd

from perseus.viral import ONCOGENIC_VIRUSES, read_fasta, viral_peptides


def test_oncogenic_viruses_has_expected_keys():
    expected = {"hpv16", "hpv18", "ebv", "htlv1", "hbv", "hcv", "kshv", "mcpyv", "hiv1"}
    assert expected == set(ONCOGENIC_VIRUSES.keys())


def test_all_viruses_have_required_fields():
    required = {
        "name",
        "full_name",
        "cancers",
        "uniprot_proteome",
        "taxonomy_id",
        "key_oncoproteins",
    }
    for key, info in ONCOGENIC_VIRUSES.items():
        assert required.issubset(info.keys()), f"{key} missing: {required - info.keys()}"


def test_all_viruses_have_cancers():
    for key, info in ONCOGENIC_VIRUSES.items():
        assert len(info["cancers"]) > 0, f"{key} has no cancers listed"


def test_viral_peptides_from_dict():
    proteins = {"test_protein": "ACDEFGHIKLMNPQRSTVWY"}  # exactly 20 aa
    df = viral_peptides(proteins=proteins, lengths=(8, 9))
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert "peptide" in df.columns
    assert "protein_id" in df.columns
    assert all(df["length"].isin([8, 9]))
    # 20 - 8 + 1 = 13 eight-mers, 20 - 9 + 1 = 12 nine-mers
    assert len(df) == 13 + 12


def test_viral_peptides_skips_non_standard_aa():
    proteins = {"test": "ACDEFXGHIK"}  # X is non-standard
    df = viral_peptides(proteins=proteins, lengths=(5,))
    # Peptides containing X should be excluded
    assert all("X" not in p for p in df["peptide"])


def test_read_fasta(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">prot1 some description\nACDEFG\nHIKLMN\n>prot2\nPQRSTV\n")
    proteins = read_fasta(fasta)
    assert len(proteins) == 2
    assert proteins["prot1"] == "ACDEFGHIKLMN"
    assert proteins["prot2"] == "PQRSTV"


def test_viral_peptides_from_fasta(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">E6 HPV16 E6 oncoprotein\nACDEFGHIKLMN\n")
    df = viral_peptides(fasta_path=fasta, lengths=(8,))
    assert len(df) == 5  # 12 - 8 + 1
    assert all(df["virus"] == "custom")
