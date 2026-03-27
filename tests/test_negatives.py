from perseus.negatives import healthy_tissue_peptide_details, healthy_tissue_peptides


def test_healthy_tissue_peptides_no_data():
    result = healthy_tissue_peptides(iedb_path=None, cedar_path=None)
    assert isinstance(result, set)
    assert len(result) == 0


def test_healthy_tissue_peptides_nonexistent_path():
    result = healthy_tissue_peptides(iedb_path="/nonexistent/path.csv")
    assert isinstance(result, set)
    assert len(result) == 0


def test_healthy_tissue_details_no_data():
    df = healthy_tissue_peptide_details(iedb_path=None)
    assert len(df) == 0
    assert "peptide" in df.columns
    assert "in_somatic" in df.columns
    assert "in_thymus" in df.columns


def test_invalid_tissue_type():
    import pytest

    with pytest.raises(ValueError, match="tissue_type"):
        healthy_tissue_peptides(tissue_type="invalid")
