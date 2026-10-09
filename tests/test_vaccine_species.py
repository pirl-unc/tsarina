"""Species synonyms select the same policy while unsupported references fail."""

import pytest

from tsarina.vaccine_construct import VaccineConfig


@pytest.mark.parametrize("name", ["human", "Human", "Homo sapiens", "homo_sapiens"])
def test_human_species_config_aliases(name):
    config = VaccineConfig(species=name)
    config.validate()
    assert config.species == "human"
    assert config.panel == "global54_abc"
    assert config.predictor == "mhcflurry"


@pytest.mark.parametrize(
    "name",
    [
        "canine",
        "CANINE",
        "dog",
        " dog ",
        "Canis familiaris",
        "Canis lupus familiaris",
        "canis_lupus_familiaris",
    ],
)
def test_dog_species_config_aliases(name):
    config = VaccineConfig.for_canine("osteosarcoma", species=name)
    config.validate()
    assert config.species == "canine"
    assert config.panel == "bundle"
    assert config.predictor == "frozen"


@pytest.mark.parametrize("name", ["canis lupis", "Canis lupus", "wolf", "mouse", "cat", "", None])
def test_unknown_or_unsupported_species_rejected(name):
    with pytest.raises(ValueError, match="species"):
        VaccineConfig(species=name)


def test_cli_unknown_species_explains_domestic_dog_name(tmp_path, capsys, monkeypatch):
    import sys

    from tsarina.cli import main

    out = tmp_path / "output"
    monkeypatch.setattr(
        sys, "argv", ["tsarina", "vaccine", "--species", "canis lupis", "-o", str(out)]
    )
    with pytest.raises(SystemExit) as exc:
        main()
    assert exc.value.code == 2
    assert "Canis lupus familiaris" in capsys.readouterr().err
    assert not out.exists()
