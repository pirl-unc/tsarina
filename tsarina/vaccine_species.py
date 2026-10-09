"""Resolve vaccine policy species using PyEnsembl's species-name registry."""


def resolve_vaccine_species(value):
    """Return a supported policy identifier without changing reference identity.

    Human/canine are the original policy labels. Other common/scientific names
    are owned by PyEnsembl, which distinguishes domestic dog from wolf and other
    genomes. The frozen evidence bundle still determines taxon and assembly.
    """
    if not isinstance(value, str) or not value.strip():
        raise ValueError("Vaccine species requires a nonempty common or scientific name")
    label = value.strip().casefold()
    if label in {"human", "canine"}:
        return label

    try:
        from pyensembl.species import dog, find_species_by_name, human
    except ImportError as error:
        raise ImportError("Species aliases require PyEnsembl; install tsarina[vaccine]") from error
    try:
        species = find_species_by_name(value)
    except ValueError as error:
        raise ValueError(
            f"Unknown vaccine species {value!r}; use human (Homo sapiens) or dog "
            "(Canis familiaris / Canis lupus familiaris). The canine policy label is also accepted."
        ) from error
    if species is human:
        return "human"
    if species is dog:
        return "canine"
    raise ValueError(
        f"Unsupported vaccine species {value!r} ({species.latin_name}); "
        "design policies currently cover human and domestic dog"
    )
