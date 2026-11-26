# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Species and reference-context resolution using PyEnsembl, OncoRef, and MHCgnomes (#217).

Keeps three biological axes strictly separate across target generation,
host-proteome self-exclusion, MS evidence lookup, and MHC prediction:

1. **host_species** — the patient/animal proteome against which tumor specificity
   (self-exclusion) is evaluated (e.g., ``Homo sapiens``, ``Mus musculus``,
   ``Canis lupus familiaris``).
2. **antigen_source_species** — the origin of candidate antigen sequences
   (host-derived CTAs/mutations, or an exogenous pathogen such as HPV16/EBV).
3. **mhc_species** — the restricting MHC species (typically equal to
   ``host_species``, or distinct for transgenic/xenograft contexts such as
   human ``HLA-A*02:01`` expressed in a mouse host).

Common/scientific species names are owned by PyEnsembl's species registry;
canonical taxon IDs for human/mouse/dog come from OncoRef; MHC nomenclature
and species parsing are delegated to MHCgnomes.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import pandas as pd

_PYENSEMBL_TO_ONCOREF_NAME: dict[str, str] = {
    "homo_sapiens": "human",
    "mus_musculus": "mouse",
    "canis_familiaris": "dog",
}

_CANONICAL_SHORT_ALIASES: dict[str, str] = {
    "canine": "dog",
    "murine": "mouse",
}

_DEFAULT_ENSEMBL_RELEASE = object()


def _default_ensembl_release_for_species(short_name: str) -> int:
    """Resolve the default Ensembl release for *short_name*.

    Human retains ``112`` for backward compatibility with existing reference
    indexes. Mouse (``115``, GRCm39) and dog (``98``, CanFam3.1) align with
    OncoRef's SHA-256-pinned Ensembl host references when registered.
    """
    if short_name in ("mouse", "dog"):
        try:
            from oncoref.species import for_species as _oncoref_for_species

            sp = _oncoref_for_species(short_name)
            for ref in sp.references.registered_references():
                ident = ref.get("identity", {})
                if ident.get("annotation_source") == "Ensembl":
                    release = str(ident.get("annotation_release", "")).strip()
                    if release.isdigit():
                        return int(release)
        except Exception:
            pass
        return 115 if short_name == "mouse" else 98
    return 112


@dataclass(frozen=True)
class SpeciesContext:
    """Immutable reference identity separating host, antigen source, and MHC taxa."""

    host_species: str
    host_name: str
    host_taxon: int | None
    pyensembl_species: str
    ensembl_release: int | None = 112
    assembly: str | None = None
    reference_source: str = "pyensembl"
    reference_key: str | None = None
    antigen_source_species: str | None = None
    antigen_source_taxon: int | None = None
    mhc_species: str | None = None
    mhc_taxon: int | None = None

    def __post_init__(self) -> None:
        if self.antigen_source_species is None:
            object.__setattr__(self, "antigen_source_species", self.host_species)
            if self.antigen_source_taxon is None:
                object.__setattr__(self, "antigen_source_taxon", self.host_taxon)
        if self.mhc_species is None:
            object.__setattr__(self, "mhc_species", self.host_species)
            if self.mhc_taxon is None:
                object.__setattr__(self, "mhc_taxon", self.host_taxon)

    @property
    def is_human_host(self) -> bool:
        return self.pyensembl_species == "homo_sapiens"

    @property
    def is_cross_species_mhc(self) -> bool:
        return bool(self.mhc_species and self.mhc_species != self.host_species)

    def to_provenance(self) -> dict[str, object]:
        """Return a JSON-serializable provenance dictionary for attrs/manifests."""
        return {
            "host_species": self.host_species,
            "host_name": self.host_name,
            "host_taxon": self.host_taxon,
            "pyensembl_species": self.pyensembl_species,
            "ensembl_release": self.ensembl_release,
            "assembly": self.assembly,
            "reference_source": self.reference_source,
            "reference_key": self.reference_key,
            "antigen_source_species": self.antigen_source_species,
            "antigen_source_taxon": self.antigen_source_taxon,
            "mhc_species": self.mhc_species,
            "mhc_taxon": self.mhc_taxon,
            "is_cross_species_mhc": self.is_cross_species_mhc,
        }


def _resolve_pyensembl_species(name: str):
    """Resolve a common or scientific species name via PyEnsembl's registry."""
    if not isinstance(name, str) or not name.strip():
        raise ValueError("Species name must be a non-empty string")
    raw = name.strip()
    lookup = _CANONICAL_SHORT_ALIASES.get(raw.casefold(), raw)

    try:
        from pyensembl.species import find_species_by_name
    except ImportError as error:
        raise ImportError(
            "Species resolution requires PyEnsembl; install tsarina[peptides] or tsarina[vaccine]"
        ) from error

    try:
        sp = find_species_by_name(lookup)
    except ValueError as error:
        raise ValueError(
            f"Unknown species {name!r}; pass a PyEnsembl-registered common or scientific name "
            "(e.g. 'human'/'Homo sapiens', 'mouse'/'Mus musculus', 'dog'/'Canis lupus familiaris')"
        ) from error

    latin_key = sp.latin_name
    scientific_name = latin_key.replace("_", " ").capitalize()
    short_name = _PYENSEMBL_TO_ONCOREF_NAME.get(latin_key, latin_key)
    taxon_id: int | None = None

    if short_name in {"human", "mouse", "dog"}:
        try:
            from oncoref.species import resolve_species as _oncoref_resolve

            onco_sp = _oncoref_resolve(short_name)
            taxon_id = int(onco_sp.taxon_id)
            scientific_name = str(onco_sp.scientific_name)
        except Exception:
            fallback = {
                "human": (9606, "Homo sapiens"),
                "mouse": (10090, "Mus musculus"),
                "dog": (9615, "Canis lupus familiaris"),
            }
            taxon_id, scientific_name = fallback[short_name]

    return sp, short_name, scientific_name, taxon_id


def resolve_species_context(
    host_species: str | SpeciesContext = "human",
    *,
    antigen_source_species: str | None = None,
    antigen_source_taxon: int | None = None,
    mhc_species: str | None = None,
    mhc_taxon: int | None = None,
    host_taxon: int | None = None,
    ensembl_release: int | object | None = _DEFAULT_ENSEMBL_RELEASE,
    assembly: str | None = None,
    reference_source: str = "pyensembl",
    reference_key: str | None = None,
) -> SpeciesContext:
    """Resolve a :class:`SpeciesContext` from common or scientific species names.

    Parameters
    ----------
    host_species
        Host organism (common or scientific name, or an existing
        :class:`SpeciesContext`).
    antigen_source_species, antigen_source_taxon
        Optional distinct antigen source organism/taxon (e.g. a viral species).
        Defaults to the host organism.
    mhc_species, mhc_taxon
        Optional distinct restricting MHC species/taxon (e.g. ``"Homo sapiens"``
        for an HLA-transgenic mouse host). Defaults to the host organism.
    host_taxon
        Optional explicit NCBI taxon ID; validated against the resolved species'
        known taxon ID when both are present.
    ensembl_release
        Ensembl release number (default 112; set ``None`` for non-Ensembl
        references such as a frozen ORFdata snapshot).
    assembly
        Genome assembly name (auto-resolved from PyEnsembl when ``None`` and
        ``ensembl_release`` is provided).
    reference_source
        Reference provider label (e.g. ``"pyensembl"`` or ``"orfdata"``).
    reference_key
        Optional explicit reference snapshot/bundle identity.
    """
    if isinstance(host_species, SpeciesContext):
        if (
            antigen_source_species is None
            and antigen_source_taxon is None
            and mhc_species is None
            and mhc_taxon is None
            and host_taxon is None
            and ensembl_release is _DEFAULT_ENSEMBL_RELEASE
            and assembly is None
            and reference_key is None
        ):
            return host_species
        return resolve_species_context(
            host_species.host_species,
            antigen_source_species=antigen_source_species or host_species.antigen_source_species,
            antigen_source_taxon=(
                antigen_source_taxon
                if antigen_source_taxon is not None
                else host_species.antigen_source_taxon
            ),
            mhc_species=mhc_species or host_species.mhc_species,
            mhc_taxon=mhc_taxon if mhc_taxon is not None else host_species.mhc_taxon,
            host_taxon=host_taxon if host_taxon is not None else host_species.host_taxon,
            ensembl_release=(
                host_species.ensembl_release
                if ensembl_release is _DEFAULT_ENSEMBL_RELEASE
                else ensembl_release
            ),
            assembly=assembly or host_species.assembly,
            reference_source=reference_source or host_species.reference_source,
            reference_key=reference_key or host_species.reference_key,
        )

    py_sp, short_name, scientific_name, resolved_taxon = _resolve_pyensembl_species(host_species)
    effective_release: int | None = (
        _default_ensembl_release_for_species(short_name)
        if ensembl_release is _DEFAULT_ENSEMBL_RELEASE
        else int(ensembl_release)
        if ensembl_release is not None
        else None
    )
    if host_taxon is not None and resolved_taxon is not None and int(host_taxon) != resolved_taxon:
        raise ValueError(
            f"Conflicting host_taxon={host_taxon} for species {scientific_name!r} "
            f"(expected taxon {resolved_taxon})"
        )
    effective_host_taxon = resolved_taxon if host_taxon is None else int(host_taxon)

    resolved_assembly = assembly
    if resolved_assembly is None and effective_release is not None:
        try:
            resolved_assembly = py_sp.which_reference(int(effective_release))
        except Exception:
            resolved_assembly = getattr(py_sp, "current_reference_name", None)

    if mhc_species is not None:
        _, _, resolved_mhc_species, resolved_mhc_taxon = _resolve_pyensembl_species(mhc_species)
        if (
            mhc_taxon is not None
            and resolved_mhc_taxon is not None
            and int(mhc_taxon) != resolved_mhc_taxon
        ):
            raise ValueError(
                f"Conflicting mhc_taxon={mhc_taxon} for MHC species {resolved_mhc_species!r} "
                f"(expected taxon {resolved_mhc_taxon})"
            )
        effective_mhc_taxon = resolved_mhc_taxon if mhc_taxon is None else int(mhc_taxon)
    else:
        resolved_mhc_species = scientific_name
        effective_mhc_taxon = effective_host_taxon if mhc_taxon is None else int(mhc_taxon)

    if antigen_source_species is not None:
        resolved_antigen_species = str(antigen_source_species).strip()
        if not resolved_antigen_species:
            raise ValueError("antigen_source_species must be non-empty when specified")
        effective_antigen_taxon = antigen_source_taxon
    else:
        resolved_antigen_species = scientific_name
        effective_antigen_taxon = (
            effective_host_taxon if antigen_source_taxon is None else int(antigen_source_taxon)
        )

    return SpeciesContext(
        host_species=scientific_name,
        host_name=short_name,
        host_taxon=effective_host_taxon,
        pyensembl_species=py_sp.latin_name,
        ensembl_release=effective_release,
        assembly=resolved_assembly,
        reference_source=reference_source,
        reference_key=reference_key,
        antigen_source_species=resolved_antigen_species,
        antigen_source_taxon=effective_antigen_taxon,
        mhc_species=resolved_mhc_species,
        mhc_taxon=effective_mhc_taxon,
    )


def _mhc_species_matches(parsed_species, expected_mhc_species: str) -> bool:
    """Check whether an MHCgnomes Species matches an expected canonical scientific name.

    Handles ``mhcgnomes``'s genus-level ``Canis sp.`` (``DLA``) nomenclature
    for domestic dog (``Canis lupus familiaris`` / ``Canis familiaris``) without
    aliasing non-dog species to dog.
    """
    if parsed_species is None:
        return False
    expected = expected_mhc_species.strip().casefold()
    latin = str(getattr(parsed_species, "latin_name", "") or "").strip().casefold()
    prefix = str(getattr(parsed_species, "mhc_prefix", "") or "").strip().upper()
    if latin == expected:
        return True
    if expected in {"canis lupus familiaris", "canis familiaris", "dog", "canine"}:
        return prefix == "DLA" or latin in {
            "canis familiaris",
            "canis lupus familiaris",
            "canis sp.",
        }
    if expected in {"mus musculus", "mouse"}:
        return prefix in {"H2", "H-2"} or latin == "mus musculus"
    if expected in {"homo sapiens", "human"}:
        return prefix == "HLA" or latin == "homo sapiens"
    return False


def validate_alleles_for_species(
    alleles: Iterable[str],
    mhc_species: str,
    *,
    allow_cross_species: bool = False,
) -> list[dict[str, object]]:
    """Parse MHC alleles via MHCgnomes and validate species consistency.

    Raises :class:`ValueError` if an allele cannot be parsed or belongs to a
    different MHC species than *mhc_species* (unless ``allow_cross_species=True``).
    """
    import mhcgnomes

    rows: list[dict[str, object]] = []
    for raw in alleles:
        token = str(raw).strip()
        if not token:
            continue
        parsed = mhcgnomes.parse(token, raise_on_error=False)
        if parsed is None:
            raise ValueError(f"Could not parse MHC allele {token!r} via mhcgnomes")
        parsed_species = getattr(parsed, "species", None)
        matches = _mhc_species_matches(parsed_species, mhc_species)
        if not matches and not allow_cross_species:
            observed_latin = getattr(parsed_species, "latin_name", "unknown")
            raise ValueError(
                f"MHC allele {token!r} ({observed_latin}) disagrees with expected "
                f"mhc_species={mhc_species!r}; set mhc_species explicitly for "
                "transgenic/cross-species contexts"
            )
        rows.append(
            {
                "allele": token,
                "normalized_allele": parsed.to_string(),
                "mhc_object_type": type(parsed).__name__,
                "observed_mhc_species": getattr(parsed_species, "latin_name", ""),
                "mhc_prefix": getattr(parsed_species, "mhc_prefix", ""),
                "matches_expected_species": matches,
            }
        )
    return rows


def allele_prediction_capability(
    alleles: Iterable[str],
    *,
    predictor: str = "mhcflurry",
    mhc_species: str = "Homo sapiens",
    capability_table: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Report per-allele prediction support and calibration without assuming cross-species validity.

    Parameters
    ----------
    alleles
        MHC allele names to audit.
    predictor
        Predictor backend name (e.g. ``"mhcflurry"``, ``"netmhcpan"``, ``"frozen"``).
    mhc_species
        Expected restricting MHC species scientific name.
    capability_table
        Optional reviewed capability DataFrame (with columns ``allele`` and
        ``tier``, plus optional ``reason`` / ``calibrated``) such as from a
        species-scoped frozen bundle.

    Returns
    -------
    pd.DataFrame
        Columns: ``allele``, ``normalized_allele``, ``mhc_species``,
        ``predictor``, ``capability_tier``, ``calibrated``, ``notes``.
    """
    import mhcgnomes

    columns = [
        "allele",
        "normalized_allele",
        "mhc_species",
        "predictor",
        "capability_tier",
        "calibrated",
        "notes",
    ]
    allele_list = [str(a).strip() for a in alleles if str(a).strip()]
    if not allele_list:
        return pd.DataFrame(columns=columns)

    reviewed_map: dict[str, dict] = {}
    if capability_table is not None and not capability_table.empty:
        for rec in capability_table.to_dict("records"):
            raw_a = str(rec.get("allele", "")).strip()
            if not raw_a:
                continue
            parsed = mhcgnomes.parse(raw_a, raise_on_error=False)
            key = parsed.to_string() if parsed is not None else raw_a
            reviewed_map[key] = rec
            reviewed_map[raw_a] = rec

    mhcflurry_supported: set[str] | None = None
    if predictor == "mhcflurry" and not reviewed_map:
        try:
            from mhcflurry import Class1PresentationPredictor

            model = Class1PresentationPredictor.load()
            mhcflurry_supported = {
                str(a) for a in getattr(model.affinity_predictor, "supported_alleles", ())
            }
        except Exception:
            mhcflurry_supported = None

    rows = []
    for raw in allele_list:
        parsed = mhcgnomes.parse(raw, raise_on_error=False)
        norm = parsed.to_string() if parsed is not None else raw
        parsed_species = getattr(parsed, "species", None)
        obs_species = str(getattr(parsed_species, "latin_name", "") or mhc_species)
        prefix = str(getattr(parsed_species, "mhc_prefix", "") or "").upper()

        if norm in reviewed_map or raw in reviewed_map:
            rec = reviewed_map.get(norm, reviewed_map.get(raw, {}))
            tier = str(rec.get("tier", "unverified"))
            calibrated = bool(rec.get("calibrated", tier == "empirically_supported"))
            notes = str(rec.get("reason", rec.get("notes", "reviewed_capability_table")))
        elif parsed is None:
            tier = "unsupported"
            calibrated = False
            notes = "unparseable_mhc_allele"
        elif predictor == "frozen":
            tier = "unverified"
            calibrated = False
            notes = "no_reviewed_capability_entry"
        elif prefix == "HLA":
            if mhcflurry_supported is not None:
                if norm in mhcflurry_supported or raw in mhcflurry_supported:
                    tier = "empirically_supported"
                    calibrated = True
                    notes = "supported_by_installed_mhcflurry_model"
                else:
                    tier = "unsupported"
                    calibrated = False
                    notes = "absent_from_installed_mhcflurry_supported_alleles"
            else:
                tier = "unverified"
                calibrated = False
                notes = f"predictor_{predictor}_calibration_not_verified_in_environment"
        else:
            # Non-human MHC (e.g. H-2 mouse, DLA dog): never assume human
            # predictor calibration without explicit species-scoped capability review.
            if mhcflurry_supported is not None and (
                norm in mhcflurry_supported or raw in mhcflurry_supported
            ):
                tier = "unverified_nonhuman_model_support"
                calibrated = False
                notes = (
                    f"{obs_species} allele present in predictor vocabulary, but "
                    "non-human calibration/presentation validity requires explicit review"
                )
            else:
                tier = "unsupported"
                calibrated = False
                notes = (
                    f"no verified {predictor} capability entry for {obs_species} allele {norm}; "
                    "observed MHC restriction remains separate from prediction"
                )

        rows.append(
            {
                "allele": raw,
                "normalized_allele": norm,
                "mhc_species": obs_species,
                "predictor": predictor,
                "capability_tier": tier,
                "calibrated": calibrated,
                "notes": notes,
            }
        )
    return pd.DataFrame(rows, columns=columns)
