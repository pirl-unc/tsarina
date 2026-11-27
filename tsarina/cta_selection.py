"""Current OncoRef target-selection decisions projected to gene members.

Human, mouse, and canine selection all read OncoRef's unified set-based target
catalog (:mod:`oncoref.target_selection` / ``oncoref.<species>.target_selection()``,
policy ``oncoref.target-selection.*``), which classifies every protein group as
``CTA`` (with ``cta_tier``), ``Cancer Antigen only``, ``Reproductive Antigen
only``, ``Other`` (with subtype), or ``Unsafe``.
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from functools import lru_cache

import pandas as pd
from oncoref import antigen_sets, gene_identity_catalog, target_selection
from oncoref import host_context as _oncoref_host_context
from oncoref.host_context import SEX_RESTRICTION_CODES, HostContext
from oncoref.species import for_species as _oncoref_for_species

_HUMAN_SPECIES_ALIASES = frozenset({"human", "homo_sapiens", "homo sapiens", "9606"})
CTA_TIERS = target_selection.CTA_TIERS

_SEX_RESTRICTION_CODE_BY_LOWER: dict[str, str] = {
    code.lower(): code for code in SEX_RESTRICTION_CODES
}
_ECTOPIC_ACTIVATION_MODES: frozenset[str] = frozenset({"ectopic", "ectopic_and_lineage_retained"})


def _parse_json_tuple(value) -> tuple[str, ...]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ()
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return ()
        if stripped.startswith("["):
            parsed = json.loads(stripped)
            return tuple(str(item) for item in parsed if item)
        return (stripped,)
    if isinstance(value, Iterable):
        return tuple(str(item) for item in value if item)
    return (str(value),)


def _clean_text(value) -> str | None:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).strip()
    return text or None


@lru_cache(maxsize=1)
def _gene_identity_index() -> dict[str, dict]:
    """Map canonical Ensembl gene ID to symbol and synonym tuple."""
    catalog = gene_identity_catalog()
    index: dict[str, dict] = {}
    for row in catalog.itertuples(index=False):
        gene_id = str(row.canonical_gene_id)
        symbol = _clean_text(getattr(row, "canonical_symbol", None))
        synonyms = _parse_json_tuple(getattr(row, "gene_name_synonyms", None))
        index[gene_id] = {
            "gene_id": gene_id,
            "symbol": symbol,
            "synonyms": synonyms,
            "identity_status": "resolved",
        }
    return index


def _row_reasons(row) -> tuple[str, ...]:
    reasons: list[str] = []
    decision_reason = _clean_text(getattr(row, "decision_reason", None))
    if decision_reason:
        reasons.append(decision_reason)
    unsafe_reasons = _clean_text(getattr(row, "Unsafe_reasons", None))
    if unsafe_reasons:
        for part in unsafe_reasons.split(";"):
            part = part.strip()
            if part and part not in reasons:
                reasons.append(part)
    return tuple(reasons)


def _row_evidence_flags(row) -> tuple[str, ...]:
    flags: list[str] = []
    for prefix, attr in (
        ("cta_tier", "cta_tier"),
        ("rare_cell", "rare_cell_status"),
        ("nuclear_caution", "nuclear_rna_caution"),
        ("clinical_clearance", "clinical_clearance_scope"),
        ("literature", "literature_effect"),
        ("expression_bounds", "expression_bounds_status"),
    ):
        value = _clean_text(getattr(row, attr, None))
        if value and value not in {"not_extracted", "not_reviewed", "none", "unassessed"}:
            flags.append(f"{prefix}:{value}")
    return tuple(flags)


@lru_cache(maxsize=1)
def _human_target_selection() -> pd.DataFrame:
    """Normalized view of OncoRef's rebuilt human target_selection catalog."""
    raw = target_selection.results().copy()
    meta = target_selection.metadata()
    identities_by_id = _gene_identity_index()
    reference_id = f"homo_sapiens:GRCh38:Ensembl:{antigen_sets.ENSEMBL_RELEASE}"
    snapshot_id = f"sha256:{meta['data_sha256']}"

    gene_ids_col: list[tuple[str, ...]] = []
    gene_identities_col: list[tuple[dict, ...]] = []
    identity_status_col: list[str] = []
    reasons_col: list[tuple[str, ...]] = []
    flags_col: list[tuple[str, ...]] = []

    for row in raw.itertuples(index=False):
        members = _parse_json_tuple(row.member_gene_ids)
        identities = tuple(
            identities_by_id.get(
                gid,
                {
                    "gene_id": gid,
                    "symbol": None,
                    "synonyms": (),
                    "identity_status": "registry_member_outside_canonical_catalog",
                },
            )
            for gid in members
        )
        all_resolved = all(item["identity_status"] == "resolved" for item in identities)
        gene_ids_col.append(members)
        gene_identities_col.append(identities)
        identity_status_col.append(
            "resolved" if all_resolved else "registry_member_outside_canonical_catalog"
        )
        reasons_col.append(_row_reasons(row))
        flags_col.append(_row_evidence_flags(row))

    raw["species"] = "human"
    raw["entity_kind"] = "protein_group"
    raw["label"] = raw["protein_group_label"]
    raw["gene_ids"] = gene_ids_col
    raw["gene_identities"] = gene_identities_col
    raw["identity_status"] = identity_status_col
    raw["policy_id"] = raw["policy"]
    raw["reference_id"] = reference_id
    raw["source_snapshot_id"] = snapshot_id
    raw["reasons"] = reasons_col
    raw["evidence_flags"] = flags_col
    return raw


@lru_cache(maxsize=8)
def _nonhuman_target_selection(species_key: str) -> pd.DataFrame:
    """Normalized view of OncoRef's cross-species target_selection catalog."""
    sp = _oncoref_for_species(species_key)
    raw = sp.antigens.target_selection().copy()
    coverage = sp.antigens.target_selection_coverage()
    snapshot_id = str(sp.antigens.cta_snapshot().as_dict().get("snapshot_id", ""))
    default_ref_id = str(coverage.get("reference", {}).get("reference_id", ""))

    gene_ids_col: list[tuple[str, ...]] = []
    identity_status_col: list[str] = []
    reasons_col: list[tuple[str, ...]] = []
    flags_col: list[tuple[str, ...]] = []

    for row in raw.itertuples(index=False):
        members = _parse_json_tuple(getattr(row, "member_gene_ids", ()))
        gene_ids_col.append(members)
        complete = bool(getattr(row, "identity_translation_pass", False)) and bool(members)
        identity_status_col.append("resolved" if complete else "unmapped_or_incomplete")
        reasons_col.append(_row_reasons(row))
        flags_col.append(_row_evidence_flags(row))

    raw["species"] = sp.species.name
    raw["entity_kind"] = "protein_group"
    raw["label"] = raw["protein_group_label"]
    raw["gene_ids"] = gene_ids_col
    raw["identity_status"] = identity_status_col
    policy_col = "policy" if "policy" in raw.columns else "target_policy_id"
    raw["policy_id"] = raw[policy_col] if policy_col in raw.columns else coverage.get("policy_id")
    if "reference_id" not in raw.columns:
        raw["reference_id"] = default_ref_id
    raw["source_snapshot_id"] = snapshot_id
    raw["reasons"] = reasons_col
    raw["evidence_flags"] = flags_col
    return raw


def CTA_selection_coverage(*, species: str = "human") -> dict:
    """Return OncoRef's target-selection readiness, anatomy, and category counts."""
    return dict(_oncoref_for_species(species).antigens.target_selection_coverage())


_SEX_ELIGIBILITY_MODES = ("all_adult", "sex_organ_only")


def _normalize_patient_sex(sex: str) -> str:
    token = str(sex).strip().lower()
    if token not in ("female", "male"):
        raise ValueError(f"Unsupported patient sex {sex!r}; expected 'female' or 'male'.")
    return token


def _normalize_sex_eligibility(eligibility: str) -> str:
    token = str(eligibility).strip().lower()
    if token not in _SEX_ELIGIBILITY_MODES:
        raise ValueError(
            f"Unsupported sex eligibility {eligibility!r}; expected one of {_SEX_ELIGIBILITY_MODES}."
        )
    return token


def _normalize_sex_codes(sex_code: str | Iterable[str]) -> tuple[str, ...]:
    """Normalize one or more 4-way OncoRef sex-restriction codes (case-insensitive)."""
    raw_items: list[str] = []
    if isinstance(sex_code, str):
        raw_items.extend(part.strip() for part in sex_code.split(",") if part.strip())
    else:
        for item in sex_code:
            for part in str(item).split(","):
                token = part.strip()
                if token:
                    raw_items.append(token)
    if not raw_items:
        raise ValueError(f"sex_code must specify at least one code from {SEX_RESTRICTION_CODES}.")
    normalized: list[str] = []
    for item in raw_items:
        canonical = _SEX_RESTRICTION_CODE_BY_LOWER.get(item.lower())
        if canonical is None:
            raise ValueError(
                f"Unsupported sex_restriction_code {item!r}; expected one of {SEX_RESTRICTION_CODES}."
            )
        if canonical not in normalized:
            normalized.append(canonical)
    return tuple(normalized)


def host_compartments(species: str = "human") -> tuple[str, ...]:
    """Return canonical OncoRef HostContext reproductive compartments for ``species``."""
    return tuple(_oncoref_host_context.compartments(species=species))


def build_host_context(
    *,
    patient_sex: str | None = None,
    pregnant: bool = False,
    lactating: bool = False,
    protect_tissues: Iterable[str] = (),
    accept_loss_tissues: Iterable[str] = (),
    absent_tissues: Iterable[str] = (),
    allow_incomplete_scope: bool = False,
    species: str = "human",
) -> HostContext | None:
    """Build a composable OncoRef :class:`HostContext` for human, mouse, or dog.

    Translates high-level biological state and compartment choices directly via
    :func:`oncoref.host_context.for_host`. Returns ``None`` when no host
    constraints are requested.
    """
    profile = "all"
    if patient_sex is not None:
        norm_sex = _normalize_patient_sex(patient_sex)
        profile = "female-patient" if norm_sex == "female" else "male-patient"

    def _flat(values: Iterable[str]) -> tuple[str, ...]:
        if isinstance(values, str):
            values = [values]
        out: list[str] = []
        for item in values:
            for part in str(item).split(","):
                token = part.strip()
                if token:
                    out.append(token)
        return tuple(out)

    return _oncoref_host_context.for_host(
        profile=profile,
        pregnant=bool(pregnant),
        lactating=bool(lactating),
        allow_tissues=_flat(accept_loss_tissues),
        protect_tissues=_flat(protect_tissues),
        absent_tissues=_flat(absent_tissues),
        allow_incomplete_scope=bool(allow_incomplete_scope),
        species=species,
    )


@lru_cache(maxsize=8)
def _sex_restricted_table(species_key: str) -> pd.DataFrame:
    if species_key in _HUMAN_SPECIES_ALIASES:
        raw = _oncoref_for_species("human").antigens.sex_restricted_antigens().copy()
        sel_full = _human_target_selection()
    else:
        sp = _oncoref_for_species(species_key)
        cov = sp.antigens.target_selection_coverage()
        status = cov.get("selection_status")
        if status != "evaluated":
            raise ValueError(
                f"Sex-restricted antigen eligibility is not yet available for species "
                f"{sp.species.name!r} (selection_status={status!r})."
            )
        raw = sp.antigens.sex_restricted_antigens().copy()
        sel_full = _nonhuman_target_selection(species_key)

    extra_cols = [
        c
        for c in (
            "label",
            "gene_ids",
            "category",
            "cta_tier",
            "tumor_activation_mode",
            "tumor_support",
        )
        if c in sel_full.columns and c not in raw.columns
    ]
    if not extra_cols:
        return raw
    return raw.merge(sel_full[["entity_id", *extra_cols]], on="entity_id", how="left")


def CTA_sex_restricted_antigens(
    *,
    sex: str | None = None,
    eligibility: str = "all_adult",
    include_reproductive_only: bool = True,
    sex_code: str | Iterable[str] | None = None,
    species: str = "human",
) -> pd.DataFrame:
    """Return OncoRef's sex-restricted antigen catalog, optionally filtered for a patient sex.

    Works for any OncoRef-evaluated species (``"human"``, ``"mouse"``, ``"dog"``).

    Parameters
    ----------
    sex
        ``"female"``, ``"male"``, or ``None`` (return all classified
        reproductive/CTA entities with both female and male eligibility flags).
    eligibility
        ``"all_adult"`` (default — zero normal expression in any adult tissue of
        that sex, e.g. male-restricted & pregnancy/peripartum-only for female
        patients) or ``"sex_organ_only"`` (permits non-vital sex-organ expression).
    include_reproductive_only
        Include ``Reproductive Antigen only`` entities alongside ``CTA`` entities
        (default ``True``).
    sex_code
        Optional 4-way OncoRef ``sex_restriction_code`` filter (single code or
        iterable/comma-separated codes from :data:`SEX_RESTRICTION_CODES`, e.g.
        ``"CTA-F"``, ``("CTA-M", "CTA-P")``, or ``"RA-M-only"``).
    species
        Species name resolved via OncoRef (default ``"human"``; also ``"mouse"``, ``"dog"``).
    """
    elig = _normalize_sex_eligibility(eligibility)
    species_key = str(species).strip().lower()
    df = _sex_restricted_table(species_key)
    if not include_reproductive_only:
        df = df.loc[df["is_cta"].astype(bool)]
    if sex is not None:
        norm_sex = _normalize_patient_sex(sex)
        col = f"{norm_sex}_{elig}_eligible"
        df = df.loc[df[col].astype(bool)]
    if sex_code is not None:
        codes = set(_normalize_sex_codes(sex_code))
        if "sex_restriction_code" in df.columns:
            df = df.loc[df["sex_restriction_code"].isin(codes)]
        else:
            df = df.iloc[:0]
    return df.copy().reset_index(drop=True)


def CTA_sex_eligible_gene_ids(
    sex: str,
    *,
    eligibility: str = "all_adult",
    include_reproductive_only: bool = True,
    sex_code: str | Iterable[str] | None = None,
    species: str = "human",
) -> set[str]:
    """Return Ensembl gene IDs eligible for a patient of ``sex`` under ``eligibility``."""
    df = CTA_sex_restricted_antigens(
        sex=sex,
        eligibility=eligibility,
        include_reproductive_only=include_reproductive_only,
        sex_code=sex_code,
        species=species,
    )
    if str(species).strip().lower() in _HUMAN_SPECIES_ALIASES:
        entity_ids = set(df["entity_id"])
        frame = _selection_gene_table()
        return set(frame.loc[frame["selection_entity_id"].isin(entity_ids), "Ensembl_Gene_ID"])
    gene_ids: set[str] = set()
    for members in df.get("gene_ids", ()):
        if isinstance(members, (tuple, list)):
            gene_ids.update(str(g) for g in members if g)
    return gene_ids


def CTA_sex_eligible_gene_names(
    sex: str,
    *,
    eligibility: str = "all_adult",
    include_reproductive_only: bool = True,
    sex_code: str | Iterable[str] | None = None,
    species: str = "human",
) -> set[str]:
    """Return gene symbols eligible for a patient of ``sex`` under ``eligibility``."""
    df = CTA_sex_restricted_antigens(
        sex=sex,
        eligibility=eligibility,
        include_reproductive_only=include_reproductive_only,
        sex_code=sex_code,
        species=species,
    )
    if str(species).strip().lower() in _HUMAN_SPECIES_ALIASES:
        entity_ids = set(df["entity_id"])
        frame = _selection_gene_table()
        return set(frame.loc[frame["selection_entity_id"].isin(entity_ids), "Symbol"])
    names: set[str] = set()
    for label in df.get("label", ()):
        text = _clean_text(label)
        if text:
            for part in text.split("/"):
                part = part.strip()
                if part:
                    names.add(part)
    return names


def CTA_selection(
    category: str | Iterable[str] | None = None,
    *,
    species: str = "human",
    sex_code: str | Iterable[str] | None = None,
    ectopic_only: bool = False,
    replicated_only: bool = False,
    host_context: HostContext | None = None,
    filter_host: bool = False,
):
    """Return detached current OncoRef selection decisions.

    Uses OncoRef's unified cross-species :mod:`oncoref.target_selection` catalog
    across human, mouse, and canine, enriched with normalized gene-ID tuples and
    provenance receipts.

    Parameters
    ----------
    category
        Optional OncoRef category or iterable of categories (e.g. ``"CTA"``,
        ``("CTA", "Reproductive Antigen only")``).
    species
        Species name resolved via OncoRef (default ``"human"``; also ``"mouse"``, ``"dog"``).
    sex_code
        Optional 4-way ``sex_restriction_code`` filter (e.g. ``"CTA-F"``,
        ``("CTA-M", "CTA-P")``, ``"RA-M-only"``).
    ectopic_only
        When ``True``, keep only rows with ectopic tumor activation
        (``tumor_activation_mode`` in ``{"ectopic", "ectopic_and_lineage_retained"}``),
        excluding purely ``lineage_retained`` targets.
    replicated_only
        When ``True``, keep only rows with multi-patient replicated tumor support
        (``tumor_support == "replicated"``), excluding ``single_sample`` targets.
    host_context
        Optional :class:`oncoref.host_context.HostContext` to evaluate via
        :func:`oncoref.host_context.assess`, attaching ``host_status``,
        ``host_protected_max_tpm``, ``host_limiting_compartment``, ``host_reasons``,
        ``host_unknown_compartments``, and ``pass_host`` columns.
    filter_host
        When ``True`` and ``host_context`` is supplied, keep only rows with
        ``pass_host == True``.
    """
    key = str(species).strip().lower()
    if key in _HUMAN_SPECIES_ALIASES:
        frame = _human_target_selection()
    else:
        frame = _nonhuman_target_selection(key)

    if category is not None:
        if isinstance(category, str):
            frame = frame.loc[frame["category"].eq(category)]
        else:
            frame = frame.loc[frame["category"].isin(set(category))]
    if sex_code is not None:
        codes = set(_normalize_sex_codes(sex_code))
        if "sex_restriction_code" in frame.columns:
            frame = frame.loc[frame["sex_restriction_code"].isin(codes)]
        else:
            frame = frame.iloc[:0]
    if ectopic_only:
        if "tumor_activation_mode" in frame.columns:
            frame = frame.loc[frame["tumor_activation_mode"].isin(_ECTOPIC_ACTIVATION_MODES)]
        else:
            frame = frame.iloc[:0]
    if replicated_only:
        if "tumor_support" in frame.columns:
            frame = frame.loc[frame["tumor_support"].eq("replicated")]
        else:
            frame = frame.iloc[:0]

    out = frame.copy().reset_index(drop=True)
    if host_context is not None:
        if not isinstance(host_context, HostContext):
            raise ValueError("host_context must be an oncoref.host_context.HostContext instance")
        resolved_species = _oncoref_for_species(key).species.name
        if host_context.species != resolved_species:
            raise ValueError(
                f"host_context.species={host_context.species!r} does not match "
                f"species={resolved_species!r}"
            )
        assessment = _oncoref_host_context.assess(out, host_context)
        for col in assessment.columns:
            out[col] = assessment[col].to_numpy()
        if filter_host:
            out = out.loc[out["pass_host"].astype(bool)].reset_index(drop=True)
        out.attrs["host_context"] = host_context.as_dict()
    return out


@lru_cache(maxsize=1)
def _selection_gene_table() -> pd.DataFrame:
    rows = []
    for decision in _human_target_selection().to_dict("records"):
        identities = {item["gene_id"]: item for item in decision["gene_identities"]}
        cleared_organs = _parse_json_tuple(decision.get("clinical_cleared_organs"))
        cleared_epitopes = _parse_json_tuple(decision.get("clinical_cleared_epitopes"))
        cleared_alleles = _parse_json_tuple(decision.get("clinical_hla_alleles"))
        cleared_transcripts = tuple(
            t.split(".")[0] for t in _parse_json_tuple(decision.get("clinical_cleared_transcripts"))
        )
        unverified_isoforms = tuple(
            t.split(".")[0] for t in _parse_json_tuple(decision.get("clinical_unverified_isoforms"))
        )
        rep_transcripts = tuple(
            t.split(".")[0] for t in _parse_json_tuple(decision.get("representative_transcripts"))
        )
        qualifying_cancers = _parse_json_tuple(decision.get("qualifying_cancer_codes"))
        for gene_id in decision["gene_ids"]:
            identity = identities[gene_id]
            rows.append(
                {
                    "Ensembl_Gene_ID": gene_id,
                    "Symbol": identity["symbol"] or gene_id,
                    "selection_synonyms": tuple(identity["synonyms"] or ()),
                    **{
                        "selection_" + key: decision[key]
                        for key in (
                            "entity_id",
                            "entity_kind",
                            "label",
                            "gene_ids",
                            "category",
                            "subtype",
                            "classification",
                            "cta_tier",
                            "is_cancer_antigen",
                            "is_reproductive_antigen",
                            "is_target_candidate",
                            "assessment_status",
                            "identity_status",
                            "policy_id",
                            "reference_id",
                            "reasons",
                            "evidence_flags",
                            "source_snapshot_id",
                        )
                    },
                    "selection_sex_restriction_class": _clean_text(
                        decision.get("sex_restriction_class")
                    ),
                    "selection_sex_restriction_code": _clean_text(
                        decision.get("sex_restriction_code")
                    ),
                    "selection_tumor_activation_mode": _clean_text(
                        decision.get("tumor_activation_mode")
                    ),
                    "selection_tumor_support": _clean_text(decision.get("tumor_support")),
                    "selection_qualifying_cancer_codes": qualifying_cancers,
                    "selection_rare_cell_status": _clean_text(decision.get("rare_cell_status")),
                    "selection_nuclear_rna_caution": _clean_text(
                        decision.get("nuclear_rna_caution")
                    ),
                    "selection_clinical_clearance_scope": _clean_text(
                        decision.get("clinical_clearance_scope")
                    ),
                    "selection_clinical_cleared_organs": cleared_organs,
                    "selection_clinical_cleared_epitopes": cleared_epitopes,
                    "selection_clinical_hla_alleles": cleared_alleles,
                    "selection_clinical_cleared_transcripts": cleared_transcripts,
                    "selection_clinical_unverified_isoforms": unverified_isoforms,
                    "selection_representative_transcripts": rep_transcripts,
                    "selection_exact_unique_9mers": decision.get("exact_unique_9mers"),
                    "selection_safe_2mm_9mers": decision.get("safe_2mm_9mers"),
                    "selection_unsafe_near_self_1mm_pct": decision.get("unsafe_near_self_1mm_pct"),
                }
            )
    frame = pd.DataFrame(rows)
    if frame.Ensembl_Gene_ID.duplicated().any():
        raise ValueError("OncoRef current selection assigns a gene to multiple entities")
    frame["is_current_cta"] = frame.selection_category.eq("CTA").fillna(False).astype(bool)
    for col in (
        "selection_is_cancer_antigen",
        "selection_is_reproductive_antigen",
        "selection_is_target_candidate",
    ):
        frame[col] = frame[col].fillna(False).astype(bool)
    return frame


@lru_cache(maxsize=8)
def _selection_representative_transcripts_for_species(
    species_key: str,
) -> dict[str, tuple[str, ...]]:
    if species_key in _HUMAN_SPECIES_ALIASES:
        frame = _selection_gene_table()
        mapping: dict[str, tuple[str, ...]] = {}
        for row in frame.itertuples(index=False):
            preferred: list[str] = []
            for tid in getattr(row, "selection_clinical_cleared_transcripts", ()) or ():
                if tid and tid not in preferred:
                    preferred.append(tid)
            for tid in getattr(row, "selection_representative_transcripts", ()) or ():
                if tid and tid not in preferred:
                    preferred.append(tid)
            if preferred:
                mapping[row.Ensembl_Gene_ID] = tuple(preferred)
        return mapping

    raw = _nonhuman_target_selection(species_key)
    mapping = {}
    for row in raw.itertuples(index=False):
        members = getattr(row, "gene_ids", ()) or ()
        rep_tx = tuple(
            t.split(".")[0]
            for t in _parse_json_tuple(getattr(row, "representative_transcripts", None))
            if t
        )
        if not members or not rep_tx:
            continue
        for gid in members:
            if gid:
                mapping[str(gid)] = rep_tx
    return mapping


def selection_representative_transcripts(species: str = "human") -> dict[str, tuple[str, ...]]:
    """Map Ensembl gene ID to unversioned OncoRef representative transcripts.

    Supports human, mouse, and dog catalogs from OncoRef >=1.55.0. On human,
    clinically cleared transcripts (e.g. ``CTAG2`` LAGE-1a ``ENST00000369585``
    carrying ``SLLMWITQC``) take precedence ahead of general representative
    transcripts so downstream peptide enumeration avoids unverified longer
    isoforms (such as ``CTAG2`` LAGE-1b ``ENST00000247306``).
    """
    key = str(species).strip().lower()
    if key not in _HUMAN_SPECIES_ALIASES:
        key = _oncoref_for_species(key).species.name
    return _selection_representative_transcripts_for_species(key)


def selection_gene_table() -> pd.DataFrame:
    """Detached gene rows retaining the current protein decision and provenance."""
    return _selection_gene_table().copy()


def _filter_selection_frame(
    category: str | Iterable[str] | None = "CTA",
    *,
    cta_tier: str | Iterable[str] | None = None,
) -> pd.DataFrame:
    frame = _selection_gene_table()
    if category is not None:
        if isinstance(category, str):
            frame = frame.loc[frame.selection_category.eq(category).fillna(False)]
        else:
            frame = frame.loc[frame.selection_category.isin(set(category))]
    if cta_tier is not None:
        tiers = {cta_tier} if isinstance(cta_tier, str) else set(cta_tier)
        unknown = tiers - set(CTA_TIERS)
        if unknown:
            raise ValueError(
                f"Unknown CTA tier(s): {sorted(unknown)}; expected subset of {CTA_TIERS}"
            )
        frame = frame.loc[frame.selection_cta_tier.isin(tiers)]
    return frame


def selection_gene_ids(
    category: str | Iterable[str] | None = "CTA",
    *,
    cta_tier: str | Iterable[str] | None = None,
) -> set[str]:
    return set(_filter_selection_frame(category=category, cta_tier=cta_tier).Ensembl_Gene_ID)


def selection_gene_names(
    category: str | Iterable[str] | None = "CTA",
    *,
    cta_tier: str | Iterable[str] | None = None,
) -> set[str]:
    return set(_filter_selection_frame(category=category, cta_tier=cta_tier).Symbol)


def selection_flag_gene_ids(flag_column: str) -> set[str]:
    frame = _selection_gene_table()
    return set(frame.loc[frame[flag_column].astype(bool), "Ensembl_Gene_ID"])


def selection_flag_gene_names(flag_column: str) -> set[str]:
    frame = _selection_gene_table()
    return set(frame.loc[frame[flag_column].astype(bool), "Symbol"])


def selection_proteoform_groups() -> dict[str, tuple[str, ...]]:
    """Selected reference protein groups with every member gene symbol."""
    frame = _selection_gene_table()
    groups = {}
    for _, members in frame.loc[frame.is_current_cta].groupby("selection_entity_id", sort=False):
        symbols = tuple(sorted(members.Symbol))
        if len(symbols) > 1:
            groups["/".join(symbols)] = symbols
    return groups


@lru_cache(maxsize=1)
def _identity_aliases():
    symbols, synonyms = {}, {}
    for row in _selection_gene_table().itertuples(index=False):
        symbols[row.Symbol.casefold()] = row.Symbol
        for synonym in row.selection_synonyms:
            synonyms.setdefault(synonym.casefold(), set()).add(row.Symbol)
    return symbols, {
        key: next(iter(values)) for key, values in synonyms.items() if len(values) == 1
    }


def selection_symbol_for_alias(name):
    symbols, synonyms = _identity_aliases()
    token = str(name).strip().casefold()
    return symbols.get(token) or synonyms.get(token)


def selection_provenance(*, ensembl_release: int | None = None) -> dict:
    """Selection receipts, independent of a consumer's sequence reference."""
    frame = _selection_gene_table()
    meta = target_selection.metadata()
    return {
        "cta_selection_policy_ids": sorted(frame.selection_policy_id.dropna().unique()),
        "cta_selection_reference_ids": sorted(frame.selection_reference_id.dropna().unique()),
        "cta_selection_source_snapshot_ids": sorted(
            frame.selection_source_snapshot_id.dropna().unique()
        ),
        "cta_selection_release_ready": bool(meta.get("release_ready", False)),
        "cta_selection_open_reviews": list(meta.get("open_reviews", ())),
        "background_policy": "tsarina.current-selection-alternate-haplotypes.v1",
        **({"ensembl_release": ensembl_release} if ensembl_release is not None else {}),
    }
