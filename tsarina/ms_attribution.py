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

"""Reusable restriction-aware MS evidence attribution across target families (#216).

Provides one provenance-preserving evidence contract shared across CTA, viral,
recurrent-mutation, personalized-target, and downstream vaccine workflows
(e.g. Vaxrank).  Distinguishes:

1. ``exact_requested_mhc`` — positive MS observation with measured exact /
   monoallelic restriction matching a requested molecular MHC allele.
2. ``donor_set_ambiguous`` — positive MS observation from a multi-allelic
   sample or donor allele set containing a requested MHC allele (ambiguous
   which allele presented the peptide).
3. ``predicted_restriction`` — positive MS observation without measured exact
   or donor-set match to the requested allele, paired with an explicit
   computational prediction assignment (never conflated with measured MS).
4. ``class_or_serotype_only`` — positive MS observation reporting only MHC
   class or serotype (preserved without inventing a 4-digit allele claim).
5. ``any_mhc_peptide`` — positive MS observation of the **exact peptide** with
   any/unknown/other-allele restriction (never transferred from a longer
   nested ligand).
6. ``no_ms_match`` — evidence was loaded/queried and no positive exact-peptide
   MS observation matched.
7. ``evidence_unavailable`` — MS evidence was not loaded or skipped.
"""

from __future__ import annotations

from collections.abc import Iterable
from functools import lru_cache

import pandas as pd

from .ms_evidence import filter_positive_ms_observations
from .species import SpeciesContext, resolve_species_context

EXACT_REQUESTED_MHC = "exact_requested_mhc"
DONOR_SET_AMBIGUOUS = "donor_set_ambiguous"
PREDICTED_RESTRICTION = "predicted_restriction"
CLASS_OR_SEROTYPE_ONLY = "class_or_serotype_only"
ANY_MHC_PEPTIDE = "any_mhc_peptide"
NO_MS_MATCH = "no_ms_match"
EVIDENCE_UNAVAILABLE = "evidence_unavailable"

ATTRIBUTION_TIER_ORDER: tuple[str, ...] = (
    EXACT_REQUESTED_MHC,
    DONOR_SET_AMBIGUOUS,
    PREDICTED_RESTRICTION,
    CLASS_OR_SEROTYPE_ONLY,
    ANY_MHC_PEPTIDE,
    NO_MS_MATCH,
    EVIDENCE_UNAVAILABLE,
)
_TIER_RANK: dict[str, int] = {tier: idx for idx, tier in enumerate(ATTRIBUTION_TIER_ORDER)}

VALID_EVIDENCE_SCOPES: frozenset[str] = frozenset({"all", "requested_mhc", "any_mhc"})

_MULTIALLELIC_PROVENANCES: frozenset[str] = frozenset(
    {
        "sample_allele_match",
        "pmid_class_pool",
        "donor_set",
        "multiallelic",
        "multiallelic_sample",
        "pooled_sample",
    }
)

_ATTRIBUTION_COLUMNS: tuple[str, ...] = (
    "ms_evidence_status",
    "ms_attribution_tier",
    "matched_requested_alleles",
    "mhc_restriction_raw",
    "mhc_restriction_normalized",
    "mhc_resolution",
    "mhc_allele_provenance",
    "donor_set_alleles",
    "observed_mhc_species",
    "host_species",
    "antigen_source_species",
    "mhc_species",
    "assay_modality",
    "ms_outcome_status",
    "is_ms_observation",
    "source_tissue",
    "disease",
    "src_cancer",
    "src_healthy_tissue",
    "src_healthy_reproductive",
    "src_healthy_reproductive_female",
    "src_healthy_reproductive_male",
    "src_cell_line",
    "src_ex_vivo",
    "pmid",
    "assay_iri",
    "provenance_id",
    "source",
    "dataset_version",
)


def _is_truthy(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return False
    return str(value).strip().lower() in {"true", "1", "1.0", "yes"}


def _clean_str(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


@lru_cache(maxsize=4096)
def parse_mhc_token(token: str) -> tuple[str, str, str, str]:
    """Parse one MHC label via MHCgnomes.

    Returns ``(normalized, resolution, species_latin, mhc_prefix)``.
    ``resolution`` is one of ``exact_allele``, ``gene``, ``serotype``,
    ``class_only``, ``haplotype``, or ``unknown``.
    """
    import mhcgnomes

    cleaned = token.strip()
    if not cleaned or cleaned.lower() in {"unknown", "none", "nan", "n/a"}:
        return "", "unknown", "", ""

    parsed = mhcgnomes.parse(cleaned, raise_on_error=False)
    if parsed is None:
        return cleaned, "unknown", "", ""

    sp = getattr(parsed, "species", None)
    species_latin = str(getattr(sp, "latin_name", "") or "")
    prefix = str(getattr(sp, "mhc_prefix", "") or "")
    norm = parsed.to_string()
    cls_name = type(parsed).__name__

    if cls_name in {"Allele", "Pair"}:
        resolution = "exact_allele"
    elif cls_name == "Serotype":
        resolution = "serotype"
    elif cls_name in {"MhcClass", "Class2Locus", "Species"}:
        resolution = "class_only"
    elif cls_name == "Haplotype":
        resolution = "haplotype"
    elif cls_name == "Gene":
        resolution = "gene"
    else:
        resolution = "unknown"
    return norm, resolution, species_latin, prefix


def _normalize_requested_alleles(alleles: Iterable[str] | None) -> dict[str, str]:
    """Return mapping ``normalized_allele -> original_requested_allele``."""
    if alleles is None:
        return {}
    mapping: dict[str, str] = {}
    for raw in alleles:
        token = _clean_str(raw)
        if not token:
            continue
        norm, resolution, _, _ = parse_mhc_token(token)
        key = norm if resolution == "exact_allele" and norm else token
        mapping.setdefault(key, token)
    return mapping


def _parse_token_field(raw: str) -> list[tuple[str, str, str]]:
    """Split semicolon/comma-separated MHC tokens and parse each via MHCgnomes."""
    if not raw:
        return []
    out: list[tuple[str, str, str]] = []
    for chunk in raw.replace(",", ";").split(";"):
        piece = chunk.strip()
        if not piece:
            continue
        norm, res, sp, _ = parse_mhc_token(piece)
        out.append((norm or piece, res, sp))
    return out


def _inspect_observation_mhc(row: dict[str, object]) -> dict[str, object]:
    """Classify one positive-MS observation's MHC restriction and donor-set context."""
    raw_restriction = _clean_str(row.get("mhc_restriction", ""))
    raw_allele_set = _clean_str(row.get("mhc_allele_set", ""))
    provenance = _clean_str(row.get("mhc_allele_provenance", "")).lower()
    is_mono = _is_truthy(row.get("is_monoallelic", False))

    restriction_tokens = _parse_token_field(raw_restriction)
    allele_set_tokens = _parse_token_field(raw_allele_set)

    exact_from_restriction = [t[0] for t in restriction_tokens if t[1] == "exact_allele"]
    exact_from_set = [t[0] for t in allele_set_tokens if t[1] == "exact_allele"]
    non_exact_tokens = [
        (t[0], t[1])
        for t in restriction_tokens
        if t[1] in {"serotype", "class_only", "haplotype", "gene"}
    ]
    observed_species = sorted({t[2] for t in (*restriction_tokens, *allele_set_tokens) if t[2]})

    is_multiallelic_context = (
        provenance in _MULTIALLELIC_PROVENANCES
        or len(set(exact_from_restriction)) > 1
        or (len(set(exact_from_set)) > 1 and not is_mono)
    )

    exact_measured: set[str] = set()
    donor_set: set[str] = set()
    if (is_mono and len(set(exact_from_restriction)) == 1) or (
        not is_multiallelic_context and len(set(exact_from_restriction)) == 1
    ):
        exact_measured = set(exact_from_restriction)
    else:
        donor_set = set(exact_from_restriction) | set(exact_from_set)

    if exact_measured:
        resolution = "exact_allele"
        normalized = ";".join(sorted(exact_measured))
    elif donor_set:
        resolution = "donor_set"
        normalized = ";".join(sorted(donor_set))
    elif non_exact_tokens:
        resolutions = {r for _, r in non_exact_tokens}
        if "serotype" in resolutions:
            resolution = "serotype"
        elif "class_only" in resolutions:
            resolution = "class_only"
        elif "haplotype" in resolutions:
            resolution = "haplotype"
        else:
            resolution = "gene"
        normalized = ";".join(sorted({n for n, _ in non_exact_tokens}))
    else:
        resolution = "unknown"
        normalized = ""

    return {
        "raw_restriction": raw_restriction,
        "normalized": normalized,
        "resolution": resolution,
        "provenance": provenance,
        "exact_measured": exact_measured,
        "donor_set": donor_set,
        "observed_mhc_species": ";".join(observed_species),
    }


def _extract_predicted_map(
    predicted_assignments: pd.DataFrame | None,
    requested_norm_map: dict[str, str],
) -> dict[str, set[str]]:
    """Build ``peptide -> set(normalized_requested_alleles)`` from explicit predictions."""
    if predicted_assignments is None or predicted_assignments.empty:
        return {}
    if "peptide" not in predicted_assignments.columns:
        raise ValueError("predicted_assignments must include a 'peptide' column")
    allele_col = (
        "allele"
        if "allele" in predicted_assignments.columns
        else "best_allele"
        if "best_allele" in predicted_assignments.columns
        else None
    )
    if allele_col is None:
        raise ValueError("predicted_assignments must include an 'allele' or 'best_allele' column")

    out: dict[str, set[str]] = {}
    for row in predicted_assignments[["peptide", allele_col]].itertuples(index=False):
        pep = _clean_str(row[0])
        raw_allele = _clean_str(row[1])
        if not pep or not raw_allele:
            continue
        norm, res, _, _ = parse_mhc_token(raw_allele)
        key = norm if res == "exact_allele" and norm else raw_allele
        if requested_norm_map and key not in requested_norm_map:
            continue
        out.setdefault(pep, set()).add(key)
    return out


def attribute_ms_evidence(
    candidates: pd.DataFrame | Iterable[str],
    observations: pd.DataFrame | None = None,
    *,
    requested_alleles: Iterable[str] | None = None,
    predicted_assignments: pd.DataFrame | None = None,
    species_context: SpeciesContext | None = None,
    host_species: str = "human",
    mhc_species: str | None = None,
    antigen_source_species: str | None = None,
    evidence_scope: str = "all",
    evidence_loaded: bool = True,
    dataset_version: str | None = None,
) -> pd.DataFrame:
    """Attribute positive MS observations to candidate targets with explicit MHC tiers (#216).

    Matches **exact peptides only** (longer ligands never transfer to nested
    unobserved peptides), enforces Hitlist positive-MS admission, preserves all
    candidate target/protein/variant columns, and keeps exact measured MHC,
    ambiguous multi-allelic donor-set context, computational predicted
    restriction, class/serotype-only restriction, and any-MHC peptide evidence
    strictly distinct.

    Parameters
    ----------
    candidates
        Candidate DataFrame with a ``peptide`` column, or an iterable of
        peptide strings. All input candidate columns are preserved.
    observations
        Row-level Hitlist/public MS observation DataFrame (or ``None`` when
        evidence was not loaded).
    requested_alleles
        Optional requested molecular MHC alleles (e.g. patient HLA / mouse H-2
        / dog DLA alleles). Required when ``evidence_scope="requested_mhc"``.
    predicted_assignments
        Optional explicit computational prediction DataFrame (``peptide`` and
        ``allele``/``best_allele``) used solely to label
        ``predicted_restriction`` when measured exact/donor-set restriction is
        absent. Never conflated with measured MHC restriction.
    species_context
        Optional :class:`tsarina.species.SpeciesContext`. If omitted, resolved
        from *host_species*, *mhc_species*, and *antigen_source_species*.
    evidence_scope
        One of ``"all"`` (default, retain all candidates and observation tiers),
        ``"requested_mhc"`` (retain only observations supporting a requested
        allele via exact, donor-set, or explicit predicted restriction), or
        ``"any_mhc"`` (retain all positive exact-peptide MS observations).
    evidence_loaded
        Set ``False`` when MS evidence was intentionally skipped or unavailable
        so rows are marked ``evidence_unavailable`` (``not_loaded``) rather
        than ``no_ms_match`` (``no_match``).
    dataset_version
        Optional explicit MS dataset version string (falls back to
        ``observations.attrs`` when present).
    """
    if evidence_scope not in VALID_EVIDENCE_SCOPES:
        raise ValueError(
            f"evidence_scope must be one of {sorted(VALID_EVIDENCE_SCOPES)}, got {evidence_scope!r}"
        )

    requested_norm_map = _normalize_requested_alleles(requested_alleles)
    if evidence_scope == "requested_mhc" and not requested_norm_map:
        raise ValueError("evidence_scope='requested_mhc' requires non-empty requested_alleles")

    if isinstance(candidates, pd.DataFrame):
        if "peptide" not in candidates.columns:
            raise ValueError("candidates DataFrame must contain a 'peptide' column")
        cand_df = candidates.copy().reset_index(drop=True)
    else:
        peps = [_clean_str(p) for p in candidates if _clean_str(p)]
        cand_df = pd.DataFrame({"peptide": peps})

    ctx = species_context or resolve_species_context(
        host_species,
        antigen_source_species=antigen_source_species,
        mhc_species=mhc_species,
    )
    cand_cols = list(cand_df.columns)
    all_cols = cand_cols + [c for c in _ATTRIBUTION_COLUMNS if c not in cand_cols]

    resolved_dataset_version = dataset_version
    if resolved_dataset_version is None and observations is not None:
        resolved_dataset_version = str(
            observations.attrs.get("dataset_version")
            or observations.attrs.get("evidence_store_package_id")
            or observations.attrs.get("hitlist_version")
            or ""
        )

    def _finalize(frame: pd.DataFrame) -> pd.DataFrame:
        out = frame.reindex(columns=all_cols).reset_index(drop=True)
        out.attrs.update(
            {
                "species_context": ctx.to_provenance(),
                "evidence_scope": evidence_scope,
                "evidence_loaded": bool(evidence_loaded and observations is not None),
                "requested_alleles": list(requested_norm_map.values()),
                "dataset_version": resolved_dataset_version or "",
            }
        )
        return out

    if cand_df.empty:
        return _finalize(pd.DataFrame(columns=all_cols))

    def _empty_attr_row(status: str, tier: str) -> dict[str, object]:
        return {
            "ms_evidence_status": status,
            "ms_attribution_tier": tier,
            "matched_requested_alleles": "",
            "mhc_restriction_raw": "",
            "mhc_restriction_normalized": "",
            "mhc_resolution": "",
            "mhc_allele_provenance": "",
            "donor_set_alleles": "",
            "observed_mhc_species": "",
            "host_species": ctx.host_species,
            "antigen_source_species": ctx.antigen_source_species,
            "mhc_species": ctx.mhc_species,
            "assay_modality": "",
            "ms_outcome_status": "",
            "is_ms_observation": False,
            "source_tissue": "",
            "disease": "",
            "src_cancer": False,
            "src_healthy_tissue": False,
            "src_healthy_reproductive": False,
            "src_healthy_reproductive_female": False,
            "src_healthy_reproductive_male": False,
            "src_cell_line": False,
            "src_ex_vivo": False,
            "pmid": "",
            "assay_iri": "",
            "provenance_id": "",
            "source": "",
            "dataset_version": resolved_dataset_version or "",
        }

    if not evidence_loaded or observations is None:
        if evidence_scope != "all":
            return _finalize(pd.DataFrame(columns=all_cols))
        rows = [
            {**rec, **_empty_attr_row("not_loaded", EVIDENCE_UNAVAILABLE)}
            for rec in cand_df.to_dict("records")
        ]
        return _finalize(pd.DataFrame(rows))

    # Enforce Hitlist positive-MS admission and exact peptide match.
    obs = observations.copy()
    if not obs.empty:
        obs = filter_positive_ms_observations(obs, include_binding=False, include_other=False)
        if not obs.empty and "peptide" in obs.columns:
            cand_peptides = set(cand_df["peptide"].astype(str))
            obs = obs[obs["peptide"].astype(str).isin(cand_peptides)].copy()

    # Attach normalized Hitlist assay metadata if not already present.
    if not obs.empty:
        from hitlist.assays import ASSAY_COLUMNS, annotate_assays

        if (
            "assay_modality" not in obs.columns
            or "ms_outcome_status" not in obs.columns
            or "is_ms_observation" not in obs.columns
        ):
            has_assay_cols = bool(set(ASSAY_COLUMNS) & set(obs.columns))
            if has_assay_cols:
                ann = annotate_assays(obs.reindex(columns=ASSAY_COLUMNS))
                for col in ("assay_modality", "ms_outcome_status", "is_ms_observation"):
                    if col not in obs.columns:
                        obs[col] = ann[col].values
            else:
                if "assay_modality" not in obs.columns:
                    obs["assay_modality"] = "ms"
                if "ms_outcome_status" not in obs.columns:
                    obs["ms_outcome_status"] = "positive_result"
                if "is_ms_observation" not in obs.columns:
                    obs["is_ms_observation"] = True

    predicted_map = _extract_predicted_map(predicted_assignments, requested_norm_map)

    obs_by_peptide: dict[str, list[dict[str, object]]] = {}
    if not obs.empty and "peptide" in obs.columns:
        for rec in obs.to_dict("records"):
            pep = _clean_str(rec.get("peptide", ""))
            if pep:
                obs_by_peptide.setdefault(pep, []).append(rec)

    requested_keys = set(requested_norm_map.keys())
    out_rows: list[dict[str, object]] = []

    for cand in cand_df.to_dict("records"):
        pep = _clean_str(cand.get("peptide", ""))
        matched_obs = obs_by_peptide.get(pep, [])
        if not matched_obs:
            if evidence_scope == "all":
                out_rows.append({**cand, **_empty_attr_row("no_match", NO_MS_MATCH)})
            continue

        pred_alleles = predicted_map.get(pep, set())
        emitted_for_cand = False
        for rec in matched_obs:
            mhc_info = _inspect_observation_mhc(rec)
            exact_hits = mhc_info["exact_measured"] & requested_keys if requested_keys else set()
            donor_hits = mhc_info["donor_set"] & requested_keys if requested_keys else set()

            if exact_hits:
                tier = EXACT_REQUESTED_MHC
                matched_req = sorted(exact_hits)
            elif donor_hits:
                tier = DONOR_SET_AMBIGUOUS
                matched_req = sorted(donor_hits)
            elif requested_keys and pred_alleles:
                tier = PREDICTED_RESTRICTION
                matched_req = sorted(pred_alleles)
            elif mhc_info["resolution"] in {"serotype", "class_only", "haplotype"}:
                tier = CLASS_OR_SEROTYPE_ONLY
                matched_req = []
            else:
                tier = ANY_MHC_PEPTIDE
                matched_req = []

            if evidence_scope == "requested_mhc" and tier not in {
                EXACT_REQUESTED_MHC,
                DONOR_SET_AMBIGUOUS,
                PREDICTED_RESTRICTION,
            }:
                continue

            emitted_for_cand = True
            row_dataset_version = (
                _clean_str(rec.get("dataset_version", "")) or resolved_dataset_version or ""
            )
            out_rows.append(
                {
                    **cand,
                    "ms_evidence_status": "matched",
                    "ms_attribution_tier": tier,
                    "matched_requested_alleles": ";".join(matched_req),
                    "mhc_restriction_raw": mhc_info["raw_restriction"],
                    "mhc_restriction_normalized": mhc_info["normalized"],
                    "mhc_resolution": mhc_info["resolution"],
                    "mhc_allele_provenance": mhc_info["provenance"],
                    "donor_set_alleles": ";".join(sorted(mhc_info["donor_set"])),
                    "observed_mhc_species": mhc_info["observed_mhc_species"],
                    "host_species": _clean_str(rec.get("host_species", "")) or ctx.host_species,
                    "antigen_source_species": (
                        _clean_str(rec.get("antigen_source_species", ""))
                        or ctx.antigen_source_species
                    ),
                    "mhc_species": _clean_str(rec.get("mhc_species", "")) or ctx.mhc_species,
                    "assay_modality": _clean_str(rec.get("assay_modality", "")),
                    "ms_outcome_status": _clean_str(rec.get("ms_outcome_status", "")),
                    "is_ms_observation": bool(_is_truthy(rec.get("is_ms_observation", True))),
                    "source_tissue": _clean_str(rec.get("source_tissue", "")),
                    "disease": _clean_str(rec.get("disease", "")),
                    "src_cancer": _is_truthy(rec.get("src_cancer", False)),
                    "src_healthy_tissue": _is_truthy(rec.get("src_healthy_tissue", False)),
                    "src_healthy_reproductive": _is_truthy(
                        rec.get("src_healthy_reproductive", False)
                    ),
                    "src_healthy_reproductive_female": _is_truthy(
                        rec.get("src_healthy_reproductive_female", False)
                    ),
                    "src_healthy_reproductive_male": _is_truthy(
                        rec.get("src_healthy_reproductive_male", False)
                    ),
                    "src_cell_line": _is_truthy(rec.get("src_cell_line", False)),
                    "src_ex_vivo": _is_truthy(rec.get("src_ex_vivo", False)),
                    "pmid": _clean_str(rec.get("pmid", "")),
                    "assay_iri": _clean_str(rec.get("assay_iri", "")),
                    "provenance_id": _clean_str(rec.get("provenance_id", "")),
                    "source": _clean_str(rec.get("source", "")),
                    "dataset_version": row_dataset_version,
                }
            )
        if not emitted_for_cand and evidence_scope == "all":
            out_rows.append({**cand, **_empty_attr_row("no_match", NO_MS_MATCH)})

    return _finalize(pd.DataFrame(out_rows))


def summarize_ms_attribution(
    attributed: pd.DataFrame,
    *,
    group_cols: tuple[str, ...] = ("peptide",),
) -> pd.DataFrame:
    """Summarize observation-level attributed MS evidence to one row per target key (#216).

    Preserves separate counts and normalized allele lists for exact measured
    requested-MHC matches, ambiguous multi-allelic donor-set matches, explicit
    predicted restrictions, class/serotype-only observations, and any-MHC
    exact-peptide observations.
    """
    summary_cols = [
        *group_cols,
        "ms_evidence_status",
        "ms_attribution_tier",
        "ms_hit_count",
        "ms_exact_mhc_hits",
        "ms_exact_mhc_alleles",
        "ms_donor_set_hits",
        "ms_donor_set_alleles",
        "ms_predicted_restriction_alleles",
        "ms_class_or_serotype_labels",
        "ms_any_mhc_hits",
        "ms_alleles",
        "ms_allele_count",
        "ms_pmids",
        "ms_in_cancer",
        "ms_in_healthy_tissue",
        "ms_in_healthy_reproductive",
        "ms_in_healthy_reproductive_female",
        "ms_in_healthy_reproductive_male",
        "ms_in_cell_line",
        "ms_in_ex_vivo",
    ]
    if attributed.empty:
        out = pd.DataFrame(columns=summary_cols)
        out.attrs.update(attributed.attrs)
        return out

    rows: list[dict[str, object]] = []
    grouped = attributed.groupby(list(group_cols), sort=False, dropna=False)
    for keys, grp in grouped:
        if not isinstance(keys, tuple):
            keys = (keys,)
        key_dict = dict(zip(group_cols, keys))

        matched = grp[grp["ms_evidence_status"].eq("matched")]
        if matched.empty:
            statuses = set(grp["ms_evidence_status"].astype(str))
            status = "not_loaded" if "not_loaded" in statuses else "no_match"
            tier = EVIDENCE_UNAVAILABLE if status == "not_loaded" else NO_MS_MATCH
            rows.append(
                {
                    **key_dict,
                    "ms_evidence_status": status,
                    "ms_attribution_tier": tier,
                    "ms_hit_count": 0,
                    "ms_exact_mhc_hits": 0,
                    "ms_exact_mhc_alleles": "",
                    "ms_donor_set_hits": 0,
                    "ms_donor_set_alleles": "",
                    "ms_predicted_restriction_alleles": "",
                    "ms_class_or_serotype_labels": "",
                    "ms_any_mhc_hits": 0,
                    "ms_alleles": "",
                    "ms_allele_count": 0,
                    "ms_pmids": "",
                    "ms_in_cancer": False,
                    "ms_in_healthy_tissue": False,
                    "ms_in_healthy_reproductive": False,
                    "ms_in_healthy_reproductive_female": False,
                    "ms_in_healthy_reproductive_male": False,
                    "ms_in_cell_line": False,
                    "ms_in_ex_vivo": False,
                }
            )
            continue

        tiers = list(matched["ms_attribution_tier"].astype(str))
        best_tier = min(tiers, key=lambda t: _TIER_RANK.get(t, 999))

        exact_rows = matched[matched["ms_attribution_tier"].eq(EXACT_REQUESTED_MHC)]
        donor_rows = matched[matched["ms_attribution_tier"].eq(DONOR_SET_AMBIGUOUS)]
        pred_rows = matched[matched["ms_attribution_tier"].eq(PREDICTED_RESTRICTION)]
        coarse_rows = matched[matched["ms_attribution_tier"].eq(CLASS_OR_SEROTYPE_ONLY)]

        def _split_set(series: pd.Series) -> set[str]:
            vals: set[str] = set()
            for item in series:
                text = _clean_str(item)
                if text:
                    vals.update(p.strip() for p in text.split(";") if p.strip())
            return vals

        exact_alleles = _split_set(exact_rows["matched_requested_alleles"])
        donor_alleles = _split_set(donor_rows["matched_requested_alleles"])
        pred_alleles = _split_set(pred_rows["matched_requested_alleles"])
        coarse_labels = _split_set(coarse_rows["mhc_restriction_normalized"])

        # All exact/donor normalized alleles reported across matched rows.
        all_alleles: set[str] = set()
        for row in matched.itertuples(index=False):
            if row.mhc_resolution in {"exact_allele", "donor_set"}:
                all_alleles.update(
                    p.strip()
                    for p in _clean_str(row.mhc_restriction_normalized).split(";")
                    if p.strip()
                )
            all_alleles.update(
                p.strip() for p in _clean_str(row.donor_set_alleles).split(";") if p.strip()
            )
        pmids = _split_set(matched["pmid"])

        rows.append(
            {
                **key_dict,
                "ms_evidence_status": "matched",
                "ms_attribution_tier": best_tier,
                "ms_hit_count": len(matched),
                "ms_exact_mhc_hits": len(exact_rows),
                "ms_exact_mhc_alleles": ";".join(sorted(exact_alleles)),
                "ms_donor_set_hits": len(donor_rows),
                "ms_donor_set_alleles": ";".join(sorted(donor_alleles)),
                "ms_predicted_restriction_alleles": ";".join(sorted(pred_alleles)),
                "ms_class_or_serotype_labels": ";".join(sorted(coarse_labels)),
                "ms_any_mhc_hits": len(matched),
                "ms_alleles": ";".join(sorted(all_alleles)),
                "ms_allele_count": len(all_alleles),
                "ms_pmids": ";".join(sorted(pmids)),
                "ms_in_cancer": bool(matched["src_cancer"].astype(bool).any()),
                "ms_in_healthy_tissue": bool(matched["src_healthy_tissue"].astype(bool).any()),
                "ms_in_healthy_reproductive": bool(
                    matched["src_healthy_reproductive"].astype(bool).any()
                ),
                "ms_in_healthy_reproductive_female": bool(
                    matched["src_healthy_reproductive_female"].astype(bool).any()
                )
                if "src_healthy_reproductive_female" in matched.columns
                else False,
                "ms_in_healthy_reproductive_male": bool(
                    matched["src_healthy_reproductive_male"].astype(bool).any()
                )
                if "src_healthy_reproductive_male" in matched.columns
                else False,
                "ms_in_cell_line": bool(matched["src_cell_line"].astype(bool).any()),
                "ms_in_ex_vivo": bool(matched["src_ex_vivo"].astype(bool).any()),
            }
        )

    out = pd.DataFrame(rows, columns=summary_cols)
    out.attrs.update(attributed.attrs)
    return out
