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

"""Shared public-MS evidence loading and peptide-level aggregation."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .tiers import SAFETY_TISSUE_GROUPS, VITAL_TISSUE_MS_NAMES

SOURCE_FLAG_OUTPUTS: dict[str, str] = {
    "src_cancer": "ms_cancer",
    "src_healthy_tissue": "ms_healthy_tissue",
    "src_healthy_thymus": "ms_healthy_thymus",
    "src_healthy_reproductive": "ms_healthy_reproductive",
    "src_healthy_reproductive_female": "ms_healthy_reproductive_female",
    "src_healthy_reproductive_male": "ms_healthy_reproductive_male",
    "src_cell_line": "ms_cell_line",
    "src_ebv_lcl": "ms_ebv_lcl",
    "src_ex_vivo": "ms_ex_vivo",
}


def _peptide_set(peptides: set[str] | list[str] | tuple[str, ...]) -> set[str]:
    return set(peptides)


def _nonempty_strings(values) -> set[str]:
    strings = set()
    for value in values:
        if not isinstance(value, str):
            continue
        stripped = value.strip()
        if stripped:
            strings.add(stripped)
    return strings


#: Class-only / unresolved MHC labels that carry no allele information at
#: all, dropped entirely rather than counted as if they were a real (if
#: coarse) allele the way a serotype like "HLA-B7" is.
_UNINFORMATIVE_ALLELE_LABELS = frozenset(
    {"hla class i", "hla class ii", "class i", "class ii", "unknown"}
)


def _individual_alleles(values) -> set[str]:
    """Split each row's MHC restriction string into individual allele
    tokens and deduplicate across rows.

    A single row's ``mhc_restriction`` is sometimes itself a semicolon-
    joined multi-donor context (a pooled-sample study reporting its whole
    cohort's HLA typing) rather than one allele. Joining those raw strings
    across many rows without splitting them first (the previous behavior)
    pastes the same alleles together many times over -- a peptide observed
    across ten mostly-overlapping ten-allele donor panels reported "ten
    copies of nearly the same ten alleles" instead of the handful of truly
    distinct ones.
    """
    alleles: set[str] = set()
    for value in values:
        if not isinstance(value, str):
            continue
        for token in value.split(";"):
            token = token.strip()
            if token and token.lower() not in _UNINFORMATIVE_ALLELE_LABELS:
                alleles.add(token)
    return alleles


def filter_positive_ms_observations(
    hits: pd.DataFrame,
    *,
    include_binding: bool = False,
    include_other: bool = False,
) -> pd.DataFrame:
    """Filter raw or cached observations using Hitlist's positive-MS contract (#199).

    Uses ``is_ms_observation`` when present, otherwise recomputes admission
    through :func:`hitlist.assays.annotate_assays` when assay metadata columns
    are available. Structural, unknown, conflicting, and negative-MS records
    are excluded from MS-only queries rather than admitted by the complement
    of ``is_binding_assay``. Binding inclusion stays separate from other
    non-MS modalities.
    """
    if hits.empty or (include_binding and include_other):
        return hits

    from hitlist.assays import ASSAY_COLUMNS, annotate_assays

    has_assay_cols = bool(set(ASSAY_COLUMNS) & set(hits.columns))
    annotated = annotate_assays(hits.reindex(columns=ASSAY_COLUMNS)) if has_assay_cols else None

    if "is_ms_observation" in hits.columns:
        ms_mask = hits["is_ms_observation"].fillna(False).astype(bool)
    elif annotated is not None:
        ms_mask = annotated["is_ms_observation"].astype(bool)
    elif "is_binding_assay" in hits.columns:
        ms_mask = ~hits["is_binding_assay"].fillna(False).astype(bool)
    else:
        ms_mask = pd.Series(True, index=hits.index)

    if annotated is not None:
        binding_mask = annotated["is_binding_assay"].astype(bool)
    elif "is_binding_assay" in hits.columns:
        binding_mask = hits["is_binding_assay"].fillna(False).astype(bool)
    else:
        binding_mask = pd.Series(False, index=hits.index)

    if include_binding:
        keep = ms_mask | binding_mask
    elif include_other:
        keep = ms_mask | ~binding_mask
    else:
        keep = ms_mask
    return hits[keep].copy()


def load_public_ms_hits(
    peptides: set[str] | list[str] | tuple[str, ...],
    *,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    evidence_package: str | Path | None = None,
    mhc_class: str | None = "I",
    mhc_species: str | None = "Homo sapiens",
    classify_source: bool = True,
    min_allele_resolution: str | None = None,
    drop_binding_assays: bool = True,
    include_other_assays: bool = False,
    auto_build: bool = True,
) -> pd.DataFrame:
    """Load public MS evidence for peptides from the canonical tsarina path.

    The default path uses a configured Hitlist 2.x ``EvidenceStore`` package
    (via *evidence_package* or ``TSARINA_EVIDENCE_PACKAGE`` /
    ``HITLIST_EVIDENCE_PACKAGE``) or hitlist's cached observations parquet.
    Explicit ``iedb_path`` or ``cedar_path`` switches to the raw scanner for
    ad-hoc snapshots while enforcing Hitlist's positive-MS admission contract
    (excluding binding, structural, unknown, and negative-MS assays by default).
    """
    peptide_set = _peptide_set(peptides)
    if iedb_path is not None or cedar_path is not None:
        from .datasources import resolve_dataset_paths
        from .iedb import scan_public_ms

        resolved_iedb, resolved_cedar = resolve_dataset_paths(
            iedb_path,
            cedar_path,
            require_iedb=iedb_path is not None,
        )
        scan_kwargs = {
            "peptides": peptide_set,
            "iedb_path": resolved_iedb,
            "cedar_path": resolved_cedar,
            "mhc_species": mhc_species,
            "mhc_class": mhc_class,
            "classify_source": classify_source,
            "min_allele_resolution": min_allele_resolution,
        }
        if not drop_binding_assays or include_other_assays:
            scan_kwargs["positive_ms_only"] = False
        hits = scan_public_ms(**scan_kwargs)
    else:
        from .indexing import load_ms_evidence

        load_kwargs = {
            "peptides": peptide_set,
            "mhc_class": mhc_class,
            "mhc_species": mhc_species,
            "drop_binding_assays": drop_binding_assays,
            "auto_build": auto_build,
        }
        if evidence_package is not None:
            load_kwargs["evidence_package"] = evidence_package
        hits = load_ms_evidence(**load_kwargs)

    attrs = dict(hits.attrs)
    hits = filter_positive_ms_observations(
        hits,
        include_binding=not drop_binding_assays,
        include_other=include_other_assays,
    )
    out = hits.reset_index(drop=True)
    out.attrs.update(attrs)
    return out


def aggregate_ms_hits_by_peptide(
    hits: pd.DataFrame,
    *,
    source_flag_outputs: dict[str, str] | None = None,
    cell_line_output: str | None = "ms_cell_lines",
) -> pd.DataFrame:
    """Aggregate hitlist evidence rows to one row per peptide."""
    columns = ["peptide", "ms_hit_count", "ms_alleles", "ms_allele_count"]
    source_flag_outputs = (
        SOURCE_FLAG_OUTPUTS if source_flag_outputs is None else source_flag_outputs
    )
    columns.extend(source_flag_outputs.values())
    if cell_line_output is not None:
        columns.append(cell_line_output)
    if hits.empty or "peptide" not in hits.columns:
        return pd.DataFrame(columns=columns)

    agg: dict[str, tuple] = {
        "ms_hit_count": ("peptide", "size"),
        "ms_alleles": ("mhc_restriction", lambda x: ";".join(sorted(_individual_alleles(x)))),
        "ms_allele_count": ("mhc_restriction", lambda x: len(_individual_alleles(x))),
    }
    for source_col, output_col in source_flag_outputs.items():
        if source_col in hits.columns:
            agg[output_col] = (source_col, "any")
    if cell_line_output is not None and "cell_line_name" in hits.columns:
        agg[cell_line_output] = (
            "cell_line_name",
            lambda x: ";".join(sorted(_nonempty_strings(x))),
        )

    return hits.groupby("peptide", as_index=False).agg(**agg)


def aggregate_ms_hits_for_iedb_columns(hits: pd.DataFrame) -> pd.DataFrame:
    """Aggregate hits to the legacy ``iedb_*`` peptide columns."""
    aggregated = aggregate_ms_hits_by_peptide(
        hits,
        source_flag_outputs={},
        cell_line_output=None,
    )
    if aggregated.empty:
        return pd.DataFrame(columns=["peptide", "iedb_hit_count", "iedb_alleles"])
    return aggregated.rename(
        columns={
            "ms_hit_count": "iedb_hit_count",
            "ms_alleles": "iedb_alleles",
        }
    )[["peptide", "iedb_hit_count", "iedb_alleles"]]


#: Exact vital-organ tissue names, taken verbatim from the RNA-side safety
#: vocabulary so this MS screen can never drift from it: every brain subregion
#: in :data:`tiers.SAFETY_TISSUE_GROUPS` (plus heart/lung/liver/pancreas) and the
#: MS-specific spellings in :data:`tiers.VITAL_TISSUE_MS_NAMES` (e.g.
#: ``"central nervous system (cns)"``).  Any one of these is a vital organ.
_VITAL_ORGAN_EXACT_NAMES: frozenset[str] = frozenset().union(*SAFETY_TISSUE_GROUPS.values()) | {
    name.lower() for name in VITAL_TISSUE_MS_NAMES
}

#: Stems flagging an MS source tissue / cell name as a vital organ
#: (brain/CNS, heart, lung, liver, pancreas) -- the MS-side counterpart of
#: :data:`tiers.SAFETY_TISSUE_GROUPS`.  Stems (not whole words) so spelling
#: variants and inflections are caught: ``hepat`` -> hepatic/hepatocyte,
#: ``cerebr`` -> cerebral/cerebrum, ``thalam`` -> thalamus/hypothalamus,
#: ``myocard`` -> myocardium/myocardial, ``pulmonar`` -> pulmonary.  Matched as
#: plain substrings so compound names match too (``brain`` -> midbrain/
#: brainstem; ``medulla oblongata`` avoids bare ``medulla``, which also names
#: the non-vital adrenal/renal medulla).  Brain/CNS gets the widest net because
#: a missed vital-organ hit (false negative) is worse for a safety screen than
#: an over-flag.
_VITAL_ORGAN_SUBSTRINGS: tuple[str, ...] = (
    # brain / CNS
    "brain",
    "cerebr",
    "cerebell",
    "central nervous system",
    "spinal cord",
    "hippocamp",
    "amygdala",
    "thalam",
    "basal ganglia",
    "choroid plexus",
    "medulla oblongata",
    "white matter",
    "retina",
    # heart
    "heart",
    "cardiac",
    "myocard",
    # lung
    "lung",
    "pulmonar",
    # liver
    "hepat",
    # pancreas
    "pancrea",
    "islet",
)

#: Short/ambiguous vital-organ stems matched at a word boundary so they cannot
#: appear spuriously inside an unrelated tissue: ``\bliver`` does not match
#: ``deliver``; ``cns``/``pons`` stay whole words.
_VITAL_ORGAN_WORD_BOUNDED: tuple[str, ...] = ("cns", "liver", "pons")

_VITAL_ORGAN_REGEX = re.compile(
    "|".join(
        [
            *(re.escape(stem) for stem in _VITAL_ORGAN_SUBSTRINGS),
            *(rf"\b{re.escape(stem)}" for stem in _VITAL_ORGAN_WORD_BOUNDED),
        ]
    )
)

#: Organ qualifiers that make a ``"... cortex"`` non-CNS, and therefore NOT a
#: vital organ.  ``"renal"`` also covers ``"adrenal"`` (substring), but both are
#: listed for clarity; ``"ovar"`` covers ovary/ovarian, ``"thym"`` thymic/thymus,
#: ``"lymph"``/``"nodal"`` lymph-node cortex.  Everything else containing
#: ``"cortex"`` (cerebral, frontal, temporal, motor, visual, ...) is treated as
#: brain -- an inclusive bias appropriate for a safety screen.
_NON_CNS_CORTEX_QUALIFIERS: tuple[str, ...] = (
    "adrenal",
    "renal",
    "kidney",
    "ovar",
    "thym",
    "lymph",
    "nodal",
)


def _is_vital_organ_tissue(tissue: str) -> bool:
    """Whether an MS source-tissue / cell name denotes a vital organ.

    Vital organs are brain/CNS, heart, lung, liver, and pancreas -- the
    MS-side counterpart of :data:`tiers.SAFETY_TISSUE_GROUPS`.  Matching is
    case-insensitive and whitespace-normalized.
    """
    text = " ".join(str(tissue).lower().split())
    if not text:
        return False
    if text in _VITAL_ORGAN_EXACT_NAMES:
        return True
    if "cortex" in text and not any(q in text for q in _NON_CNS_CORTEX_QUALIFIERS):
        return True
    return _VITAL_ORGAN_REGEX.search(text) is not None


CTA_HEALTHY_TISSUE_MS_COLUMNS: list[str] = [
    "peptide",
    "tissue",
    "allele",
    "allele_set",
    "provenance",
    "pmid",
    "vital_organ",
]


def cta_healthy_tissue_ms_hits(
    gene: str,
    *,
    mhc_class: str | None = "I",
    mhc_species: str | None = "Homo sapiens",
    evidence_package: str | Path | None = None,
    patient_sex: str | None = None,
    sex_eligibility: str = "all_adult",
) -> pd.DataFrame:
    """Per-peptide healthy *somatic*-tissue MS hits for a CTA gene.

    The gene-level ``ms_restriction`` tier collapses detail that matters for
    safety; this surfaces the underlying immunopeptidome hits at
    peptide x tissue x allele granularity so consumers can do allele-aware,
    peptide-level screening (tsarina#76).

    Core CTA reproductive and thymic hits (``src_healthy_reproductive``:
    testis, placenta, embryo, epididymis) are excluded. When ``patient_sex`` is
    provided with ``sex_eligibility="all_adult"``, Hitlist's sex-specific
    normal reproductive tract flags (``src_healthy_reproductive_female`` for
    female patients; ``src_healthy_reproductive_male`` for male patients) are
    also included as same-sex normal-organ hits while opposite-sex reproductive
    hits remain excluded. Under ``sex_eligibility="sex_organ_only"``, only
    non-reproductive somatic hits (``src_healthy_tissue``) are returned.

    Parameters
    ----------
    gene
        Gene symbol (or alias — resolved via the hitlist peptide mappings).
    mhc_class, mhc_species, evidence_package
        Passed through to :func:`tsarina.indexing.load_ms_evidence`.
    patient_sex
        Optional biological sex (``"female"`` or ``"male"``).
    sex_eligibility
        ``"all_adult"`` (default) or ``"sex_organ_only"``.

    Returns
    -------
    pd.DataFrame
        One row per healthy-somatic (plus same-sex adult reproductive when
        requested) MS observation, columns:
        ``peptide``, ``tissue``, ``allele`` (recorded ``mhc_restriction``),
        ``allele_set`` (candidate set for multiallelic samples),
        ``provenance``, ``pmid``, and ``vital_organ`` (bool — tissue matches a
        ``SAFETY_TISSUE_GROUPS`` vital organ).  Empty if no healthy MS hits
        match.
    """
    from .indexing import load_ms_evidence

    norm_sex: str | None = None
    norm_elig = "all_adult"
    if patient_sex is not None:
        from .cta_selection import _normalize_patient_sex, _normalize_sex_eligibility

        norm_sex = _normalize_patient_sex(patient_sex)
        norm_elig = _normalize_sex_eligibility(sex_eligibility)

    load_kwargs = {"gene_name": gene, "mhc_class": mhc_class, "mhc_species": mhc_species}
    if evidence_package is not None:
        load_kwargs["evidence_package"] = evidence_package
    df = load_ms_evidence(**load_kwargs)
    if df.empty or "src_healthy_tissue" not in df.columns:
        return pd.DataFrame(columns=CTA_HEALTHY_TISSUE_MS_COLUMNS)

    mask = df["src_healthy_tissue"].fillna(False).astype(bool)
    if norm_sex is not None and norm_elig == "all_adult":
        sex_col = (
            "src_healthy_reproductive_female"
            if norm_sex == "female"
            else "src_healthy_reproductive_male"
        )
        if sex_col in df.columns:
            mask = mask | df[sex_col].fillna(False).astype(bool)

    hits = df[mask].copy()
    if hits.empty:
        return pd.DataFrame(columns=CTA_HEALTHY_TISSUE_MS_COLUMNS)

    tissue = hits.get("source_tissue", pd.Series("", index=hits.index)).fillna("").astype(str)
    if "cell_name" in hits.columns:
        cell = hits["cell_name"].fillna("").astype(str)
        tissue = tissue.where(tissue.str.strip() != "", cell)
    vital = tissue.map(_is_vital_organ_tissue)

    out = pd.DataFrame(
        {
            "peptide": hits["peptide"].astype(str).values,
            "tissue": tissue.values,
            "allele": hits.get("mhc_restriction", pd.Series("", index=hits.index)).values,
            "allele_set": hits.get("mhc_allele_set", pd.Series("", index=hits.index)).values,
            "provenance": hits.get("mhc_allele_provenance", pd.Series("", index=hits.index)).values,
            "pmid": hits.get("pmid", pd.Series("", index=hits.index)).values,
            "vital_organ": vital.values,
        }
    )
    return out.reset_index(drop=True)


C_TERM_CYS_MS_NOTE = (
    "C-terminal Cys: unmodified elution MS frequently misses oxidized/cysteinylated "
    "C-term Cys ligands (e.g. NY-ESO-1 SLLMWITQC)"
)
INTERNAL_CYS_MS_NOTE = (
    "Contains Cys: unmodified elution MS systematically under-detects cysteine "
    "oxidation/disulfide peptidoforms"
)


def annotate_ms_detectability(
    df: pd.DataFrame,
    *,
    peptide_col: str = "peptide",
    ms_hit_col: str = "ms_hit_count",
) -> pd.DataFrame:
    """Annotate peptides with cysteine/PTM unmodified-MS detectability caveats.

    Standard HLA-I immunopeptidomics searches look for unmodified peptides and
    systematically under-detect cysteine-containing ligands — especially
    peptides ending in a C-terminal cysteine (e.g. the clinically validated
    HLA-A*02:01 NY-ESO-1 epitope ``SLLMWITQC``), which readily undergo
    cysteinylation, oxidation, or disulfide formation during elution.

    Adds three columns to a copy of ``df`` without changing row order, scores,
    or tiers:

    - ``has_cysteine`` (bool): peptide contains ``"C"``.
    - ``c_terminal_cysteine`` (bool): peptide ends with ``"C"``.
    - ``ms_detectability_note`` (str): populated only when a peptide contains
      cysteine and has no observed MS hits (``ms_hit_col`` absent, NaN, or 0);
      empty string otherwise.
    """
    out = df.copy()
    if out.empty or peptide_col not in out.columns:
        out["has_cysteine"] = pd.Series(dtype=bool)
        out["c_terminal_cysteine"] = pd.Series(dtype=bool)
        out["ms_detectability_note"] = pd.Series(dtype=str)
        return out

    peps = out[peptide_col].fillna("").astype(str).str.strip().str.upper()
    has_cys = peps.str.contains("C", regex=False)
    c_term_cys = peps.str.endswith("C")

    if ms_hit_col in out.columns:
        zero_ms = pd.to_numeric(out[ms_hit_col], errors="coerce").fillna(0).le(0)
    else:
        zero_ms = pd.Series(True, index=out.index)

    notes = pd.Series("", index=out.index, dtype=object)
    notes = notes.mask(has_cys & zero_ms, INTERNAL_CYS_MS_NOTE)
    notes = notes.mask(c_term_cys & zero_ms, C_TERM_CYS_MS_NOTE)

    out["has_cysteine"] = has_cys.to_numpy(dtype=bool)
    out["c_terminal_cysteine"] = c_term_cys.to_numpy(dtype=bool)
    out["ms_detectability_note"] = notes.astype(str).to_numpy()
    return out
