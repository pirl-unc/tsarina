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


def load_public_ms_hits(
    peptides: set[str] | list[str] | tuple[str, ...],
    *,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    mhc_class: str | None = "I",
    mhc_species: str | None = "Homo sapiens",
    classify_source: bool = True,
    min_allele_resolution: str | None = None,
    drop_binding_assays: bool = True,
) -> pd.DataFrame:
    """Load public MS evidence for peptides from the canonical tsarina path.

    The default path uses hitlist's cached observations parquet.  Explicit
    ``iedb_path`` or ``cedar_path`` switches to the raw scanner for ad-hoc
    snapshots while preserving the same species/class and binding-assay policy.
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
        hits = scan_public_ms(
            peptides=peptide_set,
            iedb_path=resolved_iedb,
            cedar_path=resolved_cedar,
            mhc_species=mhc_species,
            mhc_class=mhc_class,
            classify_source=classify_source,
            min_allele_resolution=min_allele_resolution,
        )
    else:
        from .indexing import load_ms_evidence

        hits = load_ms_evidence(
            peptides=peptide_set,
            mhc_class=mhc_class,
            mhc_species=mhc_species,
            drop_binding_assays=drop_binding_assays,
        )

    if drop_binding_assays and "is_binding_assay" in hits.columns:
        hits = hits[~hits["is_binding_assay"]].copy()
    return hits.reset_index(drop=True)


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
        "ms_alleles": ("mhc_restriction", lambda x: ";".join(sorted(_nonempty_strings(x)))),
        "ms_allele_count": ("mhc_restriction", lambda x: len(_nonempty_strings(x))),
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
) -> pd.DataFrame:
    """Per-peptide healthy *somatic*-tissue MS hits for a CTA gene.

    The gene-level ``ms_restriction`` tier collapses detail that matters for
    safety; this surfaces the underlying immunopeptidome hits at
    peptide x tissue x allele granularity so consumers can do allele-aware,
    peptide-level screening (tsarina#76).

    Reproductive and thymic hits (expected for CTAs) are excluded — only
    ``src_healthy_tissue`` (non-reproductive, non-thymus) rows are returned.

    Parameters
    ----------
    gene
        Gene symbol (or alias — resolved via the hitlist peptide mappings).
    mhc_class, mhc_species
        Passed through to :func:`tsarina.indexing.load_ms_evidence`.

    Returns
    -------
    pd.DataFrame
        One row per healthy-somatic MS observation, columns:
        ``peptide``, ``tissue``, ``allele`` (recorded ``mhc_restriction``),
        ``allele_set`` (candidate set for multiallelic samples),
        ``provenance``, ``pmid``, and ``vital_organ`` (bool — tissue matches a
        ``SAFETY_TISSUE_GROUPS`` vital organ).  Empty if the gene has no
        healthy-somatic MS hits.
    """
    from .indexing import load_ms_evidence

    df = load_ms_evidence(gene_name=gene, mhc_class=mhc_class, mhc_species=mhc_species)
    if df.empty or "src_healthy_tissue" not in df.columns:
        return pd.DataFrame(columns=CTA_HEALTHY_TISSUE_MS_COLUMNS)

    hits = df[df["src_healthy_tissue"].astype(bool)].copy()
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
