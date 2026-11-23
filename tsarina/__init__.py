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

"""tsarina: shared cancer immunotherapy target selection and vaccine design.

Tsarina joins gene-level tumor specificity (``oncoref``), peptide-level public
immunopeptidomics MS evidence (``hitlist``), reference proteomes and viral ORFs
(``pyensembl`` / ``orfdata``), MHC allele parsing (``mhcgnomes``), and
presentation/processing scoring (``mhcflurry`` / ``mhctools`` / ``topiary`` /
``pepsickle`` / ``dnachisel``) into four high-level workflows.

Headline workflows
------------------
- ``personalized_targets`` — rank CTA, fusion, mutant, and viral pMHCs for a
  single patient (human or mouse, with optional ``patient_sex`` RNA + MS safety).
- ``spanning_pmhc_set`` — build a population CTA x HLA panel with MS evidence
  tiers, safety flags, and population coverage (CLI: ``tsarina panel``).
- ``design_vaccine`` / ``VaccineConfig`` — select mortality-weighted CTA
  segments and assemble an audited DNA/RNA vaccine construct.
- ``cta_healthy_tissue_ms_hits`` — inspect per-peptide somatic (and optional
  same-sex reproductive) healthy-tissue MS hits for a gene (CLI: ``tsarina hits``).

Target & allele catalogs
------------------------
- ``CTA_selection`` / ``CTA_selection_coverage`` — current OncoRef CTA decisions
  (human, mouse, dog) and readiness metadata.
- ``CTA_sex_restricted_antigens`` / ``CTA_sex_eligible_gene_names`` — biological
  sex eligibility across human and mouse reproductive/CTA sets.
- ``cancer_fusions`` / ``fusion_junction_peptides`` — OncoRef curated cancer
  fusions and host-exclusive breakpoint k-mer enumeration.
- ``HOTSPOT_MUTATIONS`` / ``mutant_peptides`` — curated hotspot driver mutations.
- ``viral_peptides`` / ``orfdata_viral_subsets`` / ``orfdata_viral_peptides`` —
  built-in UniProt oncogenic viruses and optional ``orfdata`` IARC GDC panels.
- ``panel_names`` / ``get_panel`` — built-in population HLA allele panels.

MS evidence & scoring
---------------------
- ``attribute_ms_evidence`` / ``summarize_ms_attribution`` — exact vs. donor-set
  vs. predicted MHC restriction attribution.
- ``annotate_ms_detectability`` — flag zero-MS cysteine / C-terminal Cys ligands.
- ``score_presentation`` / ``score_affinity`` — MHC presentation & affinity scoring.

Public exports resolve lazily (PEP 562) so ``import tsarina`` stays fast.
"""

import importlib
from typing import TYPE_CHECKING

from .version import __version__

if TYPE_CHECKING:
    from oncoref.cta import cta_evidence as CTA_hpa_evidence

    from .alleles import IEDB27_AB, PANEL_DEFINITIONS, get_panel, panel_names
    from .cancer_expression import (
        cta_cancer_expression_features,
        hpa_cancer_ihc_prevalence,
        hpa_cancer_rna_prevalence,
    )
    from .cta_selection import (
        CTA_TIERS,
        CTA_selection,
        CTA_selection_coverage,
        CTA_sex_eligible_gene_ids,
        CTA_sex_eligible_gene_names,
        CTA_sex_restricted_antigens,
    )
    from .evidence import CTA_detailed_evidence, CTA_evidence
    from .evidence_store import evidence_store_provenance, load_evidence_store_ms_hits
    from .fusions import cancer_fusions, fusion_junction_peptides
    from .gene_sets import (
        CTA_by_axes,
        CTA_by_tier,
        CTA_cancer_antigen_gene_ids,
        CTA_cancer_antigen_gene_names,
        CTA_clinical_target_gene_ids,
        CTA_clinical_target_gene_names,
        CTA_excluded_gene_ids,
        CTA_excluded_gene_names,
        CTA_filtered_gene_ids,
        CTA_filtered_gene_names,
        CTA_gene_ids,
        CTA_gene_names,
        CTA_never_expressed_gene_ids,
        CTA_never_expressed_gene_names,
        CTA_placental_restricted_gene_ids,
        CTA_placental_restricted_gene_names,
        CTA_relaxed_reproductive_gene_ids,
        CTA_relaxed_reproductive_gene_names,
        CTA_reproductive_antigen_gene_ids,
        CTA_reproductive_antigen_gene_names,
        CTA_target_candidate_gene_ids,
        CTA_target_candidate_gene_names,
        CTA_testis_restricted_gene_ids,
        CTA_testis_restricted_gene_names,
        CTA_unfiltered_gene_ids,
        CTA_unfiltered_gene_names,
        cta_symbol_for_alias,
    )
    from .ms_attribution import attribute_ms_evidence, summarize_ms_attribution
    from .ms_evidence import annotate_ms_detectability, cta_healthy_tissue_ms_hits
    from .mutations import HOTSPOT_MUTATIONS, mutant_peptides
    from .personalize import personalized_targets
    from .scoring import score_affinity, score_presentation
    from .spanning import spanning_pmhc_set
    from .species import (
        SpeciesContext,
        allele_prediction_capability,
        resolve_species_context,
    )
    from .targets import target_peptides
    from .tiers import (
        CONFIDENCE_VALUES,
        MS_RESTRICTION_VALUES,
        RESTRICTION_VALUES,
        RNA_RESTRICTION_LEVELS,
        SAFETY_TISSUE_GROUPS,
        VITAL_TISSUE_MS_NAMES,
    )
    from .tissues import (
        CORE_REPRODUCTIVE_TISSUES,
        EXTENDED_REPRODUCTIVE_TISSUES,
        HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS,
        PERMISSIVE_REPRODUCTIVE_TISSUES,
    )
    from .vaccine import design_vaccine
    from .vaccine_construct import VaccineConfig
    from .viral import orfdata_viral_peptides, orfdata_viral_subsets, viral_peptides

# Public names resolve lazily (PEP 562) so ``import tsarina`` -- and therefore
# ``tsarina --version`` -- does not pay for pandas/oncoref until a name is used.
# `personalized_targets` is the canonical name for the headline workflow (not
# `personalize`, which would shadow the `tsarina.personalize` submodule).
_LAZY = {
    "CONFIDENCE_VALUES": (".tiers", "CONFIDENCE_VALUES"),
    "CORE_REPRODUCTIVE_TISSUES": (".tissues", "CORE_REPRODUCTIVE_TISSUES"),
    "CTA_TIERS": (".cta_selection", "CTA_TIERS"),
    "CTA_by_axes": (".gene_sets", "CTA_by_axes"),
    "CTA_by_tier": (".gene_sets", "CTA_by_tier"),
    "CTA_cancer_antigen_gene_ids": (".gene_sets", "CTA_cancer_antigen_gene_ids"),
    "CTA_cancer_antigen_gene_names": (".gene_sets", "CTA_cancer_antigen_gene_names"),
    "CTA_clinical_target_gene_ids": (".gene_sets", "CTA_clinical_target_gene_ids"),
    "CTA_clinical_target_gene_names": (".gene_sets", "CTA_clinical_target_gene_names"),
    "CTA_detailed_evidence": (".evidence", "CTA_detailed_evidence"),
    "CTA_evidence": (".evidence", "CTA_evidence"),
    "CTA_excluded_gene_ids": (".gene_sets", "CTA_excluded_gene_ids"),
    "CTA_excluded_gene_names": (".gene_sets", "CTA_excluded_gene_names"),
    "CTA_filtered_gene_ids": (".gene_sets", "CTA_filtered_gene_ids"),
    "CTA_filtered_gene_names": (".gene_sets", "CTA_filtered_gene_names"),
    "CTA_gene_ids": (".gene_sets", "CTA_gene_ids"),
    "CTA_gene_names": (".gene_sets", "CTA_gene_names"),
    "CTA_hpa_evidence": ("oncoref.cta", "cta_evidence"),
    "CTA_never_expressed_gene_ids": (".gene_sets", "CTA_never_expressed_gene_ids"),
    "CTA_never_expressed_gene_names": (".gene_sets", "CTA_never_expressed_gene_names"),
    "CTA_placental_restricted_gene_ids": (".gene_sets", "CTA_placental_restricted_gene_ids"),
    "CTA_placental_restricted_gene_names": (".gene_sets", "CTA_placental_restricted_gene_names"),
    "CTA_relaxed_reproductive_gene_ids": (".gene_sets", "CTA_relaxed_reproductive_gene_ids"),
    "CTA_relaxed_reproductive_gene_names": (".gene_sets", "CTA_relaxed_reproductive_gene_names"),
    "CTA_reproductive_antigen_gene_ids": (".gene_sets", "CTA_reproductive_antigen_gene_ids"),
    "CTA_reproductive_antigen_gene_names": (".gene_sets", "CTA_reproductive_antigen_gene_names"),
    "CTA_selection": (".cta_selection", "CTA_selection"),
    "CTA_selection_coverage": (".cta_selection", "CTA_selection_coverage"),
    "CTA_sex_eligible_gene_ids": (".cta_selection", "CTA_sex_eligible_gene_ids"),
    "CTA_sex_eligible_gene_names": (".cta_selection", "CTA_sex_eligible_gene_names"),
    "CTA_sex_restricted_antigens": (".cta_selection", "CTA_sex_restricted_antigens"),
    "CTA_target_candidate_gene_ids": (".gene_sets", "CTA_target_candidate_gene_ids"),
    "CTA_target_candidate_gene_names": (".gene_sets", "CTA_target_candidate_gene_names"),
    "CTA_testis_restricted_gene_ids": (".gene_sets", "CTA_testis_restricted_gene_ids"),
    "CTA_testis_restricted_gene_names": (".gene_sets", "CTA_testis_restricted_gene_names"),
    "CTA_unfiltered_gene_ids": (".gene_sets", "CTA_unfiltered_gene_ids"),
    "CTA_unfiltered_gene_names": (".gene_sets", "CTA_unfiltered_gene_names"),
    "EXTENDED_REPRODUCTIVE_TISSUES": (".tissues", "EXTENDED_REPRODUCTIVE_TISSUES"),
    "HOTSPOT_MUTATIONS": (".mutations", "HOTSPOT_MUTATIONS"),
    "HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS": (".tissues", "HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS"),
    "IEDB27_AB": (".alleles", "IEDB27_AB"),
    "MS_RESTRICTION_VALUES": (".tiers", "MS_RESTRICTION_VALUES"),
    "PANEL_DEFINITIONS": (".alleles", "PANEL_DEFINITIONS"),
    "PERMISSIVE_REPRODUCTIVE_TISSUES": (".tissues", "PERMISSIVE_REPRODUCTIVE_TISSUES"),
    "RESTRICTION_VALUES": (".tiers", "RESTRICTION_VALUES"),
    "RNA_RESTRICTION_LEVELS": (".tiers", "RNA_RESTRICTION_LEVELS"),
    "SAFETY_TISSUE_GROUPS": (".tiers", "SAFETY_TISSUE_GROUPS"),
    "SpeciesContext": (".species", "SpeciesContext"),
    "VITAL_TISSUE_MS_NAMES": (".tiers", "VITAL_TISSUE_MS_NAMES"),
    "VaccineConfig": (".vaccine_construct", "VaccineConfig"),
    "allele_prediction_capability": (".species", "allele_prediction_capability"),
    "annotate_ms_detectability": (".ms_evidence", "annotate_ms_detectability"),
    "attribute_ms_evidence": (".ms_attribution", "attribute_ms_evidence"),
    "cancer_fusions": (".fusions", "cancer_fusions"),
    "cta_cancer_expression_features": (".cancer_expression", "cta_cancer_expression_features"),
    "cta_healthy_tissue_ms_hits": (".ms_evidence", "cta_healthy_tissue_ms_hits"),
    "cta_symbol_for_alias": (".gene_sets", "cta_symbol_for_alias"),
    "design_vaccine": (".vaccine", "design_vaccine"),
    "evidence_store_provenance": (".evidence_store", "evidence_store_provenance"),
    "fusion_junction_peptides": (".fusions", "fusion_junction_peptides"),
    "get_panel": (".alleles", "get_panel"),
    "hpa_cancer_ihc_prevalence": (".cancer_expression", "hpa_cancer_ihc_prevalence"),
    "hpa_cancer_rna_prevalence": (".cancer_expression", "hpa_cancer_rna_prevalence"),
    "load_evidence_store_ms_hits": (".evidence_store", "load_evidence_store_ms_hits"),
    "mutant_peptides": (".mutations", "mutant_peptides"),
    "orfdata_viral_peptides": (".viral", "orfdata_viral_peptides"),
    "orfdata_viral_subsets": (".viral", "orfdata_viral_subsets"),
    "panel_names": (".alleles", "panel_names"),
    "personalized_targets": (".personalize", "personalized_targets"),
    "resolve_species_context": (".species", "resolve_species_context"),
    "score_affinity": (".scoring", "score_affinity"),
    "score_presentation": (".scoring", "score_presentation"),
    "spanning_pmhc_set": (".spanning", "spanning_pmhc_set"),
    "summarize_ms_attribution": (".ms_attribution", "summarize_ms_attribution"),
    "target_peptides": (".targets", "target_peptides"),
    "viral_peptides": (".viral", "viral_peptides"),
}


def __getattr__(name: str):
    if name in _LAZY:
        module, attr = _LAZY[name]
        value = getattr(importlib.import_module(module, __name__), attr)
        globals()[name] = value
        return value
    try:  # submodules were implicitly importable as attributes before
        return importlib.import_module(f".{name}", __name__)
    except ModuleNotFoundError as e:
        if e.name != f"{__name__}.{name}":
            raise
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted({*globals(), *_LAZY})


__all__ = [
    "CONFIDENCE_VALUES",
    "CORE_REPRODUCTIVE_TISSUES",
    "CTA_TIERS",
    "EXTENDED_REPRODUCTIVE_TISSUES",
    "HOTSPOT_MUTATIONS",
    "HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS",
    "IEDB27_AB",
    "MS_RESTRICTION_VALUES",
    "PANEL_DEFINITIONS",
    "PERMISSIVE_REPRODUCTIVE_TISSUES",
    "RESTRICTION_VALUES",
    "RNA_RESTRICTION_LEVELS",
    "SAFETY_TISSUE_GROUPS",
    "VITAL_TISSUE_MS_NAMES",
    "CTA_by_axes",
    "CTA_by_tier",
    "CTA_cancer_antigen_gene_ids",
    "CTA_cancer_antigen_gene_names",
    "CTA_clinical_target_gene_ids",
    "CTA_clinical_target_gene_names",
    "CTA_detailed_evidence",
    "CTA_evidence",
    "CTA_excluded_gene_ids",
    "CTA_excluded_gene_names",
    "CTA_filtered_gene_ids",
    "CTA_filtered_gene_names",
    "CTA_gene_ids",
    "CTA_gene_names",
    "CTA_hpa_evidence",
    "CTA_never_expressed_gene_ids",
    "CTA_never_expressed_gene_names",
    "CTA_placental_restricted_gene_ids",
    "CTA_placental_restricted_gene_names",
    "CTA_relaxed_reproductive_gene_ids",
    "CTA_relaxed_reproductive_gene_names",
    "CTA_reproductive_antigen_gene_ids",
    "CTA_reproductive_antigen_gene_names",
    "CTA_selection",
    "CTA_selection_coverage",
    "CTA_sex_eligible_gene_ids",
    "CTA_sex_eligible_gene_names",
    "CTA_sex_restricted_antigens",
    "CTA_target_candidate_gene_ids",
    "CTA_target_candidate_gene_names",
    "CTA_testis_restricted_gene_ids",
    "CTA_testis_restricted_gene_names",
    "CTA_unfiltered_gene_ids",
    "CTA_unfiltered_gene_names",
    "SpeciesContext",
    "VaccineConfig",
    "__version__",
    "allele_prediction_capability",
    "annotate_ms_detectability",
    "attribute_ms_evidence",
    "cancer_fusions",
    "cta_cancer_expression_features",
    "cta_healthy_tissue_ms_hits",
    "cta_symbol_for_alias",
    "design_vaccine",
    "evidence_store_provenance",
    "fusion_junction_peptides",
    "get_panel",
    "hpa_cancer_ihc_prevalence",
    "hpa_cancer_rna_prevalence",
    "load_evidence_store_ms_hits",
    "mutant_peptides",
    "orfdata_viral_peptides",
    "orfdata_viral_subsets",
    "panel_names",
    "personalized_targets",
    "resolve_species_context",
    "score_affinity",
    "score_presentation",
    "spanning_pmhc_set",
    "summarize_ms_attribution",
    "target_peptides",
    "viral_peptides",
]
