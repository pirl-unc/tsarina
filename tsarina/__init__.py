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

"""tsarina: personalized cancer immunotherapy target selection from curated shared antigen data."""

from .alleles import (
    IEDB27_AB,
    PANEL_DEFINITIONS,
    get_panel,
    panel_names,
)
from .cancer_expression import (
    cta_cancer_expression_features,
    hpa_cancer_ihc_prevalence,
    hpa_cancer_rna_prevalence,
)
from .evidence import CTA_detailed_evidence, CTA_evidence
from .gene_sets import (
    CTA_by_axes,
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
    CTA_testis_restricted_gene_ids,
    CTA_testis_restricted_gene_names,
    CTA_unfiltered_gene_ids,
    CTA_unfiltered_gene_names,
    cta_symbol_for_alias,
)
from .ms_evidence import cta_healthy_tissue_ms_hits
from .mutations import HOTSPOT_MUTATIONS, mutant_peptides

# NB: `personalize` (the function) is intentionally NOT re-exported here — that
# would shadow the `tsarina.personalize` submodule. Import it from its module:
# ``from tsarina.personalize import personalize``.
from .scoring import score_affinity, score_presentation
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
from .version import __version__
from .viral import viral_peptides

__all__ = [
    "CONFIDENCE_VALUES",
    "CORE_REPRODUCTIVE_TISSUES",
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
    "CTA_detailed_evidence",
    "CTA_evidence",
    "CTA_excluded_gene_ids",
    "CTA_excluded_gene_names",
    "CTA_filtered_gene_ids",
    "CTA_filtered_gene_names",
    "CTA_gene_ids",
    "CTA_gene_names",
    "CTA_never_expressed_gene_ids",
    "CTA_never_expressed_gene_names",
    "CTA_placental_restricted_gene_ids",
    "CTA_placental_restricted_gene_names",
    "CTA_relaxed_reproductive_gene_ids",
    "CTA_relaxed_reproductive_gene_names",
    "CTA_testis_restricted_gene_ids",
    "CTA_testis_restricted_gene_names",
    "CTA_unfiltered_gene_ids",
    "CTA_unfiltered_gene_names",
    "__version__",
    "cta_cancer_expression_features",
    "cta_healthy_tissue_ms_hits",
    "cta_symbol_for_alias",
    "get_panel",
    "hpa_cancer_ihc_prevalence",
    "hpa_cancer_rna_prevalence",
    "mutant_peptides",
    "panel_names",
    "score_affinity",
    "score_presentation",
    "target_peptides",
    "viral_peptides",
]
