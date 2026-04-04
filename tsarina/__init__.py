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
    CTA_testis_restricted_gene_ids,
    CTA_testis_restricted_gene_names,
    CTA_unfiltered_gene_ids,
    CTA_unfiltered_gene_names,
)
from .tiers import (
    CONFIDENCE_VALUES,
    MS_RESTRICTION_VALUES,
    RESTRICTION_VALUES,
    RNA_RESTRICTION_LEVELS,
)
from .tissues import (
    CORE_REPRODUCTIVE_TISSUES,
    EXTENDED_REPRODUCTIVE_TISSUES,
    HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS,
    PERMISSIVE_REPRODUCTIVE_TISSUES,
)
from .version import __version__

__all__ = [
    "CONFIDENCE_VALUES",
    "CORE_REPRODUCTIVE_TISSUES",
    "EXTENDED_REPRODUCTIVE_TISSUES",
    "HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS",
    "IEDB27_AB",
    "MS_RESTRICTION_VALUES",
    "PANEL_DEFINITIONS",
    "PERMISSIVE_REPRODUCTIVE_TISSUES",
    "RESTRICTION_VALUES",
    "RNA_RESTRICTION_LEVELS",
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
    "CTA_testis_restricted_gene_ids",
    "CTA_testis_restricted_gene_names",
    "CTA_unfiltered_gene_ids",
    "CTA_unfiltered_gene_names",
    "__version__",
    "get_panel",
    "panel_names",
]
