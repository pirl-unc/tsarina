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

"""MS evidence curation -- delegates to hitlist.

All source classification, PMID overrides, and peptide evidence aggregation
are implemented in hitlist. This module re-exports for convenience.
"""

from hitlist.aggregate import aggregate_per_peptide as aggregate_peptide_evidence
from hitlist.aggregate import aggregate_per_pmhc
from hitlist.curation import (
    classify_ms_row,
    is_cancer_specific,
    load_pmid_overrides,
    load_tissue_categories,
)

PMID_OVERRIDES = load_pmid_overrides()

__all__ = [
    "PMID_OVERRIDES",
    "aggregate_peptide_evidence",
    "aggregate_per_pmhc",
    "classify_ms_row",
    "is_cancer_specific",
    "load_pmid_overrides",
    "load_tissue_categories",
]
