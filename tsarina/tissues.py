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

"""Tissue utilities over oncoref's canonical CTA tissue definitions."""

from __future__ import annotations

from oncoref import cta_tissues as _oncoref_tissues

CORE_REPRODUCTIVE_TISSUES = _oncoref_tissues.CORE_REPRODUCTIVE_TISSUES
EXTENDED_REPRODUCTIVE_TISSUES = _oncoref_tissues.EXTENDED_REPRODUCTIVE_TISSUES
HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS = _oncoref_tissues.HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS
HPA_EXPRESSION_FLOOR_NTPM = _oncoref_tissues.HPA_EXPRESSION_FLOOR_NTPM
PERMISSIVE_REPRODUCTIVE_TISSUES = _oncoref_tissues.PERMISSIVE_REPRODUCTIVE_TISSUES
PROTEIN_RELIABILITY_ORDER = _oncoref_tissues.PROTEIN_RELIABILITY_ORDER
adaptive_rna_threshold = _oncoref_tissues.adaptive_rna_threshold

# This marker-enrichment heuristic belongs to Tsarina's ad-hoc HPA annotation
# helper; it is not part of canonical CTA membership.
HPA_MARKER_STRICT_MAX_RNA_TISSUES: int = 7


def is_tissue_restricted(
    tissues: set[str] | frozenset[str],
    allowed: frozenset[str] = CORE_REPRODUCTIVE_TISSUES,
    exclude_thymus: bool = True,
) -> bool:
    """Whether detected tissues are confined to an allowed reproductive set."""
    check = tissues - {"thymus"} if exclude_thymus else tissues
    return check.issubset(allowed)


def nonreproductive_tissues(
    tissues: set[str] | frozenset[str],
    definition: frozenset[str] = CORE_REPRODUCTIVE_TISSUES,
    exclude_thymus: bool = True,
) -> frozenset[str]:
    """Return detected tissues outside a reproductive-tissue definition."""
    excluded = tissues - definition
    if exclude_thymus:
        excluded -= {"thymus"}
    return frozenset(excluded)
