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

from pathlib import Path

import pandas as pd

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
        "ms_alleles": ("mhc_restriction", lambda x: ";".join(sorted({v for v in x if v}))),
        "ms_allele_count": ("mhc_restriction", lambda x: len({v for v in x if v})),
    }
    for source_col, output_col in source_flag_outputs.items():
        if source_col in hits.columns:
            agg[output_col] = (source_col, "any")
    if cell_line_output is not None and "cell_line_name" in hits.columns:
        agg[cell_line_output] = (
            "cell_line_name",
            lambda x: ";".join(sorted({v for v in x if v})),
        )

    return hits.groupby("peptide", as_index=False).agg(**agg)
