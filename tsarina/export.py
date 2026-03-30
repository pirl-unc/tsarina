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

"""Bulk pMHC export pipeline and dual MS support map aggregation.

Builds per-peptide and per-(peptide, allele) mass spec evidence maps from
IEDB/CEDAR data, and exports comprehensive pMHC occurrence tables with
all evidence columns.

Typical usage::

    from tsarina.export import build_ms_support_maps, export_gene_properties

    pep_map, pmhc_map, summary_df = build_ms_support_maps(
        peptides=my_peptides,
        iedb_path="mhc_ligand_full.csv",
    )
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from hitlist.aggregate import aggregate_per_peptide, aggregate_per_pmhc
from hitlist.scanner import scan


def build_ms_support_maps(
    peptides: set[str],
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    mhc_class: str | None = "I",
) -> tuple[dict[str, dict], dict[tuple[str, str], dict], pd.DataFrame]:
    """Build dual peptide-level and pMHC-level MS evidence maps.

    Parameters
    ----------
    peptides
        Set of peptide sequences to look up.
    iedb_path
        Path to IEDB MHC ligand export.
    cedar_path
        Path to CEDAR MHC ligand export.
    mhc_class
        MHC class filter (default ``"I"``).

    Returns
    -------
    tuple[dict, dict, pd.DataFrame]
        - **peptide_map**: ``{peptide: {ms_hit_count, ms_alleles, found_in_cancer, ...}}``
        - **pmhc_map**: ``{(peptide, allele): {ms_hit_count, ms_ref_count, ...}}``
        - **summary_df**: Per-peptide summary DataFrame.
    """
    hits = scan(
        peptides=peptides,
        iedb_path=iedb_path,
        cedar_path=cedar_path,
        mhc_class=mhc_class,
        classify_source=True,
        human_only=False,
        hla_only=True,
    )

    # Per-peptide aggregation
    summary_df = aggregate_per_peptide(hits)
    peptide_map: dict[str, dict] = {}
    for row in summary_df.itertuples(index=False):
        peptide_map[row.peptide] = row._asdict()

    # Per-(peptide, allele) aggregation
    pmhc_df = aggregate_per_pmhc(hits)
    pmhc_map: dict[tuple[str, str], dict] = {}
    for row in pmhc_df.itertuples(index=False):
        pmhc_map[(row.peptide, row.mhc_restriction)] = row._asdict()

    return peptide_map, pmhc_map, summary_df


def export_gene_properties(
    gene_df: pd.DataFrame,
    output_path: str | Path,
    mtec_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Export a comprehensive gene property table.

    Parameters
    ----------
    gene_df
        Gene DataFrame (e.g. from ``CTA_evidence()`` or HPA-enriched).
    output_path
        Path to write CSV output.
    mtec_df
        Optional mTEC expression table to merge.

    Returns
    -------
    pd.DataFrame
        The exported DataFrame.
    """
    df = gene_df.copy()

    if mtec_df is not None:
        merge_col = "gene_name" if "gene_name" in mtec_df.columns else None
        if merge_col and merge_col in df.columns:
            mtec_cols = ["gene_name", "mean_gene_tpm", "threshold_bin", "mtec_le_1", "mtec_le_4"]
            available = [c for c in mtec_cols if c in mtec_df.columns]
            df = df.merge(
                mtec_df[available].drop_duplicates(subset=["gene_name"]),
                on="gene_name",
                how="left",
            )

    df.to_csv(str(output_path), index=False)
    return df
