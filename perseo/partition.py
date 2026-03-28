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

"""Three-way partition of all protein-coding genes into CTA / never-expressed / non-CTA.

Requires the ``pyensembl`` package (install with ``pip install perseo[partition]``).
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from .evidence import CTA_evidence


@dataclass(frozen=True)
class CTAPartitionSets:
    """Three-way partition of protein-coding genes as sets.

    Attributes
    ----------
    cta : set[str]
        Expressed, reproductive-restricted CTAs.  Source of CTA pMHCs.
    cta_never_expressed : set[str]
        CTAs from source databases with no meaningful HPA expression
        (no protein data + max RNA < 2 nTPM).
    non_cta : set[str]
        All other protein-coding genes, including CTAs that fail the
        reproductive-tissue filter (somatic expression detected).
    """

    cta: set
    cta_never_expressed: set
    non_cta: set


@dataclass(frozen=True)
class CTAPartitionDataFrames:
    """Three-way partition of protein-coding genes as DataFrames.

    Attributes
    ----------
    cta : pd.DataFrame
        Expressed, reproductive-restricted CTAs with full evidence columns.
    cta_never_expressed : pd.DataFrame
        Never-expressed CTAs with full evidence columns.
    non_cta : pd.DataFrame
        All other protein-coding genes (Symbol, Ensembl_Gene_ID).
    """

    cta: pd.DataFrame
    cta_never_expressed: pd.DataFrame
    non_cta: pd.DataFrame


def _build_partition(ensembl_release: int = 112):
    """Shared logic for building the three-way partition."""
    from pyensembl import EnsemblRelease

    ensembl = EnsemblRelease(ensembl_release)
    evidence_df = CTA_evidence()

    all_pc_genes = {
        g.gene_id: g.gene_name for g in ensembl.genes() if g.biotype == "protein_coding"
    }
    all_pc_ids = set(all_pc_genes.keys())

    filtered_mask = evidence_df["filtered"].astype(str).str.lower() == "true"
    never_expr_mask = evidence_df["never_expressed"].astype(str).str.lower() == "true"

    cta_mask = filtered_mask & ~never_expr_mask
    never_expressed_mask = filtered_mask & never_expr_mask

    cta_ids = set(evidence_df.loc[cta_mask, "Ensembl_Gene_ID"])
    never_expressed_ids = set(evidence_df.loc[never_expressed_mask, "Ensembl_Gene_ID"])
    non_cta_ids = all_pc_ids - cta_ids - never_expressed_ids

    return (
        all_pc_genes,
        evidence_df,
        cta_mask,
        never_expressed_mask,
        cta_ids,
        never_expressed_ids,
        non_cta_ids,
    )


def CTA_partition_gene_ids(ensembl_release: int = 112) -> CTAPartitionSets:
    """Partition all protein-coding genes into CTA / never-expressed / non-CTA
    as sets of Ensembl gene IDs.

    CTAs that fail the reproductive-tissue filter (somatic expression)
    are included in ``non_cta``.

    Examples
    --------
    >>> p = CTA_partition_gene_ids()
    >>> "ENSG00000147381" in p.cta   # MAGEA4
    True
    >>> len(p.cta & p.non_cta)       # no overlap
    0
    """
    _, _, _, _, cta_ids, never_expressed_ids, non_cta_ids = _build_partition(ensembl_release)
    return CTAPartitionSets(
        cta=cta_ids,
        cta_never_expressed=never_expressed_ids,
        non_cta=non_cta_ids,
    )


def CTA_partition_gene_names(ensembl_release: int = 112) -> CTAPartitionSets:
    """Partition all protein-coding genes into CTA / never-expressed / non-CTA
    as sets of gene symbols.

    CTAs that fail the reproductive-tissue filter (somatic expression)
    are included in ``non_cta``.

    Examples
    --------
    >>> p = CTA_partition_gene_names()
    >>> "MAGEA4" in p.cta
    True
    >>> "TP53" in p.non_cta
    True
    """
    all_pc_genes, evidence_df, cta_mask, never_expressed_mask, _, _, _ = _build_partition(
        ensembl_release
    )
    all_pc_names = set(all_pc_genes.values())

    cta_names = set(evidence_df.loc[cta_mask, "Symbol"])
    never_expressed_names = set(evidence_df.loc[never_expressed_mask, "Symbol"])
    non_cta_names = all_pc_names - cta_names - never_expressed_names

    return CTAPartitionSets(
        cta=cta_names,
        cta_never_expressed=never_expressed_names,
        non_cta=non_cta_names,
    )


def CTA_partition_dataframes(ensembl_release: int = 112) -> CTAPartitionDataFrames:
    """Partition all protein-coding genes into CTA / never-expressed / non-CTA
    as DataFrames.

    The ``cta`` and ``cta_never_expressed`` DataFrames include all CTA
    evidence columns.  The ``non_cta`` DataFrame has Symbol and
    Ensembl_Gene_ID columns.

    CTAs that fail the reproductive-tissue filter (somatic expression)
    are included in ``non_cta``.

    Examples
    --------
    >>> p = CTA_partition_dataframes()
    >>> "rna_deflated_reproductive_frac" in p.cta.columns
    True
    """
    all_pc_genes, evidence_df, cta_mask, never_expressed_mask, _, _, non_cta_ids = _build_partition(
        ensembl_release
    )

    non_cta_records = [
        {"Symbol": all_pc_genes[gid], "Ensembl_Gene_ID": gid}
        for gid in sorted(non_cta_ids)
        if gid in all_pc_genes
    ]

    return CTAPartitionDataFrames(
        cta=evidence_df.loc[cta_mask].copy().reset_index(drop=True),
        cta_never_expressed=evidence_df.loc[never_expressed_mask].copy().reset_index(drop=True),
        non_cta=pd.DataFrame(non_cta_records),
    )
