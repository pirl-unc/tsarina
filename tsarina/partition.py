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

"""Current CTA member genes and all other coding genes in a sequence reference.

The legacy never-expressed partition stays empty under the current policy.
Other/unassigned/Cancer-Antigen-only/Reproductive-Antigen-only selection
decisions remain in the non-CTA background.

Alternate-haplotype Ensembl annotations on ``HSCHR*`` contigs that represent
a curated CTA gene (either via OncoRef's ``alternate_ensembl_gene_ids`` or an
exact non-empty symbol match to a curated primary-contig CTA gene) are excluded
from ``non_cta`` so they do not self-veto their primary CTA gene (#187), while
keeping ``cta`` equal to the curated CTA gene IDs. Ordinary non-CTA loci on
primary chromosomes always remain in ``non_cta``.

Requires the ``pyensembl`` package (install with ``pip install tsarina[partition]``).
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from dataclasses import dataclass, field
from functools import lru_cache

import pandas as pd

from .evidence import CTA_evidence
from .gene_sets import CTA_gene_ids, is_coding_gene

ALIAS_COLUMNS = (
    "gene_id",
    "gene_name",
    "contig",
    "cta_reference_gene_ids",
    "resolution_basis",
)


@dataclass(frozen=True)
class CTAPartitionSets:
    """Three-way partition of protein-coding genes as sets.

    Attributes
    ----------
    cta : set[str]
        Current selected member genes. Source of CTA pMHCs.
    cta_never_expressed : set[str]
        Empty compatibility field; current categories do not imply no expression.
    non_cta : set[str]
        All other coding genes from the requested sequence reference, including
        Unsafe, Other, Cancer Antigen only, Reproductive Antigen only, unassigned
        and genes outside the selection table, excluding audited ``HSCHR*``
        alternate-haplotype annotations of curated CTA genes.
    cta_alternate_haplotypes : set[str]
        Audited ``HSCHR*`` alternate-haplotype gene IDs (or symbols) resolved to
        a curated primary-contig CTA gene and removed from ``non_cta``.
    """

    cta: set
    cta_never_expressed: set
    non_cta: set
    cta_alternate_haplotypes: set = field(default_factory=set)


@dataclass(frozen=True)
class CTAPartitionDataFrames:
    """Three-way partition of protein-coding genes as DataFrames.

    Attributes
    ----------
    cta : pd.DataFrame
        Current selected members with full evidence columns.
    cta_never_expressed : pd.DataFrame
        Empty compatibility frame retaining the evidence columns.
    non_cta : pd.DataFrame
        All other protein-coding genes (Symbol, Ensembl_Gene_ID).
    """

    cta: pd.DataFrame
    cta_never_expressed: pd.DataFrame
    non_cta: pd.DataFrame


@lru_cache(maxsize=1)
def _oncoref_alternate_id_to_canonical() -> dict[str, set[str]]:
    """Map OncoRef ``alternate_ensembl_gene_ids`` back to canonical gene IDs."""
    from oncoref import gene_identity_catalog

    catalog = gene_identity_catalog()
    mapping: dict[str, set[str]] = {}
    if "alternate_ensembl_gene_ids" not in catalog.columns:
        return mapping
    for row in catalog[["canonical_gene_id", "alternate_ensembl_gene_ids"]].itertuples(index=False):
        raw = row.alternate_ensembl_gene_ids
        if not isinstance(raw, str) or not raw.strip():
            continue
        try:
            alts = json.loads(raw) if raw.strip().startswith("[") else [raw.strip()]
        except json.JSONDecodeError:
            continue
        canonical = str(row.canonical_gene_id)
        for alt in alts:
            if alt:
                mapping.setdefault(str(alt), set()).add(canonical)
    return mapping


def _gene_attr(gene: object, primary: str, fallback: str) -> str:
    value = getattr(gene, primary, None)
    if value is None:
        value = getattr(gene, fallback, "")
    return "" if value is None else str(value).strip()


def _hschr_chromosome(contig: str) -> str:
    """Return the primary chromosome token from an ``HSCHR<chrom>_...`` contig name."""
    if not contig.startswith("HSCHR"):
        return ""
    body = contig[5:]
    if "_" not in body:
        return ""
    return body.split("_", 1)[0].strip()


def resolve_cta_alternate_haplotype_ids(
    genes: Iterable[object],
    cta_ids: set[str],
    *,
    catalog_alternates: dict[str, set[str]] | None = None,
) -> tuple[set[str], pd.DataFrame]:
    """Identify audited ``HSCHR*`` alternate-haplotype annotations of curated CTAs (#187).

    A coding gene outside ``cta_ids`` is resolved as an alternate-haplotype CTA
    occurrence only when its contig starts with ``HSCHR`` AND either:

    1. its Ensembl gene ID is explicitly listed in OncoRef's
       ``alternate_ensembl_gene_ids`` for a gene in ``cta_ids``, or
    2. its non-empty gene symbol matches a curated non-HSCHR gene in ``cta_ids``
       on the corresponding chromosome (when the ``HSCHR<chrom>_...`` contig
       names a chromosome and the primary CTA records a contig).

    Primary-chromosome loci and unnamed genes (empty symbol) are never matched.
    """
    gene_list = list(genes)
    if catalog_alternates is None:
        all_alts = _oncoref_alternate_id_to_canonical()
        alt_to_cta = {alt_id: refs & cta_ids for alt_id, refs in all_alts.items() if refs & cta_ids}
    else:
        alt_to_cta = {
            alt_id: set(refs) & cta_ids
            for alt_id, refs in catalog_alternates.items()
            if set(refs) & cta_ids
        }

    primary_by_name: dict[str, list[tuple[str, str]]] = {}
    for gene in gene_list:
        gid = _gene_attr(gene, "gene_id", "id")
        gname = _gene_attr(gene, "gene_name", "name")
        contig = _gene_attr(gene, "contig", "seqname")
        if gid in cta_ids and gname and not contig.startswith("HSCHR"):
            primary_by_name.setdefault(gname, []).append((gid, contig))

    aliases: dict[str, dict] = {}
    for gene in gene_list:
        gid = _gene_attr(gene, "gene_id", "id")
        gname = _gene_attr(gene, "gene_name", "name")
        contig = _gene_attr(gene, "contig", "seqname")
        if gid in cta_ids or not contig.startswith("HSCHR"):
            continue
        if gid in alt_to_cta:
            refs = alt_to_cta[gid]
            basis = "oncoref_alternate_id"
        elif gname and gname in primary_by_name:
            hschr_chrom = _hschr_chromosome(contig)
            matched_refs = {
                ref_gid
                for ref_gid, ref_contig in primary_by_name[gname]
                if not hschr_chrom or not ref_contig or ref_contig == hschr_chrom
            }
            if not matched_refs:
                continue
            refs = matched_refs
            basis = "hschr_same_symbol"
        else:
            continue
        aliases[gid] = {
            "gene_id": gid,
            "gene_name": gname,
            "contig": contig,
            "cta_reference_gene_ids": ";".join(sorted(refs)),
            "resolution_basis": basis,
        }

    table = pd.DataFrame(
        [aliases[key] for key in sorted(aliases)],
        columns=list(ALIAS_COLUMNS),
    )
    return set(aliases.keys()), table


def _build_partition(ensembl_release: int = 112):
    """Shared logic for building the three-way partition."""
    from pyensembl import EnsemblRelease

    ensembl = EnsemblRelease(ensembl_release)
    evidence_df = CTA_evidence()

    coding_genes = [g for g in ensembl.genes() if is_coding_gene(g)]
    all_pc_genes = {g.gene_id: g.gene_name for g in coding_genes}
    all_pc_ids = set(all_pc_genes.keys())

    evidence_ids = evidence_df["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
    cta_mask = evidence_ids.isin(CTA_gene_ids())
    # Current Other includes diffuse expression; null means unassigned.
    # Neither category justifies removing a gene from the sequence background.
    never_expressed_mask = pd.Series(False, index=evidence_df.index)

    cta_ids = set(evidence_df.loc[cta_mask, "Ensembl_Gene_ID"])
    never_expressed_ids = set(evidence_df.loc[never_expressed_mask, "Ensembl_Gene_ID"])
    alt_haplotype_ids, alt_aliases = resolve_cta_alternate_haplotype_ids(coding_genes, cta_ids)
    non_cta_ids = all_pc_ids - cta_ids - alt_haplotype_ids - never_expressed_ids

    return (
        all_pc_genes,
        evidence_df,
        cta_mask,
        never_expressed_mask,
        cta_ids,
        never_expressed_ids,
        non_cta_ids,
        alt_haplotype_ids,
        alt_aliases,
    )


def CTA_partition_alternate_haplotype_aliases(ensembl_release: int = 112) -> pd.DataFrame:
    """Return the audited ``HSCHR*`` alternate-haplotype CTA genes removed from ``non_cta``."""
    _, _, _, _, _, _, _, _, alt_aliases = _build_partition(ensembl_release)
    return alt_aliases


def CTA_partition_gene_ids(ensembl_release: int = 112) -> CTAPartitionSets:
    """Partition all protein-coding genes into CTA / never-expressed / non-CTA
    as sets of Ensembl gene IDs.

    All nonselected coding genes remain in ``non_cta``, except audited ``HSCHR*``
    alternate-haplotype annotations of curated CTA genes (recorded in
    ``cta_alternate_haplotypes``). The historical ``cta_never_expressed`` field
    is empty; Other/unassigned are not removed. Selection reference receipts are
    distinct from ``ensembl_release``.

    Examples
    --------
    >>> p = CTA_partition_gene_ids()
    >>> "ENSG00000147381" in p.cta   # MAGEA4
    True
    >>> len(p.cta & p.non_cta)       # no overlap
    0
    """
    _, _, _, _, cta_ids, never_expressed_ids, non_cta_ids, alt_ids, _ = _build_partition(
        ensembl_release
    )
    return CTAPartitionSets(
        cta=cta_ids,
        cta_never_expressed=never_expressed_ids,
        non_cta=non_cta_ids,
        cta_alternate_haplotypes=alt_ids,
    )


def CTA_partition_gene_names(ensembl_release: int = 112) -> CTAPartitionSets:
    """Partition all protein-coding genes into CTA / never-expressed / non-CTA
    as sets of gene symbols.

    All nonselected coding genes remain in ``non_cta``, except audited ``HSCHR*``
    alternate-haplotype annotations of curated CTA genes. The historical
    ``cta_never_expressed`` field is empty; Other/unassigned are not removed.
    Selection reference receipts are distinct from ``ensembl_release``.

    Examples
    --------
    >>> p = CTA_partition_gene_names()
    >>> "MAGEA4" in p.cta
    True
    >>> "TP53" in p.non_cta
    True
    """
    all_pc_genes, evidence_df, cta_mask, never_expressed_mask, _, _, non_cta_ids, alt_ids, _ = (
        _build_partition(ensembl_release)
    )

    cta_names = set(evidence_df.loc[cta_mask, "Symbol"])
    never_expressed_names = set(evidence_df.loc[never_expressed_mask, "Symbol"])
    non_cta_names = {all_pc_genes[gene_id] for gene_id in non_cta_ids}
    alt_names = {all_pc_genes[gene_id] for gene_id in alt_ids if gene_id in all_pc_genes}

    return CTAPartitionSets(
        cta=cta_names,
        cta_never_expressed=never_expressed_names,
        non_cta=non_cta_names,
        cta_alternate_haplotypes=alt_names,
    )


def CTA_partition_dataframes(ensembl_release: int = 112) -> CTAPartitionDataFrames:
    """Partition all protein-coding genes into CTA / never-expressed / non-CTA
    as DataFrames.

    The ``cta`` and ``cta_never_expressed`` DataFrames include all CTA
    evidence columns.  The ``non_cta`` DataFrame has Symbol and
    Ensembl_Gene_ID columns.

    All nonselected coding genes remain in ``non_cta``, except audited ``HSCHR*``
    alternate-haplotype annotations of curated CTA genes. The historical
    ``cta_never_expressed`` field is empty; Other/unassigned are not removed.
    Selection reference receipts are distinct from ``ensembl_release``.

    Examples
    --------
    >>> p = CTA_partition_dataframes()
    >>> "rna_deflated_reproductive_frac" in p.cta.columns
    True
    """
    all_pc_genes, evidence_df, cta_mask, never_expressed_mask, _, _, non_cta_ids, _, _ = (
        _build_partition(ensembl_release)
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
