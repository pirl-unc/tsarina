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

"""Greedy pMHC-guided gene panel selection.

Iteratively selects genes that maximize expected peptide-MHC coverage
across a patient cohort, with multi-level tie-breaking by MS evidence
strength, expression level, and mTEC expression.

Typical usage::

    from tsarina.selection import greedy_select_genes

    winners = greedy_select_genes(
        candidate_df=my_candidates,
        n_genes=10,
    )
"""

from __future__ import annotations

import pandas as pd


def greedy_select_genes(
    candidate_df: pd.DataFrame,
    n_genes: int = 10,
    coverage_col: str = "candidate_pmhc_count",
    ms_col: str | None = "ms_hit_count",
    expression_col: str | None = "source_tpm",
    mtec_col: str | None = "mean_mtec_tpm",
    priority_genes: list[str] | None = None,
) -> pd.DataFrame:
    """Greedily select genes to maximize pMHC coverage.

    At each step, selects the gene providing the most new pMHC
    candidates, with multi-level tie-breaking.

    Parameters
    ----------
    candidate_df
        DataFrame with at least ``gene`` and ``peptide`` columns,
        plus optional scoring columns. One row per gene-peptide pair.
    n_genes
        Number of genes to select (default 10).
    coverage_col
        Column with pMHC candidate count per gene (computed if missing).
    ms_col
        Column for MS evidence count tie-breaking. None to skip.
    expression_col
        Column for expression level tie-breaking. None to skip.
    mtec_col
        Column for mTEC TPM tie-breaking (lower is better). None to skip.
    priority_genes
        Ordered list of genes to prefer in ties.

    Returns
    -------
    pd.DataFrame
        Selected genes in order, with columns: ``gene``, ``step``,
        ``new_peptides`` (peptides added at this step),
        ``cumulative_peptides``.
    """
    if candidate_df.empty or n_genes <= 0:
        return pd.DataFrame(columns=["gene", "step", "new_peptides", "cumulative_peptides"])

    # Build gene -> set of peptides
    gene_peptides: dict[str, set[str]] = {}
    for gene, group in candidate_df.groupby("gene"):
        gene_peptides[str(gene)] = set(group["peptide"])

    # Build gene stats for tie-breaking
    gene_stats: dict[str, dict] = {}
    for gene, group in candidate_df.groupby("gene"):
        stats: dict = {"gene": str(gene)}
        if ms_col and ms_col in group.columns:
            stats["ms_score"] = group[ms_col].sum()
        if expression_col and expression_col in group.columns:
            stats["expression"] = group[expression_col].max()
        if mtec_col and mtec_col in group.columns:
            stats["mtec_tpm"] = group[mtec_col].min()
        gene_stats[str(gene)] = stats

    priority_rank = {}
    if priority_genes:
        for i, g in enumerate(priority_genes):
            priority_rank[g] = i

    # Greedy selection
    selected: list[dict] = []
    covered_peptides: set[str] = set()
    remaining_genes = set(gene_peptides.keys())

    for step in range(1, n_genes + 1):
        if not remaining_genes:
            break

        best_gene = None
        best_key = None

        for gene in remaining_genes:
            new_peps = gene_peptides[gene] - covered_peptides
            n_new = len(new_peps)

            # Tie-breaking key: (new_peptides DESC, ms_score DESC,
            # expression DESC, mtec_tpm ASC, priority ASC, gene ASC)
            stats = gene_stats.get(gene, {})
            key = (
                -n_new,
                -stats.get("ms_score", 0),
                -stats.get("expression", 0),
                stats.get("mtec_tpm", 999),
                priority_rank.get(gene, 999),
                gene,
            )

            if best_gene is None or key < best_key:
                best_gene = gene
                best_key = key

        if best_gene is None:
            break

        new_peps = gene_peptides[best_gene] - covered_peptides
        covered_peptides |= new_peps
        remaining_genes.remove(best_gene)

        selected.append(
            {
                "gene": best_gene,
                "step": step,
                "new_peptides": len(new_peps),
                "cumulative_peptides": len(covered_peptides),
            }
        )

    return pd.DataFrame(selected)


def region_weighted_selection(
    candidate_df: pd.DataFrame,
    region_allele_weights: dict[str, dict[str, float]],
    n_genes: int = 10,
    **kwargs,
) -> dict[str, pd.DataFrame]:
    """Run greedy gene selection per region and as a global mixture.

    Parameters
    ----------
    candidate_df
        DataFrame with ``gene``, ``peptide``, ``allele`` columns.
    region_allele_weights
        Dict mapping region name to {allele: weight} dicts.
    n_genes
        Number of genes to select per region.
    **kwargs
        Passed to :func:`greedy_select_genes`.

    Returns
    -------
    dict[str, pd.DataFrame]
        One selection result per region key, plus ``"global_mixture"``.
    """
    results: dict[str, pd.DataFrame] = {}

    for region, weights in region_allele_weights.items():
        region_alleles = set(weights.keys())
        region_df = candidate_df[candidate_df["allele"].isin(region_alleles)].copy()
        if not region_df.empty:
            results[region] = greedy_select_genes(region_df, n_genes=n_genes, **kwargs)

    # Global mixture: use all alleles
    results["global_mixture"] = greedy_select_genes(candidate_df, n_genes=n_genes, **kwargs)

    return results
