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

"""Medullary thymic epithelial cell (mTEC) expression filtering.

mTECs express tissue-restricted antigens via AIRE for central immune
tolerance. CTAs with low mTEC expression (<=1 TPM) are better
immunotherapy targets because T cells reactive to them are less
likely to have been deleted during thymic selection.

Typical usage::

    from tsarina.mtec import load_mtec_gene_table, filter_by_mtec

    mtec = load_mtec_gene_table("mTEC-quants.gene_tpm_matrix.tsv")
    low_mtec_genes = mtec[mtec["mtec_le_1"]]["gene_name"].tolist()
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def load_mtec_gene_table(
    matrix_path: str | Path,
    gene_col: str = "gene_symbol",
) -> pd.DataFrame:
    """Load an mTEC gene TPM matrix and compute per-gene summary statistics.

    Parameters
    ----------
    matrix_path
        Path to a TSV with gene symbols as rows and mTEC samples as columns.
    gene_col
        Name of the gene symbol column (default ``"gene_symbol"``).

    Returns
    -------
    pd.DataFrame
        Columns: ``gene_name``, ``mean_gene_tpm``, ``median_gene_tpm``,
        ``stdev_gene_tpm``, ``mtec_le_1``, ``mtec_le_2``, ``mtec_le_4``,
        ``threshold_bin``.
    """
    matrix = pd.read_csv(str(matrix_path), sep="\t", low_memory=False)
    sample_cols = [c for c in matrix.columns if c != gene_col]

    gene_df = pd.DataFrame(
        {
            "gene_name": matrix[gene_col].astype(str),
            "mean_gene_tpm": matrix[sample_cols].mean(axis=1),
            "median_gene_tpm": matrix[sample_cols].median(axis=1),
            "stdev_gene_tpm": matrix[sample_cols].std(axis=1, ddof=1),
        }
    )

    gene_df["mtec_le_1"] = gene_df["mean_gene_tpm"] <= 1.0
    gene_df["mtec_le_2"] = gene_df["mean_gene_tpm"] <= 2.0
    gene_df["mtec_le_4"] = gene_df["mean_gene_tpm"] <= 4.0
    gene_df["threshold_bin"] = pd.cut(
        gene_df["mean_gene_tpm"],
        bins=[-np.inf, 1.0, 2.0, 4.0, np.inf],
        labels=["<=1 TPM", "(1,2] TPM", "(2,4] TPM", ">4 TPM"],
        right=True,
    ).astype(str)

    return gene_df.sort_values(["mean_gene_tpm", "gene_name"], ascending=[True, True]).reset_index(
        drop=True
    )


def filter_by_mtec(
    gene_names: set[str],
    mtec_df: pd.DataFrame,
    threshold: float = 1.0,
) -> set[str]:
    """Filter gene set to those with low mTEC expression.

    Parameters
    ----------
    gene_names
        Set of gene symbols to filter.
    mtec_df
        Output from :func:`load_mtec_gene_table`.
    threshold
        Maximum mean mTEC TPM (default 1.0).

    Returns
    -------
    set[str]
        Genes with mean mTEC TPM <= threshold.
    """
    low = mtec_df[mtec_df["mean_gene_tpm"] <= threshold]
    return gene_names & set(low["gene_name"])
