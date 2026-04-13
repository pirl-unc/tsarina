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

"""Negative peptide sets from normal (non-tumor) tissue mass spectrometry.

Peptides observed on healthy somatic tissue MHC molecules by mass spec
represent off-target toxicity risks -- targeting them would damage normal
tissue.  This module extracts these negative sets from IEDB/CEDAR data,
distinguishing:

- **Somatic tissue**: the primary safety exclusion set
- **Thymus**: expected for CTAs (AIRE-mediated), not a safety concern
- **Reproductive tissue**: expected for CTAs, not a safety concern

Typical usage::

    from tsarina.negatives import healthy_tissue_peptides

    # Peptides confirmed on normal somatic tissue (the danger set)
    somatic = healthy_tissue_peptides(
        iedb_path="mhc_ligand_full.csv",
        tissue_type="somatic",
    )

    # Check overlap with your targets
    my_targets = {"SLYNTVATL", "GILGFVFTL"}
    dangerous = my_targets & somatic
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def healthy_tissue_peptides(
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    tissue_type: str = "somatic",
    mhc_class: str | None = "I",
    ex_vivo_only: bool = True,
) -> set[str]:
    """Extract peptides observed on healthy tissue MHC by mass spectrometry.

    Parameters
    ----------
    iedb_path
        Path to IEDB MHC ligand export.
    cedar_path
        Path to CEDAR MHC ligand export.
    tissue_type
        Which tissue category to extract:

        - ``"somatic"`` (default): non-reproductive, non-thymus healthy tissue.
          This is the primary safety exclusion set.
        - ``"thymus"``: thymus tissue only (expected for CTAs).
        - ``"reproductive"``: reproductive tissue only (expected for CTAs).
        - ``"all_healthy"``: all healthy tissue regardless of type.
    mhc_class
        MHC class filter (default ``"I"``).
    ex_vivo_only
        If True (default), only include Direct Ex Vivo observations
        (highest confidence tissue evidence).

    Returns
    -------
    set[str]
        Set of peptide sequences.
    """
    valid_types = ("somatic", "thymus", "reproductive", "all_healthy")
    if tissue_type not in valid_types:
        raise ValueError(f"tissue_type must be one of {valid_types}, got {tissue_type!r}")

    from .iedb import profile_dataset

    df = profile_dataset(iedb_path=iedb_path, cedar_path=cedar_path)

    if df.empty:
        return set()

    if mhc_class is not None:
        df = df[df["mhc_class"] == mhc_class]

    if ex_vivo_only:
        df = df[df["src_ex_vivo"]]

    if tissue_type == "somatic":
        df = df[df["src_healthy_tissue"]]
    elif tissue_type == "thymus":
        df = df[df["src_healthy_thymus"]]
    elif tissue_type == "reproductive":
        df = df[df["src_healthy_reproductive"]]
    elif tissue_type == "all_healthy":
        df = df[
            df["src_healthy_tissue"] | df["src_healthy_thymus"] | df["src_healthy_reproductive"]
        ]

    return set(df["peptide"].dropna().unique())


def healthy_tissue_peptide_details(
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    mhc_class: str | None = "I",
) -> pd.DataFrame:
    """Extract detailed healthy tissue peptide observations.

    Returns a DataFrame with per-peptide aggregated information about
    which tissues, alleles, and source contexts each peptide was
    observed in.

    Parameters
    ----------
    iedb_path
        Path to IEDB MHC ligand export.
    cedar_path
        Path to CEDAR MHC ligand export.
    mhc_class
        MHC class filter (default ``"I"``).

    Returns
    -------
    pd.DataFrame
        Columns: ``peptide``, ``tissue_types`` (semicolon-separated),
        ``alleles``, ``allele_count``, ``hit_count``,
        ``in_somatic``, ``in_thymus``, ``in_reproductive``.
    """
    from .iedb import profile_dataset

    df = profile_dataset(iedb_path=iedb_path, cedar_path=cedar_path)

    if df.empty:
        return pd.DataFrame(
            columns=[
                "peptide",
                "tissue_types",
                "alleles",
                "allele_count",
                "hit_count",
                "in_somatic",
                "in_thymus",
                "in_reproductive",
            ]
        )

    if mhc_class is not None:
        df = df[df["mhc_class"] == mhc_class]

    # Filter to any healthy tissue
    healthy = df[
        df["src_healthy_tissue"] | df["src_healthy_thymus"] | df["src_healthy_reproductive"]
    ]

    if healthy.empty:
        return pd.DataFrame(
            columns=[
                "peptide",
                "tissue_types",
                "alleles",
                "allele_count",
                "hit_count",
                "in_somatic",
                "in_thymus",
                "in_reproductive",
            ]
        )

    result = healthy.groupby("peptide", as_index=False).agg(
        tissue_types=(
            "source_tissue",
            lambda x: ";".join(sorted({v for v in x if v and str(v) != "nan"})),
        ),
        alleles=(
            "mhc_restriction",
            lambda x: ";".join(sorted({v for v in x if v and str(v) != "nan"})),
        ),
        allele_count=("mhc_restriction", lambda x: len({v for v in x if v and str(v) != "nan"})),
        hit_count=("peptide", "size"),
        in_somatic=("src_healthy_tissue", "any"),
        in_thymus=("src_healthy_thymus", "any"),
        in_reproductive=("src_healthy_reproductive", "any"),
    )

    return result.sort_values("hit_count", ascending=False).reset_index(drop=True)
