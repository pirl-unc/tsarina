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

"""Population-spanning source protein x HLA allele panel matrices.

Generates spreadsheet-style tables with source proteins as rows, HLA alleles
as columns, and peptide counts (or best presentation scores) at intersections.
These tables support cohort-level analysis of which targets are coverable by
which alleles across a population.

Typical usage::

    from tsarina.panels import build_panel_matrix
    from tsarina.alleles import get_panel

    # CTA gene x allele matrix with peptide counts
    matrix = build_panel_matrix(
        category="cta",
        alleles=get_panel("iedb36_abc"),
        metric="peptide_count",
    )

    # All categories, best presentation percentile per cell
    matrix = build_panel_matrix(
        category="all",
        alleles=get_panel("iedb27_ab"),
        metric="best_percentile",
        iedb_path="mhc_ligand_full.csv",
    )
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .alleles import IEDB27_AB


def build_panel_matrix(
    category: str = "all",
    alleles: list[str] | None = None,
    viruses: list[str] | None = None,
    metric: str = "peptide_count",
    lengths: tuple[int, ...] = (8, 9, 10, 11),
    ensembl_release: int = 112,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    mhc_class: str = "I",
    ms_confirmed_only: bool = False,
) -> pd.DataFrame:
    """Build a source protein x HLA allele panel matrix.

    Parameters
    ----------
    category
        Which target categories to include: ``"cta"``, ``"viral"``,
        ``"mutant"``, or ``"all"`` (default).
    alleles
        HLA allele panel. Defaults to IEDB-27.
    viruses
        Virus keys for viral category (default: all 9).
    metric
        What to put in each cell:

        - ``"peptide_count"``: number of unique peptides from this source
          predicted to bind this allele (requires MHCflurry)
        - ``"ms_peptide_count"``: number of MS-confirmed peptides
          (requires ``iedb_path``)
        - ``"best_percentile"``: best MHCflurry presentation percentile
          across all peptides for this source x allele (requires MHCflurry)
        - ``"has_peptide"``: boolean -- any peptide exists for this
          source x allele combination
    lengths
        Peptide lengths (default 8-11).
    ensembl_release
        Ensembl release (default 112).
    iedb_path
        Path to IEDB export (for MS evidence metrics).
    cedar_path
        Path to CEDAR export.
    mhc_class
        MHC class filter (default ``"I"``).
    ms_confirmed_only
        If True, only include peptides with IEDB/CEDAR MS evidence.

    Returns
    -------
    pd.DataFrame
        Rows = source proteins (gene name, virus, or mutation label),
        columns include ``source``, ``category``, plus one column per
        HLA allele with the selected metric.
    """
    valid_metrics = {"peptide_count", "ms_peptide_count", "best_percentile", "has_peptide"}
    if metric not in valid_metrics:
        raise ValueError(
            f"Unknown metric '{metric}'. Supported: {', '.join(sorted(valid_metrics))}."
        )

    if alleles is None:
        alleles = list(IEDB27_AB)

    from .targets import target_peptides

    include_cta = category in ("cta", "all")
    include_viral = category in ("viral", "all")
    include_mutant = category in ("mutant", "all")

    viral_targets: list[str] | bool
    if include_viral:
        viral_targets = True if viruses is None else viruses
    else:
        viral_targets = False

    df = target_peptides(
        cta=include_cta,
        viruses=viral_targets,
        mutations=include_mutant,
        lengths=lengths,
        ensembl_release=ensembl_release,
        iedb_path=iedb_path,
        cedar_path=cedar_path,
        mhc_class=mhc_class,
        require_ms_evidence=ms_confirmed_only,
        attach_ms_evidence=metric == "ms_peptide_count",
    )

    if df.empty:
        result = pd.DataFrame(columns=["source", "category", *alleles])
        return result

    # Get unique sources
    sources = df[["source", "category"]].drop_duplicates().sort_values(["category", "source"])

    if metric in ("peptide_count", "best_percentile"):
        # Need MHCflurry scoring; import/setup failures must surface because
        # allele-agnostic fallback values look valid but are scientifically wrong.
        from .scoring import score_presentation

        unique_peps = df["peptide"].unique().tolist()
        scores = score_presentation(peptides=unique_peps, alleles=alleles)
        if scores is not None and not scores.empty:
            # Join scores back to source info
            pep_sources = df[["peptide", "source", "category"]].drop_duplicates()
            scored = scores.merge(pep_sources, on="peptide", how="inner")

            rows = []
            for _, src_row in sources.iterrows():
                src = src_row["source"]
                cat = src_row["category"]
                src_data = scored[(scored["source"] == src) & (scored["category"] == cat)]
                row: dict = {"source": src, "category": cat}
                for allele in alleles:
                    allele_data = src_data[src_data["allele"] == allele]
                    if metric == "peptide_count":
                        row[allele] = allele_data["peptide"].nunique()
                    elif metric == "best_percentile":
                        if (
                            not allele_data.empty
                            and "presentation_percentile" in allele_data.columns
                        ):
                            row[allele] = allele_data["presentation_percentile"].min()
                        else:
                            row[allele] = None
                rows.append(row)
            return pd.DataFrame(rows)

        empty_rows = []
        for _, src_row in sources.iterrows():
            row = {"source": src_row["source"], "category": src_row["category"]}
            for allele in alleles:
                row[allele] = 0 if metric == "peptide_count" else None
            empty_rows.append(row)
        return pd.DataFrame(empty_rows)

    # Simple metrics that don't need MHCflurry
    if metric == "ms_peptide_count" and "has_ms_evidence" in df.columns:
        from .mhc import mhc_restriction_matches_any, normalize_mhc_restriction_set

        ms_df = df[df["has_ms_evidence"]].copy()
        wanted_by_allele = {allele: normalize_mhc_restriction_set([allele]) for allele in alleles}
        rows = []
        for _, src_row in sources.iterrows():
            src = src_row["source"]
            cat = src_row["category"]
            src_peps = ms_df[(ms_df["source"] == src) & (ms_df["category"] == cat)]
            row = {"source": src, "category": cat}
            for allele in alleles:
                # Count peptides with MS evidence for this allele
                if "ms_alleles" in src_peps.columns:
                    wanted = wanted_by_allele[allele]
                    mask = src_peps["ms_alleles"].map(
                        lambda value, wanted=wanted: mhc_restriction_matches_any(value, wanted)
                    )
                    count = src_peps[mask]["peptide"].nunique()
                else:
                    count = 0
                row[allele] = count
            rows.append(row)
        return pd.DataFrame(rows)

    if metric == "ms_peptide_count":
        raise ValueError("metric='ms_peptide_count' requires MS evidence columns.")

    # Fallback: simple peptide presence per source (allele-agnostic).
    rows = []
    for _, src_row in sources.iterrows():
        src = src_row["source"]
        cat = src_row["category"]
        src_peps = df[(df["source"] == src) & (df["category"] == cat)]
        n = src_peps["peptide"].nunique()
        row = {"source": src, "category": cat}
        for allele in alleles:
            row[allele] = bool(n) if metric == "has_peptide" else n
        rows.append(row)
    return pd.DataFrame(rows)
