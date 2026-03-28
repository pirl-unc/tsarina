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

"""Unified target peptide selection across CTA, viral, and mutant categories.

Combines peptide generation from all three target categories with IEDB/CEDAR
mass spec evidence filtering, MHC allele restriction, and cancer-specificity
classification into a single query interface.

Typical usage::

    from perseo.targets import target_peptides
    from perseo.alleles import get_panel

    # All MS-confirmed, cancer-specific peptides across categories
    df = target_peptides(
        cta=True,
        viruses=["hpv16", "ebv"],
        mutations=True,
        iedb_path="mhc_ligand_full.csv",
        require_ms_evidence=True,
    )
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def target_peptides(
    cta: bool = True,
    viruses: list[str] | bool = False,
    mutations: bool = True,
    lengths: tuple[int, ...] = (8, 9, 10, 11),
    ensembl_release: int = 112,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    mhc_class: str | None = "I",
    require_ms_evidence: bool = False,
    classify_source: bool = True,
    cancer_specific: bool = False,
) -> pd.DataFrame:
    """Build a unified peptide table across CTA, viral, and mutant targets.

    Parameters
    ----------
    cta
        Include CTA-exclusive peptides (default True).
    viruses
        Include viral peptides. Pass True for all 9 oncogenic viruses,
        or a list of virus keys (e.g. ``["hpv16", "ebv"]``).
        Pass False to skip.
    mutations
        Include recurrent mutant-spanning peptides (default True).
    lengths
        Peptide lengths (default 8-11).
    ensembl_release
        Ensembl release (default 112).
    iedb_path
        Path to IEDB MHC ligand export for MS evidence lookup.
    cedar_path
        Path to CEDAR MHC ligand export.
    mhc_class
        MHC class filter for IEDB scanning (default ``"I"``).
    require_ms_evidence
        If True, only return peptides with at least one IEDB/CEDAR hit.
    classify_source
        If True (default), add source context columns from IEDB
        (``src_cancer``, ``src_healthy``, etc.).
    cancer_specific
        If True, filter to peptides with cancer MS evidence
        (``src_cancer=True``) and without healthy-tissue evidence
        (``src_healthy=False``).  Only applies when MS evidence is
        available.

    Returns
    -------
    pd.DataFrame
        Columns: ``peptide``, ``length``, ``category`` (``"cta"``,
        ``"viral"``, ``"mutant"``), ``source`` (gene name, virus name,
        or mutation label), plus IEDB evidence columns when
        ``iedb_path`` is provided.
    """
    frames: list[pd.DataFrame] = []

    # ── CTA peptides ────────────────────────────────────────────────────
    if cta:
        from .peptides import cta_exclusive_peptides

        cta_df = cta_exclusive_peptides(ensembl_release=ensembl_release, lengths=lengths)
        if not cta_df.empty:
            frames.append(
                pd.DataFrame(
                    {
                        "peptide": cta_df["peptide"],
                        "length": cta_df["length"],
                        "category": "cta",
                        "source": cta_df["gene_name"],
                        "source_detail": cta_df["gene_id"],
                    }
                )
            )

    # ── Viral peptides ──────────────────────────────────────────────────
    if viruses:
        from .viral import ONCOGENIC_VIRUSES, viral_peptides

        virus_keys = sorted(ONCOGENIC_VIRUSES.keys()) if viruses is True else list(viruses)
        for vk in virus_keys:
            vdf = viral_peptides(virus=vk, lengths=lengths)
            if not vdf.empty:
                frames.append(
                    pd.DataFrame(
                        {
                            "peptide": vdf["peptide"],
                            "length": vdf["length"],
                            "category": "viral",
                            "source": vdf["virus"],
                            "source_detail": vdf["protein_id"],
                        }
                    )
                )

    # ── Mutant peptides ─────────────────────────────────────────────────
    if mutations:
        from .mutations import mutant_peptides as _mutant_peptides

        mdf = _mutant_peptides(lengths=lengths, ensembl_release=ensembl_release)
        if not mdf.empty:
            frames.append(
                pd.DataFrame(
                    {
                        "peptide": mdf["peptide"],
                        "length": mdf["length"],
                        "category": "mutant",
                        "source": mdf["label"],
                        "source_detail": mdf["mutation"],
                    }
                )
            )

    if not frames:
        return pd.DataFrame(columns=["peptide", "length", "category", "source", "source_detail"])

    combined = pd.concat(frames, ignore_index=True)

    # Deduplicate: same peptide may come from multiple sources.
    # Keep all source annotations but mark unique peptides.
    combined["peptide_unique_id"] = combined["peptide"]

    # ── IEDB/CEDAR evidence lookup ──────────────────────────────────────
    has_iedb = iedb_path is not None or cedar_path is not None
    if has_iedb:
        from .iedb import scan_public_ms

        unique_peptides = set(combined["peptide"].unique())
        hits = scan_public_ms(
            peptides=unique_peptides,
            iedb_path=iedb_path,
            cedar_path=cedar_path,
            mhc_class=mhc_class,
            classify_source=classify_source,
            human_only=False,  # viral peptides may not be human-source
            hla_only=True,
        )

        if not hits.empty:
            # Aggregate per peptide
            _SRC_COLS = [
                "src_cancer",
                "src_healthy_tissue",
                "src_healthy_thymus",
                "src_healthy_reproductive",
                "src_cell_line",
                "src_ebv_lcl",
                "src_ex_vivo",
            ]
            src_agg = {}
            if classify_source:
                for sc in _SRC_COLS:
                    if sc in hits.columns:
                        src_agg[f"ms_{sc.removeprefix('src_')}"] = (sc, "any")
                # Collect cell line names
                if "cell_line_name" in hits.columns:
                    src_agg["ms_cell_lines"] = (
                        "cell_line_name",
                        lambda x: ";".join(sorted({v for v in x if v})),
                    )

            hit_agg = hits.groupby("peptide", as_index=False).agg(
                ms_hit_count=("peptide", "size"),
                ms_alleles=("mhc_restriction", lambda x: ";".join(sorted(set(x)))),
                ms_allele_count=("mhc_restriction", lambda x: len(set(x))),
                **src_agg,
            )

            combined = combined.merge(hit_agg, on="peptide", how="left")
            combined["has_ms_evidence"] = combined["ms_hit_count"].notna() & (
                combined["ms_hit_count"] > 0
            )
            combined["ms_hit_count"] = combined["ms_hit_count"].fillna(0).astype(int)
            combined["ms_alleles"] = combined["ms_alleles"].fillna("")
            combined["ms_allele_count"] = combined["ms_allele_count"].fillna(0).astype(int)

            if classify_source:
                for sc in _SRC_COLS:
                    col = f"ms_{sc.removeprefix('src_')}"
                    if col in combined.columns:
                        combined[col] = combined[col].fillna(False)
                if "ms_cell_lines" in combined.columns:
                    combined["ms_cell_lines"] = combined["ms_cell_lines"].fillna("")
        else:
            combined["has_ms_evidence"] = False
            combined["ms_hit_count"] = 0
            combined["ms_alleles"] = ""
            combined["ms_allele_count"] = 0

        # ── Filters ─────────────────────────────────────────────────────
        if require_ms_evidence:
            combined = combined[combined["has_ms_evidence"]].reset_index(drop=True)

        if cancer_specific and "ms_cancer" in combined.columns:
            combined = combined[
                combined["ms_cancer"] & ~combined.get("ms_healthy_tissue", False)
            ].reset_index(drop=True)

    combined = combined.drop(columns=["peptide_unique_id"], errors="ignore")
    return combined


def target_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Summarize a target peptide table by category.

    Parameters
    ----------
    df
        Output from :func:`target_peptides`.

    Returns
    -------
    pd.DataFrame
        One row per category with counts of unique peptides, sources,
        and MS-confirmed peptides.
    """
    rows = []
    for cat, gdf in df.groupby("category"):
        row = {
            "category": cat,
            "unique_peptides": gdf["peptide"].nunique(),
            "total_occurrences": len(gdf),
            "sources": gdf["source"].nunique(),
        }
        if "has_ms_evidence" in gdf.columns:
            row["ms_confirmed_peptides"] = gdf[gdf["has_ms_evidence"]]["peptide"].nunique()
        if "ms_cancer" in gdf.columns:
            row["cancer_ms_peptides"] = gdf[gdf["ms_cancer"]]["peptide"].nunique()
        if "ms_healthy_tissue" in gdf.columns:
            row["healthy_tissue_peptides"] = gdf[gdf["ms_healthy_tissue"]]["peptide"].nunique()
        if "ms_cell_lines" in gdf.columns:
            all_lines = set()
            for v in gdf["ms_cell_lines"].dropna():
                all_lines.update(cl for cl in str(v).split(";") if cl)
            row["unique_cell_lines"] = len(all_lines)
        rows.append(row)
    return pd.DataFrame(rows)
