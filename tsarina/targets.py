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

    from tsarina.targets import target_peptides
    from tsarina.alleles import get_panel

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
    attach_ms_evidence: bool = False,
    classify_source: bool = True,
    cancer_specific: bool = False,
    require_human_exclusive_viral: bool = True,
    *,
    host_species: str = "human",
    mhc_species: str | None = None,
    requested_alleles: list[str] | tuple[str, ...] | None = None,
    evidence_scope: str = "all",
    observations: pd.DataFrame | None = None,
    custom_targets: pd.DataFrame | None = None,
    host_proteins: dict[str, str] | list[str] | None = None,
) -> pd.DataFrame:
    """Build a unified peptide table across CTA, viral, mutant, and custom targets.

    Supports species-scoped host/MHC references (#217) and restriction-aware
    MS attribution (#216). Refuses to silently transfer human CTA or human
    recurrent-mutation catalogs to non-human host species; pass
    ``cta=False, mutations=False`` with species-scoped *custom_targets* or
    *viruses* instead.
    """
    from .species import resolve_species_context

    ctx = resolve_species_context(
        host_species,
        mhc_species=mhc_species,
        ensembl_release=ensembl_release,
    )
    if not ctx.is_human_host:
        if cta:
            raise ValueError(
                f"Human CTA catalog cannot be silently transferred to host_species="
                f"{ctx.host_species!r}; set cta=False and pass species-scoped custom_targets."
            )
        if mutations:
            raise ValueError(
                f"Human HOTSPOT_MUTATIONS cannot be silently transferred to host_species="
                f"{ctx.host_species!r}; set mutations=False and pass species-scoped custom_targets."
            )

    frames: list[pd.DataFrame] = []

    if custom_targets is not None and not custom_targets.empty:
        req_cols = {"peptide", "category", "source"}
        missing = req_cols - set(custom_targets.columns)
        if missing:
            raise ValueError(f"custom_targets missing required columns: {sorted(missing)}")
        ct = custom_targets.copy()
        if "length" not in ct.columns:
            ct["length"] = ct["peptide"].astype(str).str.len()
        if "source_detail" not in ct.columns:
            ct["source_detail"] = ""
        frames.append(ct[["peptide", "length", "category", "source", "source_detail"]])

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
        from .viral import (
            ONCOGENIC_VIRUSES,
            host_exclusive_viral_peptides,
            viral_peptides,
        )

        virus_keys = sorted(ONCOGENIC_VIRUSES.keys()) if viruses is True else list(viruses)
        for vk in virus_keys:
            if require_human_exclusive_viral:
                if ctx.is_human_host and host_proteins is None:
                    from .viral import human_exclusive_viral_peptides

                    vdf = human_exclusive_viral_peptides(
                        virus=vk, lengths=lengths, ensembl_release=ensembl_release
                    )
                else:
                    vdf = host_exclusive_viral_peptides(
                        virus=vk,
                        lengths=lengths,
                        ensembl_release=ensembl_release,
                        host_species=host_species,
                        host_proteins=host_proteins,
                    )
            else:
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
        empty = pd.DataFrame(columns=["peptide", "length", "category", "source", "source_detail"])
        empty.attrs["species_context"] = ctx.to_provenance()
        return empty

    combined = pd.concat(frames, ignore_index=True)

    # Deduplicate: same peptide may come from multiple sources.
    # Keep all source annotations but mark unique peptides.
    combined["peptide_unique_id"] = combined["peptide"]

    # ── IEDB/CEDAR evidence lookup ──────────────────────────────────────
    needs_ms_evidence = attach_ms_evidence or require_ms_evidence or cancer_specific
    needs_ms_evidence = needs_ms_evidence or iedb_path is not None
    needs_ms_evidence = needs_ms_evidence or cedar_path is not None
    needs_ms_evidence = needs_ms_evidence or observations is not None
    needs_ms_evidence = needs_ms_evidence or requested_alleles is not None
    if needs_ms_evidence:
        from .ms_evidence import (
            SOURCE_FLAG_OUTPUTS,
            aggregate_ms_hits_by_peptide,
            load_public_ms_hits,
        )

        source_flag_outputs = SOURCE_FLAG_OUTPUTS if classify_source or cancer_specific else {}
        if observations is not None:
            hits = observations.copy()
        else:
            hits = load_public_ms_hits(
                peptides=set(combined["peptide"].unique()),
                iedb_path=iedb_path,
                cedar_path=cedar_path,
                mhc_class=mhc_class,
                mhc_species=ctx.mhc_species,
                classify_source=classify_source or cancer_specific,
                drop_binding_assays=True,
            )

        if requested_alleles is not None or evidence_scope != "all" or observations is not None:
            from .ms_attribution import attribute_ms_evidence, summarize_ms_attribution

            attributed = attribute_ms_evidence(
                combined[["peptide"]].drop_duplicates(),
                hits,
                requested_alleles=requested_alleles,
                species_context=ctx,
                evidence_scope=evidence_scope,
                evidence_loaded=True,
            )
            summary = summarize_ms_attribution(attributed, group_cols=("peptide",))
            combined = combined.merge(summary, on="peptide", how="left")
            combined["has_ms_evidence"] = combined["ms_hit_count"].notna() & (
                combined["ms_hit_count"] > 0
            )
            combined["ms_hit_count"] = combined["ms_hit_count"].fillna(0).astype(int)
            combined["ms_alleles"] = combined["ms_alleles"].fillna("")
            combined["ms_allele_count"] = combined["ms_allele_count"].fillna(0).astype(int)
            combined["ms_cancer"] = combined["ms_in_cancer"].fillna(False).astype(bool)
            combined["ms_healthy_tissue"] = (
                combined["ms_in_healthy_tissue"].fillna(False).astype(bool)
            )
            if evidence_scope != "all":
                combined = combined[combined["has_ms_evidence"]].reset_index(drop=True)
        else:
            hit_agg = aggregate_ms_hits_by_peptide(
                hits,
                source_flag_outputs=source_flag_outputs,
                cell_line_output="ms_cell_lines" if source_flag_outputs else None,
            )
            combined = combined.merge(hit_agg, on="peptide", how="left")
            combined["has_ms_evidence"] = combined["ms_hit_count"].notna() & (
                combined["ms_hit_count"] > 0
            )
            combined["ms_hit_count"] = combined["ms_hit_count"].fillna(0).astype(int)
            combined["ms_alleles"] = combined["ms_alleles"].fillna("")
            combined["ms_allele_count"] = combined["ms_allele_count"].fillna(0).astype(int)

            for col in source_flag_outputs.values():
                if col in combined.columns:
                    combined[col] = combined[col].where(combined[col].notna(), False).astype(bool)
            if "ms_cell_lines" in combined.columns:
                combined["ms_cell_lines"] = combined["ms_cell_lines"].fillna("")

        # ── Filters ─────────────────────────────────────────────────────
        if require_ms_evidence:
            combined = combined[combined["has_ms_evidence"]].reset_index(drop=True)

        if cancer_specific and "ms_cancer" in combined.columns:
            combined = combined[
                combined["ms_cancer"] & ~combined.get("ms_healthy_tissue", False)
            ].reset_index(drop=True)

    combined = combined.drop(columns=["peptide_unique_id"], errors="ignore")
    combined.attrs["species_context"] = ctx.to_provenance()
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
