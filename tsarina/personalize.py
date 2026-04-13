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

"""Patient-level personalization: from tumor characteristics to prioritized targets.

The main entry point for clinical use.  Given a patient's HLA type, tumor
CTA expression, detected mutations, and viral status, produces a ranked
list of peptide-MHC targets with public MS evidence and predicted
presentation scores.

Typical usage::

    from tsarina.personalize import personalize

    targets = personalize(
        hla_alleles=["HLA-A*02:01", "HLA-B*07:02"],
        cta_expression={"MAGEA4": 142.5, "PRAME": 87.3},
        mutations=["KRAS G12D"],
        viruses=["hpv16"],
        iedb_path="mhc_ligand_full.csv",
    )
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def personalize(
    hla_alleles: list[str],
    cta_expression: dict[str, float] | None = None,
    mutations: list[str] | None = None,
    viruses: list[str] | None = None,
    lengths: tuple[int, ...] = (8, 9, 10, 11),
    ensembl_release: int = 112,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    mhc_class: str = "I",
    min_cta_tpm: float = 1.0,
    score_presentation: bool = True,
    skip_ms_evidence: bool = False,
    predictor: str = "mhcflurry",
) -> pd.DataFrame:
    """Build a personalized target list for a single patient.

    Parameters
    ----------
    hla_alleles
        Patient's HLA class I alleles (e.g. ``["HLA-A*02:01", "HLA-B*07:02"]``).
    cta_expression
        Dict mapping CTA gene symbol to tumor RNA expression in TPM.
        Only CTAs above ``min_cta_tpm`` are included.  If None, all
        expressed CTAs are considered (no expression filter).
    mutations
        List of mutation labels to match against the hotspot list
        (e.g. ``["KRAS G12D", "TP53 R175H"]``).  If None, no mutant
        peptides are included.
    viruses
        List of virus keys (e.g. ``["hpv16", "ebv"]``).
        If None, no viral peptides are included.
    lengths
        Peptide lengths (default 8-11).
    ensembl_release
        Ensembl release (default 112).
    iedb_path
        Path to IEDB MHC ligand export.  If None (default), auto-resolves
        from the hitlist data registry (register with
        ``tsarina data register iedb /path/to/mhc_ligand_full.csv``).
    cedar_path
        Path to CEDAR MHC ligand export.  Optional; auto-resolves if
        registered, otherwise silently skipped.
    mhc_class
        MHC class filter for IEDB scanning (default ``"I"``).
    min_cta_tpm
        Minimum CTA expression in TPM to include (default 1.0).
    score_presentation
        If True (default), score peptide-allele pairs for predicted
        presentation via topiary + mhctools.
    skip_ms_evidence
        If True, do not look up IEDB/CEDAR evidence.  Useful for dry
        runs when the datasets are not registered.
    predictor
        Which mhctools predictor to use for scoring.  One of
        ``"mhcflurry"`` (default), ``"netmhcpan"``, ``"netmhcpan_el"``.

    Returns
    -------
    pd.DataFrame
        Prioritized target list with columns:

        - ``peptide``: peptide sequence
        - ``length``: peptide length
        - ``category``: ``"cta"``, ``"viral"``, or ``"mutant"``
        - ``source``: gene name / virus name / mutation label
        - ``source_tpm``: RNA expression in TPM (CTAs only, NaN otherwise)
        - ``ms_hit_count``: number of IEDB/CEDAR MS observations
        - ``ms_alleles``: MHC restrictions observed in public data
        - ``ms_in_cancer``: detected in cancer samples
        - ``ms_in_healthy_tissue``: detected in normal non-reproductive tissue (safety flag)
        - ``best_allele``: patient HLA allele with best predicted presentation
        - ``presentation_percentile``: MHCflurry percentile for best allele
        - ``priority_score``: composite prioritization score (higher = better)
    """
    frames: list[pd.DataFrame] = []

    # ── CTA peptides ────────────────────────────────────────────────────
    if cta_expression is not None:
        from .gene_sets import CTA_gene_names
        from .peptides import cta_peptides

        valid_ctas = CTA_gene_names()
        expressed_ctas = {
            gene: tpm
            for gene, tpm in cta_expression.items()
            if gene in valid_ctas and tpm >= min_cta_tpm
        }
        if expressed_ctas:
            all_cta_peps = cta_peptides(ensembl_release=ensembl_release, lengths=lengths)
            cta_peps = all_cta_peps[all_cta_peps["gene_name"].isin(expressed_ctas)].copy()
            if not cta_peps.empty:
                cta_peps["source_tpm"] = cta_peps["gene_name"].map(expressed_ctas)
                frames.append(
                    pd.DataFrame(
                        {
                            "peptide": cta_peps["peptide"],
                            "length": cta_peps["length"],
                            "category": "cta",
                            "source": cta_peps["gene_name"],
                            "source_detail": cta_peps["gene_id"],
                            "source_tpm": cta_peps["source_tpm"],
                        }
                    )
                )
    elif cta_expression is None:
        # Include all expressed CTAs when no expression dict provided
        pass  # caller must explicitly provide expression to include CTAs

    # ── Mutant peptides ─────────────────────────────────────────────────
    if mutations:
        from .mutations import HOTSPOT_MUTATIONS
        from .mutations import mutant_peptides as _mutant_peptides

        # Match patient mutations against hotspot list
        mutation_labels = set(mutations)
        matched = [m for m in HOTSPOT_MUTATIONS if m["label"] in mutation_labels]
        if matched:
            mdf = _mutant_peptides(
                mutations=matched, lengths=lengths, ensembl_release=ensembl_release
            )
            if not mdf.empty:
                frames.append(
                    pd.DataFrame(
                        {
                            "peptide": mdf["peptide"],
                            "length": mdf["length"],
                            "category": "mutant",
                            "source": mdf["label"],
                            "source_detail": mdf["mutation"],
                            "source_tpm": float("nan"),
                        }
                    )
                )

    # ── Viral peptides ──────────────────────────────────────────────────
    if viruses:
        from .viral import viral_peptides

        for vk in viruses:
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
                            "source_tpm": float("nan"),
                        }
                    )
                )

    if not frames:
        return pd.DataFrame(
            columns=[
                "peptide",
                "length",
                "category",
                "source",
                "source_detail",
                "source_tpm",
            ]
        )

    combined = pd.concat(frames, ignore_index=True)

    # ── IEDB/CEDAR evidence ─────────────────────────────────────────────
    if not skip_ms_evidence:
        from .indexing import load_ms_evidence

        # iedb_path / cedar_path are legacy; if provided, we still pass them
        # through to a direct scan_public_ms call for users who want to
        # bypass the cached index.
        if iedb_path is not None or cedar_path is not None:
            from .datasources import resolve_dataset_paths
            from .iedb import scan_public_ms

            resolved_iedb, resolved_cedar = resolve_dataset_paths(iedb_path, cedar_path)
            hits = scan_public_ms(
                peptides=set(combined["peptide"].unique()),
                iedb_path=resolved_iedb,
                cedar_path=resolved_cedar,
                mhc_class=mhc_class,
                classify_source=True,
            )
        else:
            hits = load_ms_evidence(
                peptides=set(combined["peptide"].unique()),
                mhc_class=mhc_class,
                mhc_species="Homo sapiens",
            )

        if not hits.empty:
            agg_cols = {
                "ms_hit_count": ("peptide", "size"),
                "ms_alleles": ("mhc_restriction", lambda x: ";".join(sorted(set(x)))),
                "ms_allele_count": ("mhc_restriction", lambda x: len(set(x))),
            }
            if "src_cancer" in hits.columns:
                agg_cols["ms_in_cancer"] = ("src_cancer", "any")
            if "src_healthy_tissue" in hits.columns:
                agg_cols["ms_in_healthy_tissue"] = ("src_healthy_tissue", "any")

            hit_agg = hits.groupby("peptide", as_index=False).agg(**agg_cols)
            combined = combined.merge(hit_agg, on="peptide", how="left")
        else:
            combined["ms_hit_count"] = 0
            combined["ms_alleles"] = ""
            combined["ms_allele_count"] = 0
            combined["ms_in_cancer"] = False
            combined["ms_in_healthy_tissue"] = False

        for col, default in [
            ("ms_hit_count", 0),
            ("ms_alleles", ""),
            ("ms_allele_count", 0),
            ("ms_in_cancer", False),
            ("ms_in_healthy_tissue", False),
        ]:
            if col in combined.columns:
                combined[col] = combined[col].fillna(default)

    # ── Topiary presentation scoring ────────────────────────────────────
    if score_presentation and hla_alleles:
        try:
            from .scoring import score_presentation as _score

            unique_peps = combined["peptide"].unique().tolist()
            scores = _score(peptides=unique_peps, alleles=hla_alleles, predictor=predictor)

            # For each peptide, find the best (lowest percentile) allele
            if not scores.empty and "presentation_percentile" in scores.columns:
                best = (
                    scores.sort_values("presentation_percentile")
                    .groupby("peptide", as_index=False)
                    .first()
                    .rename(
                        columns={
                            "allele": "best_allele",
                            "presentation_percentile": "presentation_percentile",
                        }
                    )[["peptide", "best_allele", "presentation_percentile"]]
                )
                combined = combined.merge(best, on="peptide", how="left")
        except ImportError:
            pass  # MHCflurry not installed, skip scoring

    # ── Priority score ──────────────────────────────────────────────────
    # Composite score: MS evidence + expression + presentation
    combined["priority_score"] = 0.0
    if "ms_hit_count" in combined.columns:
        combined["priority_score"] += combined["ms_hit_count"].clip(upper=10) * 10
    if "ms_in_cancer" in combined.columns:
        combined["priority_score"] += combined["ms_in_cancer"].astype(float) * 20
    if "ms_in_healthy_tissue" in combined.columns:
        combined["priority_score"] -= combined["ms_in_healthy_tissue"].astype(float) * 50
    if "source_tpm" in combined.columns:
        combined["priority_score"] += combined["source_tpm"].fillna(0).clip(upper=500) * 0.1
    if "presentation_percentile" in combined.columns:
        # Lower percentile = better, so invert
        combined["priority_score"] += (100 - combined["presentation_percentile"].fillna(100)).clip(
            lower=0
        ) * 0.5

    combined = combined.sort_values("priority_score", ascending=False).reset_index(drop=True)
    return combined
