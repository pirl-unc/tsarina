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

"""CTA evidence table access with HPA tissue-restriction columns.

The evidence table contains one row per CTA gene (358 genes from multiple
source databases), annotated with Human Protein Atlas v23 protein and RNA
tissue expression data and the three-axis tier results.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .loader import cta_dataframe


def CTA_evidence() -> pd.DataFrame:
    """Return the full CTA evidence DataFrame with HPA tissue-restriction columns.

    Columns
    -------
    Symbol, Aliases, Full_Name, Function, Ensembl_Gene_ID,
    source_databases, biotype, Canonical_Transcript_ID
        Gene identity fields.
    protein_reproductive : bool or "no data"
        True if all IHC-detected tissues (excluding thymus) are in
        {testis, ovary, placenta}.
    protein_thymus : bool or "no data"
        True if protein detected in thymus.
    rna_reproductive : bool
        True if every tissue with >=1 nTPM (excluding thymus) is in
        {testis, ovary, placenta}.
    rna_thymus : bool
        True if thymus nTPM >= 1.
    protein_reliability : str
        Best HPA antibody reliability for this gene: "Enhanced",
        "Supported", "Approved", "Uncertain", or "no data".
    protein_strict_expression : str
        Semicolon-separated list of tissues where protein is detected
        (excluding thymus), or "no data" / "not detected".
    rna_reproductive_frac : float
        Fraction of total nTPM (excluding thymus) in core reproductive
        tissues, computed from raw nTPM values.
    rna_reproductive_and_thymus_frac : float
        Same but with thymus nTPM added to numerator and denominator.
    rna_deflated_reproductive_frac : float
        (1 + repro_deflated) / (1 + total_deflated) where each tissue
        is deflated via max(0, nTPM - 1).  The +1 pseudocount prevents
        0/0 for very-low-expression genes.
    rna_deflated_reproductive_and_thymus_frac : float
        Same but with thymus deflated nTPM added to the reproductive
        numerator.
    rna_80_pct_filter, rna_90_pct_filter, rna_95_pct_filter : bool
        Whether deflated reproductive fraction >= 80/90/95%.
    filtered : bool
        Final inclusion flag with tiered RNA thresholds based on protein
        antibody reliability.
    rna_max_ntpm : float
        Maximum nTPM across all tissues for this gene.
    never_expressed : bool
        True if no HPA protein data AND maximum RNA nTPM < 2.
    protein_restriction : str
        Protein-based tissue restriction: ``TESTIS``, ``PLACENTAL``,
        ``OVARIAN``, or empty (no protein data).
    protein_testis, protein_ovary, protein_placenta : str
        Per-tissue IHC detection: ``"True"`` / ``"False"`` if protein
        data exists, empty string if no protein data.
    rna_restriction : str
        RNA-based tissue restriction: ``TESTIS``, ``PLACENTAL``,
        ``OVARIAN``, ``REPRODUCTIVE``, or empty.
    rna_restriction_level : str
        RNA restriction quality: ``STRICT`` / ``MODERATE`` / ``PERMISSIVE``
        or empty.
    rna_testis_ntpm, rna_ovary_ntpm, rna_placenta_ntpm : float
        Per-tissue RNA expression (nTPM) from HPA consensus.
    rna_max_somatic_tissue : str
        Highest-expressing non-reproductive tissue.
    rna_max_somatic_ntpm : float
        nTPM of that tissue.
    rna_somatic_detected_count : int
        Number of non-reproductive tissues with nTPM >= 1.
    ms_restriction : str
        MS-based restriction. Default ``NO_MS_DATA``; computed at
        runtime by :func:`CTA_detailed_evidence` when IEDB data available.
    restriction : str
        Synthesized restriction integrating protein + RNA + MS:
        ``TESTIS`` / ``PLACENTAL`` / ``OVARIAN`` / ``REPRODUCTIVE`` / empty.
    restriction_confidence : str
        Cross-modality confidence: ``HIGH`` / ``MODERATE`` / ``LOW`` / empty.
    """
    return cta_dataframe()


def CTA_detailed_evidence(
    hpa_bulk_path: str | Path | None = None,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    genes: set[str] | list[str] | None = None,
) -> pd.DataFrame:
    """Return CTA evidence with full per-tissue breakdown and MS safety.

    Augments the shipped CTA evidence table with optional detail columns
    computed from raw data files.  IEDB/CEDAR paths auto-resolve from
    the hitlist data registry when not provided.

    Parameters
    ----------
    hpa_bulk_path
        Path to HPA ``proteinatlas.tsv``.  If provided, adds per-tissue
        RNA nTPM columns: ``rna_testis_ntpm``, ``rna_ovary_ntpm``,
        ``rna_placenta_ntpm``, ``rna_max_somatic_tissue``,
        ``rna_max_somatic_ntpm``, ``rna_somatic_detected_count``.
    iedb_path
        Path to IEDB ``mhc_ligand_full.csv``.  If None, auto-resolves
        from hitlist data registry.
    cedar_path
        Path to CEDAR ``mhc-ligand-full.csv``.  If None, auto-resolves
        from hitlist data registry.
    genes
        Gene symbols to compute MS restriction for.  If None, uses all
        CTA genes (slower — enumerates peptides for all 358 genes).

    Returns
    -------
    pd.DataFrame
        CTA evidence DataFrame with additional detail columns.
    """
    df = cta_dataframe().copy()

    if hpa_bulk_path is not None:
        from .hpa import enrich_hpa_evidence, extract_per_tissue_detail, parse_ntpm_entries

        enriched = enrich_hpa_evidence(df, hpa_bulk_path=hpa_bulk_path)
        if "rna_tissue_ntpm_map" in enriched.columns:
            detail_rows = enriched["rna_tissue_ntpm_map"].map(extract_per_tissue_detail)
            detail_df = pd.DataFrame(detail_rows.tolist(), index=enriched.index)
            for col in detail_df.columns:
                df[col] = detail_df[col]
        elif "RNA tissue specific nTPM" in enriched.columns:
            maps = enriched["RNA tissue specific nTPM"].map(parse_ntpm_entries)
            detail_rows = maps.map(extract_per_tissue_detail)
            detail_df = pd.DataFrame(detail_rows.tolist(), index=enriched.index)
            for col in detail_df.columns:
                df[col] = detail_df[col]

    # MS restriction computation — auto-resolves paths from hitlist registry.
    # Requires hitlist and pyensembl; import errors are caught gracefully.
    import contextlib

    with contextlib.suppress(ImportError):
        df = _compute_ms_restriction(df, iedb_path, cedar_path, genes=genes)

    return df


def _compute_ms_restriction(
    df: pd.DataFrame,
    iedb_path: str | Path | None,
    cedar_path: str | Path | None,
    genes: set[str] | list[str] | None = None,
) -> pd.DataFrame:
    """Compute per-gene MS restriction from IEDB/CEDAR data."""
    import contextlib

    from hitlist.downloads import get_path
    from hitlist.scanner import scan

    from .tiers import aggregate_gene_ms_safety

    # Auto-resolve paths from hitlist registry if not provided
    if iedb_path is None:
        iedb_path = get_path("iedb")
    if cedar_path is None:
        with contextlib.suppress(KeyError):
            cedar_path = get_path("cedar")

    # Get CTA peptide sequences, filtered to genes of interest
    from .peptides import cta_peptides

    pep_df = cta_peptides()
    if pep_df.empty:
        return df

    if genes is not None:
        gene_set = set(genes)
        pep_df = pep_df[pep_df["gene_name"].isin(gene_set)]
        if pep_df.empty:
            return df

    peptide_set = set(pep_df["peptide"])

    hits = scan(
        peptides=peptide_set,
        iedb_path=str(iedb_path),
        cedar_path=str(cedar_path) if cedar_path else None,
        classify_source=True,
        human_only=True,
        mhc_class="I",
    )
    if hits.empty:
        return df

    gene_map = pep_df[["peptide", "gene_name"]].drop_duplicates()
    ms_df = aggregate_gene_ms_safety(hits, gene_map)

    if not ms_df.empty:
        ms_merge = ms_df.set_index("gene_name")
        for col in ms_merge.columns:
            default = (
                "NO_MS_DATA"
                if col == "ms_restriction"
                else ("" if ms_merge[col].dtype == object else 0)
            )
            df[col] = df["Symbol"].map(ms_merge[col]).fillna(default)

    return df
