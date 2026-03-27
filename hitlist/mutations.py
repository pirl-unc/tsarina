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

"""Recurrent cancer mutation-derived peptides.

Generates peptides spanning recurrent somatic hotspot mutations (KRAS G12,
BRAF V600E, TP53 hotspots, etc.) and cross-references with IEDB/CEDAR
mass spec data.

The mutation list focuses on high-frequency driver mutations that produce
shared neoantigens -- the same mutant peptide across many patients --
making them attractive off-the-shelf immunotherapy targets.

Typical usage::

    from hitlist.mutations import (
        HOTSPOT_MUTATIONS,
        mutant_peptides,
        mutant_iedb_overlap,
    )

    # Generate all mutant-spanning peptides
    df = mutant_peptides()

    # Check IEDB overlap
    hits = mutant_iedb_overlap(iedb_path="mhc_ligand_full.csv")
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .peptides import AA20

# ── Hotspot mutation definitions ────────────────────────────────────────────
#
# Each entry specifies:
#   gene: HGNC gene symbol
#   gene_id: Ensembl gene ID
#   transcript_id: Ensembl transcript ID (canonical)
#   protein_position: 1-based amino acid position in the canonical protein
#   ref_aa: wildtype amino acid (single letter)
#   alt_aa: mutant amino acid (single letter)
#   label: human-readable mutation name (e.g. "KRAS G12D")
#   cancer_types: list of cancer types where this mutation recurs
#   frequency_note: approximate frequency context

HOTSPOT_MUTATIONS: list[dict] = [
    # ── KRAS codon 12/13 ────────────────────────────────────────────────
    {
        "gene": "KRAS",
        "gene_id": "ENSG00000133703",
        "transcript_id": "ENST00000256078",
        "protein_position": 12,
        "ref_aa": "G",
        "alt_aa": "D",
        "label": "KRAS G12D",
        "cancer_types": ["pancreatic", "colorectal", "NSCLC"],
        "frequency_note": "Most common KRAS mutation overall (~30% of KRAS-mutant cancers)",
    },
    {
        "gene": "KRAS",
        "gene_id": "ENSG00000133703",
        "transcript_id": "ENST00000256078",
        "protein_position": 12,
        "ref_aa": "G",
        "alt_aa": "V",
        "label": "KRAS G12V",
        "cancer_types": ["pancreatic", "NSCLC"],
        "frequency_note": "~20% of KRAS-mutant cancers",
    },
    {
        "gene": "KRAS",
        "gene_id": "ENSG00000133703",
        "transcript_id": "ENST00000256078",
        "protein_position": 12,
        "ref_aa": "G",
        "alt_aa": "C",
        "label": "KRAS G12C",
        "cancer_types": ["NSCLC", "colorectal"],
        "frequency_note": "~13% of KRAS-mutant; sotorasib/adagrasib target",
    },
    {
        "gene": "KRAS",
        "gene_id": "ENSG00000133703",
        "transcript_id": "ENST00000256078",
        "protein_position": 12,
        "ref_aa": "G",
        "alt_aa": "R",
        "label": "KRAS G12R",
        "cancer_types": ["pancreatic"],
        "frequency_note": "~15% of pancreatic KRAS mutations",
    },
    {
        "gene": "KRAS",
        "gene_id": "ENSG00000133703",
        "transcript_id": "ENST00000256078",
        "protein_position": 13,
        "ref_aa": "G",
        "alt_aa": "D",
        "label": "KRAS G13D",
        "cancer_types": ["colorectal"],
        "frequency_note": "~13% of colorectal KRAS mutations",
    },
    # ── BRAF ─────────────────────────────────────────────────────────────
    {
        "gene": "BRAF",
        "gene_id": "ENSG00000157764",
        "transcript_id": "ENST00000646891",
        "protein_position": 600,
        "ref_aa": "V",
        "alt_aa": "E",
        "label": "BRAF V600E",
        "cancer_types": ["melanoma", "colorectal", "thyroid", "hairy cell leukemia"],
        "frequency_note": "~50% of melanomas, ~10% of colorectal",
    },
    {
        "gene": "BRAF",
        "gene_id": "ENSG00000157764",
        "transcript_id": "ENST00000646891",
        "protein_position": 600,
        "ref_aa": "V",
        "alt_aa": "K",
        "label": "BRAF V600K",
        "cancer_types": ["melanoma"],
        "frequency_note": "~5-6% of melanomas",
    },
    # ── TP53 hotspots ────────────────────────────────────────────────────
    {
        "gene": "TP53",
        "gene_id": "ENSG00000141510",
        "transcript_id": "ENST00000269305",
        "protein_position": 175,
        "ref_aa": "R",
        "alt_aa": "H",
        "label": "TP53 R175H",
        "cancer_types": ["breast", "colorectal", "ovarian", "many solid tumors"],
        "frequency_note": "Most common TP53 hotspot (~6% of TP53-mutant cancers)",
    },
    {
        "gene": "TP53",
        "gene_id": "ENSG00000141510",
        "transcript_id": "ENST00000269305",
        "protein_position": 248,
        "ref_aa": "R",
        "alt_aa": "W",
        "label": "TP53 R248W",
        "cancer_types": ["breast", "colorectal", "ovarian", "many solid tumors"],
        "frequency_note": "~4% of TP53-mutant cancers",
    },
    {
        "gene": "TP53",
        "gene_id": "ENSG00000141510",
        "transcript_id": "ENST00000269305",
        "protein_position": 273,
        "ref_aa": "R",
        "alt_aa": "H",
        "label": "TP53 R273H",
        "cancer_types": ["breast", "colorectal", "ovarian", "many solid tumors"],
        "frequency_note": "~4% of TP53-mutant cancers",
    },
    {
        "gene": "TP53",
        "gene_id": "ENSG00000141510",
        "transcript_id": "ENST00000269305",
        "protein_position": 245,
        "ref_aa": "G",
        "alt_aa": "S",
        "label": "TP53 G245S",
        "cancer_types": ["breast", "lung", "many solid tumors"],
        "frequency_note": "~3% of TP53-mutant cancers",
    },
    {
        "gene": "TP53",
        "gene_id": "ENSG00000141510",
        "transcript_id": "ENST00000269305",
        "protein_position": 249,
        "ref_aa": "R",
        "alt_aa": "S",
        "label": "TP53 R249S",
        "cancer_types": ["hepatocellular carcinoma"],
        "frequency_note": "Aflatoxin-associated HCC hotspot",
    },
    # ── PIK3CA ───────────────────────────────────────────────────────────
    {
        "gene": "PIK3CA",
        "gene_id": "ENSG00000121879",
        "transcript_id": "ENST00000263967",
        "protein_position": 1047,
        "ref_aa": "H",
        "alt_aa": "R",
        "label": "PIK3CA H1047R",
        "cancer_types": ["breast", "endometrial", "colorectal"],
        "frequency_note": "Most common PIK3CA hotspot",
    },
    {
        "gene": "PIK3CA",
        "gene_id": "ENSG00000121879",
        "transcript_id": "ENST00000263967",
        "protein_position": 545,
        "ref_aa": "E",
        "alt_aa": "K",
        "label": "PIK3CA E545K",
        "cancer_types": ["breast", "endometrial", "colorectal"],
        "frequency_note": "Second most common PIK3CA hotspot",
    },
    # ── IDH1 ─────────────────────────────────────────────────────────────
    {
        "gene": "IDH1",
        "gene_id": "ENSG00000138413",
        "transcript_id": "ENST00000415913",
        "protein_position": 132,
        "ref_aa": "R",
        "alt_aa": "H",
        "label": "IDH1 R132H",
        "cancer_types": ["glioma", "AML", "cholangiocarcinoma"],
        "frequency_note": "~90% of IDH1-mutant gliomas; peptidomics+ vaccine target (NOA-16)",
    },
    # ── NRAS ─────────────────────────────────────────────────────────────
    {
        "gene": "NRAS",
        "gene_id": "ENSG00000213281",
        "transcript_id": "ENST00000369535",
        "protein_position": 61,
        "ref_aa": "Q",
        "alt_aa": "R",
        "label": "NRAS Q61R",
        "cancer_types": ["melanoma", "AML"],
        "frequency_note": "Most common NRAS mutation in melanoma",
    },
    {
        "gene": "NRAS",
        "gene_id": "ENSG00000213281",
        "transcript_id": "ENST00000369535",
        "protein_position": 61,
        "ref_aa": "Q",
        "alt_aa": "K",
        "label": "NRAS Q61K",
        "cancer_types": ["melanoma", "thyroid"],
        "frequency_note": "Second most common NRAS Q61 mutation",
    },
    # ── EGFR ─────────────────────────────────────────────────────────────
    {
        "gene": "EGFR",
        "gene_id": "ENSG00000146648",
        "transcript_id": "ENST00000275493",
        "protein_position": 858,
        "ref_aa": "L",
        "alt_aa": "R",
        "label": "EGFR L858R",
        "cancer_types": ["NSCLC"],
        "frequency_note": "~40% of EGFR-mutant NSCLC",
    },
    {
        "gene": "EGFR",
        "gene_id": "ENSG00000146648",
        "transcript_id": "ENST00000275493",
        "protein_position": 790,
        "ref_aa": "T",
        "alt_aa": "M",
        "label": "EGFR T790M",
        "cancer_types": ["NSCLC"],
        "frequency_note": "Resistance mutation; osimertinib target",
    },
]


def mutant_peptides(
    mutations: list[dict] | None = None,
    lengths: tuple[int, ...] = (8, 9, 10, 11),
    ensembl_release: int = 112,
    flank_length: int = 15,
) -> pd.DataFrame:
    """Generate mutant-spanning peptides from recurrent cancer hotspot mutations.

    For each mutation, retrieves the wildtype protein sequence from Ensembl,
    introduces the point mutation, and generates all k-mer peptides that
    span the mutated position.  Only peptides that differ from the wildtype
    (i.e., contain the mutant residue) are returned.

    Requires ``pyensembl`` (install with ``pip install hitlist[peptides]``).

    Parameters
    ----------
    mutations
        List of mutation dicts (see :data:`HOTSPOT_MUTATIONS`).
        If None, uses all hotspot mutations.
    lengths
        Peptide lengths to generate (default 8-11).
    ensembl_release
        Ensembl release for protein sequences (default 112).
    flank_length
        Flanking residues to include (default 15).

    Returns
    -------
    pd.DataFrame
        Columns: ``label``, ``gene``, ``gene_id``, ``mutation``,
        ``ref_aa``, ``alt_aa``, ``protein_position``, ``peptide``,
        ``wildtype_peptide``, ``length``, ``start``, ``end``,
        ``n_flank``, ``c_flank``, ``mutation_offset`` (0-based position
        of the mutant residue within the peptide).
    """
    from pyensembl import EnsemblRelease

    if mutations is None:
        mutations = HOTSPOT_MUTATIONS

    ensembl = EnsemblRelease(ensembl_release)
    rows: list[dict] = []

    for mut in mutations:
        transcript_id = mut["transcript_id"]
        pos = mut["protein_position"]  # 1-based
        ref_aa = mut["ref_aa"]
        alt_aa = mut["alt_aa"]

        try:
            transcript = ensembl.transcript_by_id(transcript_id)
            wt_protein = transcript.protein_sequence
        except Exception:
            continue

        if not wt_protein:
            continue

        # Verify wildtype residue matches
        idx = pos - 1  # 0-based
        if idx >= len(wt_protein) or wt_protein[idx] != ref_aa:
            continue

        # Create mutant protein
        mut_protein = wt_protein[:idx] + alt_aa + wt_protein[idx + 1 :]

        for k in lengths:
            # Generate peptides spanning the mutation site
            start = max(0, idx - k + 1)
            end = min(len(mut_protein) - k + 1, idx + 1)
            for i in range(start, end):
                mut_pep = mut_protein[i : i + k]
                wt_pep = wt_protein[i : i + k]
                if not set(mut_pep).issubset(AA20):
                    continue
                if mut_pep == wt_pep:
                    continue  # should not happen, but safety check
                n_flank = mut_protein[max(0, i - flank_length) : i]
                c_flank = mut_protein[i + k : i + k + flank_length]
                rows.append(
                    {
                        "label": mut["label"],
                        "gene": mut["gene"],
                        "gene_id": mut["gene_id"],
                        "mutation": f"{ref_aa}{pos}{alt_aa}",
                        "ref_aa": ref_aa,
                        "alt_aa": alt_aa,
                        "protein_position": pos,
                        "peptide": mut_pep,
                        "wildtype_peptide": wt_pep,
                        "length": k,
                        "start": i + 1,
                        "end": i + k,
                        "n_flank": n_flank,
                        "c_flank": c_flank,
                        "mutation_offset": idx - i,
                    }
                )

    return pd.DataFrame(rows)


def mutant_iedb_overlap(
    mutations: list[dict] | None = None,
    lengths: tuple[int, ...] = (8, 9, 10, 11),
    ensembl_release: int = 112,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    mhc_class: str | None = "I",
) -> pd.DataFrame:
    """Find mutant-spanning peptides with public IEDB/CEDAR mass spec support.

    Parameters
    ----------
    mutations
        List of mutation dicts. If None, uses all hotspot mutations.
    lengths
        Peptide lengths (default 8-11).
    ensembl_release
        Ensembl release (default 112).
    iedb_path
        Path to IEDB MHC ligand export.
    cedar_path
        Path to CEDAR MHC ligand export.
    mhc_class
        MHC class filter (default ``"I"``).

    Returns
    -------
    pd.DataFrame
        Mutant peptide DataFrame with ``has_iedb_hit``, ``iedb_alleles``,
        ``iedb_hit_count`` columns.
    """
    from .iedb import scan_public_ms

    mdf = mutant_peptides(
        mutations=mutations,
        lengths=lengths,
        ensembl_release=ensembl_release,
    )
    if mdf.empty:
        mdf["has_iedb_hit"] = pd.Series(dtype=bool)
        mdf["iedb_alleles"] = pd.Series(dtype=str)
        mdf["iedb_hit_count"] = pd.Series(dtype=int)
        return mdf

    unique_peptides = set(mdf["peptide"].unique())
    hits = scan_public_ms(
        peptides=unique_peptides,
        iedb_path=iedb_path,
        cedar_path=cedar_path,
        mhc_class=mhc_class,
    )

    if hits.empty:
        mdf["has_iedb_hit"] = False
        mdf["iedb_alleles"] = ""
        mdf["iedb_hit_count"] = 0
        return mdf

    hit_summary = (
        hits.groupby("peptide")
        .agg(
            iedb_hit_count=("peptide", "size"),
            iedb_alleles=("mhc_restriction", lambda x: ";".join(sorted(set(x)))),
        )
        .reset_index()
    )

    mdf = mdf.merge(hit_summary, on="peptide", how="left")
    mdf["has_iedb_hit"] = mdf["iedb_hit_count"].notna() & (mdf["iedb_hit_count"] > 0)
    mdf["iedb_hit_count"] = mdf["iedb_hit_count"].fillna(0).astype(int)
    mdf["iedb_alleles"] = mdf["iedb_alleles"].fillna("")
    return mdf
