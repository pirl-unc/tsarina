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

"""Generate peptides from CTA protein sequences and check CTA exclusivity.

Requires the ``pyensembl`` package (install with ``pip install hitlist[peptides]``).

Typical usage::

    from hitlist.peptides import cta_peptides, cta_exclusive_peptides

    # All 8-11mer peptides from expressed CTAs
    df = cta_peptides()

    # Only peptides found exclusively in CTA proteins (not in any non-CTA protein)
    exclusive = cta_exclusive_peptides()
"""

from __future__ import annotations

from typing import Union

import pandas as pd

AA20 = set("ACDEFGHIKLMNPQRSTVWY")

DEFAULT_PEPTIDE_LENGTHS = (8, 9, 10, 11)


def cta_peptides(
    ensembl_release: int = 112,
    lengths: tuple[int, ...] = DEFAULT_PEPTIDE_LENGTHS,
    flank_length: int = 15,
) -> pd.DataFrame:
    """Generate all peptides from expressed CTA canonical protein sequences.

    Parameters
    ----------
    ensembl_release
        Ensembl release for protein sequences (default 112).
    lengths
        Peptide lengths to generate (default 8, 9, 10, 11).
    flank_length
        Number of flanking residues to include for processing prediction
        (default 15).

    Returns
    -------
    pd.DataFrame
        Columns: ``gene_name``, ``gene_id``, ``transcript_id``,
        ``peptide``, ``length``, ``start`` (1-based), ``end``,
        ``n_flank``, ``c_flank``.
    """
    from pyensembl import EnsemblRelease

    from .gene_sets import CTA_gene_ids

    ensembl = EnsemblRelease(ensembl_release)
    cta_ids = CTA_gene_ids()

    rows: list[dict] = []
    for gene_id in sorted(cta_ids):
        try:
            gene = ensembl.gene_by_id(gene_id)
        except ValueError:
            continue

        # Pick the canonical (longest) protein-coding transcript
        best_transcript = None
        best_length = 0
        for t in gene.transcripts:
            if t.biotype != "protein_coding":
                continue
            try:
                seq = t.protein_sequence
            except Exception:
                continue
            if seq and len(seq) > best_length:
                best_transcript = t
                best_length = len(seq)

        if best_transcript is None:
            continue

        protein = best_transcript.protein_sequence
        if not protein:
            continue

        for k in lengths:
            for i in range(len(protein) - k + 1):
                pep = protein[i : i + k]
                if not set(pep).issubset(AA20):
                    continue
                n_flank = protein[max(0, i - flank_length) : i]
                c_flank = protein[i + k : i + k + flank_length]
                rows.append(
                    {
                        "gene_name": gene.name,
                        "gene_id": gene_id,
                        "transcript_id": best_transcript.id,
                        "peptide": pep,
                        "length": k,
                        "start": i + 1,
                        "end": i + k,
                        "n_flank": n_flank,
                        "c_flank": c_flank,
                    }
                )

    return pd.DataFrame(rows)


def cta_exclusive_peptides(
    ensembl_release: int = 112,
    lengths: tuple[int, ...] = DEFAULT_PEPTIDE_LENGTHS,
) -> pd.DataFrame:
    """Return CTA peptides that do NOT appear in any non-CTA protein.

    This filters the output of :func:`cta_peptides` to only include
    peptides whose sequence is exclusive to CTA proteins -- i.e., the
    peptide is not found as a substring of any non-CTA protein-coding
    gene's canonical protein sequence.

    Parameters
    ----------
    ensembl_release
        Ensembl release (default 112).
    lengths
        Peptide lengths to generate (default 8, 9, 10, 11).

    Returns
    -------
    pd.DataFrame
        Same columns as :func:`cta_peptides`, filtered to exclusive peptides.
    """
    from pyensembl import EnsemblRelease

    from .partition import CTA_partition_gene_ids

    ensembl = EnsemblRelease(ensembl_release)
    partition = CTA_partition_gene_ids(ensembl_release)

    # Build set of all peptide k-mers from non-CTA proteins
    noncta_peptides: set[str] = set()
    for gene_id in partition.non_cta:
        try:
            gene = ensembl.gene_by_id(gene_id)
        except ValueError:
            continue
        for t in gene.transcripts:
            if t.biotype != "protein_coding":
                continue
            try:
                seq = t.protein_sequence
            except Exception:
                continue
            if not seq:
                continue
            for k in lengths:
                for i in range(len(seq) - k + 1):
                    pep = seq[i : i + k]
                    if set(pep).issubset(AA20):
                        noncta_peptides.add(pep)

    # Generate CTA peptides and filter
    cta_df = cta_peptides(ensembl_release=ensembl_release, lengths=lengths)
    mask = ~cta_df["peptide"].isin(noncta_peptides)
    return cta_df[mask].reset_index(drop=True)


def build_pmhc_table(
    peptide_df: pd.DataFrame,
    alleles: Union[list[str], None] = None,
) -> pd.DataFrame:
    """Combine a peptide DataFrame with an allele panel to produce a pMHC table.

    Parameters
    ----------
    peptide_df
        DataFrame with at least a ``peptide`` column (and typically
        ``gene_name``, ``gene_id``, etc.).
    alleles
        List of HLA allele names.  If None, uses the IEDB-27 panel.

    Returns
    -------
    pd.DataFrame
        Cross-product of peptides and alleles, with an ``allele`` column added.
    """
    if alleles is None:
        from .alleles import IEDB27_AB

        alleles = IEDB27_AB

    allele_df = pd.DataFrame({"allele": alleles, "_key": 1})
    peptide_cross = peptide_df.assign(_key=1)
    result = peptide_cross.merge(allele_df, on="_key").drop(columns="_key")
    return result
