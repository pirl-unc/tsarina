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

"""Oncogenic virus proteome peptide generation and IEDB overlap analysis.

Generates peptide libraries from oncogenic virus proteins, filters to
peptides exclusive to the viral proteome (not found in the human
proteome), and cross-references with public IEDB/CEDAR mass spec data.

Supported viruses::

    HPV-16, HPV-18       Cervical and other cancers
    EBV/HHV-4            Burkitt lymphoma, NPC, Hodgkin lymphoma
    HTLV-1               Adult T-cell leukemia/lymphoma
    HBV                  Hepatocellular carcinoma
    HCV                  Hepatocellular carcinoma
    KSHV/HHV-8           Kaposi sarcoma
    MCPyV                Merkel cell carcinoma
    HIV-1                Kaposi sarcoma, lymphoma (indirect)

Typical usage::

    from tsarina.viral import (
        ONCOGENIC_VIRUSES,
        viral_peptides,
        human_exclusive_viral_peptides,
        viral_iedb_overlap,
    )

    # Generate all 8-11mer peptides from HPV-16 proteins
    peps = viral_peptides("hpv16")

    # Filter to peptides not found in any human protein
    exclusive = human_exclusive_viral_peptides("hpv16")

    # Find which viral peptides have IEDB mass spec support
    overlap = viral_iedb_overlap("hpv16", iedb_path="mhc_ligand_full.csv")
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .peptides import AA20, _cta_gene_ids, _non_cta_gene_ids

# ── Virus definitions ───────────────────────────────────────────────────────

ONCOGENIC_VIRUSES: dict[str, dict] = {
    "hpv16": {
        "name": "HPV-16",
        "full_name": "Human papillomavirus type 16",
        "cancers": ["cervical", "oropharyngeal", "anal", "penile"],
        "uniprot_proteome": "UP000006729",
        "taxonomy_id": 333760,
        "key_oncoproteins": ["E6", "E7"],
    },
    "hpv18": {
        "name": "HPV-18",
        "full_name": "Human papillomavirus type 18",
        "cancers": ["cervical", "oropharyngeal"],
        "uniprot_proteome": "UP000006728",
        "taxonomy_id": 333761,
        "key_oncoproteins": ["E6", "E7"],
    },
    "ebv": {
        "name": "EBV",
        "full_name": "Epstein-Barr virus (HHV-4)",
        "cancers": ["Burkitt lymphoma", "nasopharyngeal carcinoma", "Hodgkin lymphoma", "PTLD"],
        "uniprot_proteome": "UP000153037",
        "taxonomy_id": 10376,
        "key_oncoproteins": ["LMP1", "LMP2A", "EBNA1", "EBNA2"],
    },
    "htlv1": {
        "name": "HTLV-1",
        "full_name": "Human T-lymphotropic virus 1",
        "cancers": ["adult T-cell leukemia/lymphoma"],
        "uniprot_proteome": "UP000002063",
        "taxonomy_id": 11908,
        "key_oncoproteins": ["Tax", "HBZ"],
    },
    "hbv": {
        "name": "HBV",
        "full_name": "Hepatitis B virus",
        "cancers": ["hepatocellular carcinoma"],
        "uniprot_proteome": "UP000126453",
        "taxonomy_id": 10407,
        "key_oncoproteins": ["HBx"],
    },
    "hcv": {
        "name": "HCV",
        "full_name": "Hepatitis C virus",
        "cancers": ["hepatocellular carcinoma", "B-cell lymphoma"],
        "uniprot_proteome": "UP000000518",
        "taxonomy_id": 11103,
        "key_oncoproteins": ["Core", "NS3", "NS5A"],
    },
    "kshv": {
        "name": "KSHV",
        "full_name": "Kaposi sarcoma-associated herpesvirus (HHV-8)",
        "cancers": ["Kaposi sarcoma", "primary effusion lymphoma", "multicentric Castleman"],
        "uniprot_proteome": "UP000009113",
        "taxonomy_id": 37296,
        "key_oncoproteins": ["vFLIP", "vCyclin", "LANA"],
    },
    "mcpyv": {
        "name": "MCPyV",
        "full_name": "Merkel cell polyomavirus",
        "cancers": ["Merkel cell carcinoma"],
        "uniprot_proteome": "UP000116695",
        "taxonomy_id": 493803,
        "key_oncoproteins": ["LT", "sT"],
    },
    "hiv1": {
        "name": "HIV-1",
        "full_name": "Human immunodeficiency virus 1",
        "cancers": ["Kaposi sarcoma", "non-Hodgkin lymphoma", "cervical cancer (indirect)"],
        "uniprot_proteome": "UP000002241",
        "taxonomy_id": 11676,
        "key_oncoproteins": ["Tat", "Nef"],
    },
}


def read_fasta(path: str | Path) -> dict[str, str]:
    """Read a FASTA file into a dict of {header: sequence}.

    Parameters
    ----------
    path
        Path to a FASTA file.

    Returns
    -------
    dict[str, str]
        Mapping from sequence header (first line after ``>``, up to first
        whitespace) to the full amino acid sequence.
    """
    proteins: dict[str, str] = {}
    current_id = ""
    current_seq: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    proteins[current_id] = "".join(current_seq)
                # Use full header line (minus >) as ID, split on space for short ID
                header = line[1:].strip()
                current_id = header.split()[0] if header else header
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        proteins[current_id] = "".join(current_seq)
    return proteins


def _get_proteins(virus: str) -> dict[str, str]:
    """Load viral proteins, fetching if needed."""
    from .downloads import fetch, get_path

    try:
        path = get_path(virus)
    except (KeyError, FileNotFoundError):
        path = fetch(virus)
    return read_fasta(path)


DEFAULT_PEPTIDE_LENGTHS = (8, 9, 10, 11)


def viral_peptides(
    virus: str | None = None,
    fasta_path: str | Path | None = None,
    proteins: dict[str, str] | None = None,
    lengths: tuple[int, ...] = DEFAULT_PEPTIDE_LENGTHS,
) -> pd.DataFrame:
    """Generate all peptides from a viral proteome.

    Provide exactly one of ``virus``, ``fasta_path``, or ``proteins``.

    Parameters
    ----------
    virus
        Virus key from :data:`ONCOGENIC_VIRUSES` (e.g. ``"hpv16"``).
        Auto-fetches the proteome if not already downloaded.
    fasta_path
        Path to a FASTA file with viral protein sequences.
    proteins
        Dict of ``{protein_id: sequence}``.
    lengths
        Peptide lengths to generate (default 8, 9, 10, 11).

    Returns
    -------
    pd.DataFrame
        Columns: ``virus``, ``protein_id``, ``peptide``, ``length``,
        ``start`` (1-based), ``end``.
    """
    if sum(x is not None for x in (virus, fasta_path, proteins)) != 1:
        raise ValueError("Provide exactly one of virus, fasta_path, or proteins")

    virus_label = virus or "custom"
    if virus is not None:
        if virus not in ONCOGENIC_VIRUSES:
            raise ValueError(f"Unknown virus '{virus}'. Available: {sorted(ONCOGENIC_VIRUSES)}")
        prots = _get_proteins(virus)
        virus_label = ONCOGENIC_VIRUSES[virus]["name"]
    elif fasta_path is not None:
        prots = read_fasta(fasta_path)
    else:
        prots = proteins  # type: ignore[assignment]

    rows: list[dict] = []
    for prot_id, seq in sorted(prots.items()):
        for k in lengths:
            for i in range(len(seq) - k + 1):
                pep = seq[i : i + k]
                if not set(pep).issubset(AA20):
                    continue
                rows.append(
                    {
                        "virus": virus_label,
                        "protein_id": prot_id,
                        "peptide": pep,
                        "length": k,
                        "start": i + 1,
                        "end": i + k,
                    }
                )
    return pd.DataFrame(rows)


def human_exclusive_viral_peptides(
    virus: str | None = None,
    fasta_path: str | Path | None = None,
    proteins: dict[str, str] | None = None,
    lengths: tuple[int, ...] = DEFAULT_PEPTIDE_LENGTHS,
    ensembl_release: int = 112,
) -> pd.DataFrame:
    """Return viral peptides NOT found in any human protein.

    Generates all k-mer peptides from the viral proteome and removes
    any that appear as a substring of any human protein-coding gene's
    protein sequence (from Ensembl).

    Requires ``pyensembl`` (install with ``pip install tsarina[peptides]``).

    Parameters
    ----------
    virus, fasta_path, proteins
        Viral proteome source (see :func:`viral_peptides`).
    lengths
        Peptide lengths (default 8-11).
    ensembl_release
        Ensembl release for human proteome (default 112).

    Returns
    -------
    pd.DataFrame
        Same columns as :func:`viral_peptides`, filtered to human-exclusive.

    Notes
    -----
    The human k-mer set is computed once per
    ``(ensembl_release, lengths)`` via
    :func:`hitlist.proteome.proteome_kmer_set` and cached process-wide.
    First call pays the Ensembl walk once; subsequent calls with the
    same inputs return immediately.  The cache is shared across any
    sibling pirl-unc package that calls the same primitive.
    """
    from hitlist.proteome import proteome_kmer_set

    vdf = viral_peptides(virus=virus, fasta_path=fasta_path, proteins=proteins, lengths=lengths)
    if vdf.empty:
        return vdf
    human_kmers = proteome_kmer_set(release=ensembl_release, lengths=tuple(lengths), gene_ids=None)
    mask = ~vdf["peptide"].isin(human_kmers)
    return vdf[mask].reset_index(drop=True)


def cancer_specific_viral_peptides(
    virus: str | None = None,
    fasta_path: str | Path | None = None,
    proteins: dict[str, str] | None = None,
    lengths: tuple[int, ...] = DEFAULT_PEPTIDE_LENGTHS,
    ensembl_release: int = 112,
) -> pd.DataFrame:
    """Return viral peptides NOT found in non-CTA human proteins.

    Unlike :func:`human_exclusive_viral_peptides` which excludes peptides
    found in *any* human protein, this function only excludes peptides found
    in **non-CTA** human proteins.  Peptides shared with CTA proteins are
    kept, since CTAs are themselves cancer-specific targets.

    The returned DataFrame includes an ``in_cta_protein`` column indicating
    whether each peptide also occurs in a CTA protein sequence.

    Requires ``pyensembl`` (install with ``pip install tsarina[peptides]``).

    Parameters
    ----------
    virus, fasta_path, proteins
        Viral proteome source (see :func:`viral_peptides`).
    lengths
        Peptide lengths (default 8-11).
    ensembl_release
        Ensembl release for human proteome (default 112).

    Returns
    -------
    pd.DataFrame
        Same columns as :func:`viral_peptides`, plus ``in_cta_protein``.

    Notes
    -----
    Both the CTA and non-CTA k-mer sets are computed via
    :func:`hitlist.proteome.proteome_kmer_set` and cached process-wide,
    shared with :func:`cta_exclusive_peptides` /
    :func:`human_exclusive_viral_peptides` where their gene filters
    overlap.  First call pays the Ensembl walk; subsequent calls with
    the same ``(release, lengths)`` return in sub-millisecond time.
    """
    from hitlist.proteome import proteome_kmer_set

    vdf = viral_peptides(virus=virus, fasta_path=fasta_path, proteins=proteins, lengths=lengths)
    if vdf.empty:
        vdf["in_cta_protein"] = pd.Series(dtype=bool)
        return vdf

    lengths_t = tuple(lengths)
    noncta_kmers = proteome_kmer_set(
        release=ensembl_release, lengths=lengths_t, gene_ids=_non_cta_gene_ids(ensembl_release)
    )
    cta_kmers = proteome_kmer_set(
        release=ensembl_release, lengths=lengths_t, gene_ids=_cta_gene_ids(ensembl_release)
    )

    # Keep peptides NOT in non-CTA proteins; annotate which also occur in CTA proteins.
    mask = ~vdf["peptide"].isin(noncta_kmers)
    result = vdf[mask].copy().reset_index(drop=True)
    result["in_cta_protein"] = result["peptide"].isin(cta_kmers)
    return result


def viral_iedb_overlap(
    virus: str | None = None,
    fasta_path: str | Path | None = None,
    proteins: dict[str, str] | None = None,
    lengths: tuple[int, ...] = DEFAULT_PEPTIDE_LENGTHS,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    mhc_class: str | None = "I",
    human_exclusive_only: bool = False,
    ensembl_release: int = 112,
) -> pd.DataFrame:
    """Find viral peptides with public IEDB/CEDAR mass spec support.

    Parameters
    ----------
    virus, fasta_path, proteins
        Viral proteome source (see :func:`viral_peptides`).
    lengths
        Peptide lengths (default 8-11).
    iedb_path
        Path to IEDB MHC ligand export.
    cedar_path
        Path to CEDAR MHC ligand export.
    mhc_class
        MHC class filter (default ``"I"``).
    human_exclusive_only
        If True, only include peptides not found in the human proteome.
        Requires ``pyensembl``.
    ensembl_release
        Ensembl release for human proteome filtering.

    Returns
    -------
    pd.DataFrame
        Viral peptide DataFrame merged with IEDB hit information.
        Columns include all viral_peptides columns plus ``has_iedb_hit``,
        ``iedb_alleles`` (semicolon-separated MHC restrictions),
        ``iedb_hit_count``.
    """
    from .iedb import scan_public_ms

    if human_exclusive_only:
        vdf = human_exclusive_viral_peptides(
            virus=virus,
            fasta_path=fasta_path,
            proteins=proteins,
            lengths=lengths,
            ensembl_release=ensembl_release,
        )
    else:
        vdf = viral_peptides(virus=virus, fasta_path=fasta_path, proteins=proteins, lengths=lengths)

    if vdf.empty:
        vdf["has_iedb_hit"] = pd.Series(dtype=bool)
        vdf["iedb_alleles"] = pd.Series(dtype=str)
        vdf["iedb_hit_count"] = pd.Series(dtype=int)
        return vdf

    unique_peptides = set(vdf["peptide"].unique())
    hits = scan_public_ms(
        peptides=unique_peptides,
        iedb_path=iedb_path,
        cedar_path=cedar_path,
        mhc_class=mhc_class,
    )

    if hits.empty:
        vdf["has_iedb_hit"] = False
        vdf["iedb_alleles"] = ""
        vdf["iedb_hit_count"] = 0
        return vdf

    # Aggregate IEDB hits per peptide
    hit_summary = (
        hits.groupby("peptide")
        .agg(
            iedb_hit_count=("peptide", "size"),
            iedb_alleles=("mhc_restriction", lambda x: ";".join(sorted(set(x)))),
        )
        .reset_index()
    )

    vdf = vdf.merge(hit_summary, on="peptide", how="left")
    vdf["has_iedb_hit"] = vdf["iedb_hit_count"].notna() & (vdf["iedb_hit_count"] > 0)
    vdf["iedb_hit_count"] = vdf["iedb_hit_count"].fillna(0).astype(int)
    vdf["iedb_alleles"] = vdf["iedb_alleles"].fillna("")
    return vdf
