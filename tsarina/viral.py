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


def _proteome_kmers_from_sequences(
    sequences: dict[str, str] | list[str] | tuple[str, ...],
    lengths: tuple[int, ...],
) -> frozenset[str]:
    values = sequences.values() if isinstance(sequences, dict) else sequences
    kmers: set[str] = set()
    for seq in values:
        if not seq:
            continue
        seq_str = str(seq)
        n = len(seq_str)
        for k in lengths:
            for i in range(n - k + 1):
                kmers.add(seq_str[i : i + k])
    return frozenset(kmers)


def host_exclusive_viral_peptides(
    virus: str | None = None,
    fasta_path: str | Path | None = None,
    proteins: dict[str, str] | None = None,
    lengths: tuple[int, ...] = DEFAULT_PEPTIDE_LENGTHS,
    ensembl_release: int = 112,
    *,
    host_species: str = "human",
    host_proteins: dict[str, str] | list[str] | None = None,
) -> pd.DataFrame:
    """Return viral peptides NOT found in any protein of *host_species* (#217).

    Generates all k-mer peptides from the viral proteome and removes any that
    appear as a substring of any protein-coding sequence in the host proteome
    (from *host_proteins* if supplied, else ``hitlist.proteome.proteome_kmer_set``
    scoped to *host_species* and *ensembl_release*).
    """
    from .species import resolve_species_context

    ctx = resolve_species_context(
        host_species,
        antigen_source_species=virus or "viral",
        ensembl_release=ensembl_release,
    )
    vdf = viral_peptides(virus=virus, fasta_path=fasta_path, proteins=proteins, lengths=lengths)
    if vdf.empty:
        vdf.attrs["species_context"] = ctx.to_provenance()
        return vdf

    if host_proteins is not None:
        host_kmers = _proteome_kmers_from_sequences(host_proteins, tuple(lengths))
    else:
        from hitlist.proteome import proteome_kmer_set

        kmer_kwargs = {"release": ensembl_release, "lengths": tuple(lengths), "gene_ids": None}
        if not ctx.is_human_host:
            kmer_kwargs["species"] = ctx.pyensembl_species
        host_kmers = proteome_kmer_set(**kmer_kwargs)
    mask = ~vdf["peptide"].isin(host_kmers)
    out = vdf[mask].reset_index(drop=True)
    out.attrs["species_context"] = ctx.to_provenance()
    return out


def human_exclusive_viral_peptides(
    virus: str | None = None,
    fasta_path: str | Path | None = None,
    proteins: dict[str, str] | None = None,
    lengths: tuple[int, ...] = DEFAULT_PEPTIDE_LENGTHS,
    ensembl_release: int = 112,
) -> pd.DataFrame:
    """Return viral peptides NOT found in any human protein.

    Thin wrapper around :func:`host_exclusive_viral_peptides` with
    ``host_species="human"``.
    """
    return host_exclusive_viral_peptides(
        virus=virus,
        fasta_path=fasta_path,
        proteins=proteins,
        lengths=lengths,
        ensembl_release=ensembl_release,
        host_species="human",
    )


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
    :func:`hitlist.proteome.proteome_kmer_set` and cached both
    in-process and on disk (hitlist 1.30.49+).  The on-disk
    ``ProteomeIndex`` cache amortizes the Ensembl build across fresh
    processes — relevant for per-patient cluster jobs that re-invoke
    this entry point in a shell loop.
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
    *,
    host_species: str = "human",
    mhc_species: str | None = None,
    requested_alleles: list[str] | tuple[str, ...] | None = None,
    evidence_scope: str = "all",
    observations: pd.DataFrame | None = None,
    host_proteins: dict[str, str] | list[str] | None = None,
) -> pd.DataFrame:
    """Find viral peptides with public IEDB/CEDAR mass spec support.

    Supports species-scoped host self-exclusion and MHC restriction filtering
    (#217) as well as multi-tier MS attribution (#216) when *requested_alleles*
    or *observations* are provided.
    """
    from .ms_evidence import aggregate_ms_hits_for_iedb_columns, load_public_ms_hits
    from .species import resolve_species_context

    ctx = resolve_species_context(
        host_species,
        antigen_source_species=virus or "viral",
        mhc_species=mhc_species,
        ensembl_release=ensembl_release,
    )

    if human_exclusive_only:
        vdf = host_exclusive_viral_peptides(
            virus=virus,
            fasta_path=fasta_path,
            proteins=proteins,
            lengths=lengths,
            ensembl_release=ensembl_release,
            host_species=host_species,
            host_proteins=host_proteins,
        )
    else:
        vdf = viral_peptides(virus=virus, fasta_path=fasta_path, proteins=proteins, lengths=lengths)

    if vdf.empty:
        vdf["has_iedb_hit"] = pd.Series(dtype=bool)
        vdf["iedb_alleles"] = pd.Series(dtype=str)
        vdf["iedb_hit_count"] = pd.Series(dtype=int)
        vdf.attrs["species_context"] = ctx.to_provenance()
        return vdf

    unique_peptides = set(vdf["peptide"].unique())
    if observations is not None:
        hits = observations.copy()
    else:
        hits = load_public_ms_hits(
            peptides=unique_peptides,
            iedb_path=iedb_path,
            cedar_path=cedar_path,
            mhc_class=mhc_class,
            mhc_species=ctx.mhc_species,
            classify_source=False,
            drop_binding_assays=True,
        )

    if requested_alleles is not None or evidence_scope != "all" or observations is not None:
        from .ms_attribution import attribute_ms_evidence, summarize_ms_attribution

        attributed = attribute_ms_evidence(
            vdf[["peptide"]].drop_duplicates(),
            hits,
            requested_alleles=requested_alleles,
            species_context=ctx,
            evidence_scope=evidence_scope,
            evidence_loaded=True,
        )
        summary = summarize_ms_attribution(attributed, group_cols=("peptide",))
        vdf = vdf.merge(summary, on="peptide", how="left")
        vdf["iedb_hit_count"] = (
            pd.to_numeric(vdf["ms_hit_count"], errors="coerce").fillna(0).astype(int)
        )
        vdf["iedb_alleles"] = vdf["ms_alleles"].fillna("")
        vdf["has_iedb_hit"] = vdf["iedb_hit_count"] > 0
        if evidence_scope != "all":
            vdf = vdf[vdf["has_iedb_hit"]].reset_index(drop=True)
        vdf.attrs["species_context"] = ctx.to_provenance()
        return vdf

    if hits.empty:
        vdf["has_iedb_hit"] = False
        vdf["iedb_alleles"] = ""
        vdf["iedb_hit_count"] = 0
        vdf.attrs["species_context"] = ctx.to_provenance()
        return vdf

    # Aggregate IEDB hits per peptide
    hit_summary = aggregate_ms_hits_for_iedb_columns(hits)

    vdf = vdf.merge(hit_summary, on="peptide", how="left")
    vdf["has_iedb_hit"] = vdf["iedb_hit_count"].notna() & (vdf["iedb_hit_count"] > 0)
    vdf["iedb_hit_count"] = vdf["iedb_hit_count"].fillna(0).astype(int)
    vdf["iedb_alleles"] = vdf["iedb_alleles"].fillna("")
    vdf.attrs["species_context"] = ctx.to_provenance()
    return vdf


def _require_orfdata():
    try:
        import orfdata
    except ImportError as exc:
        raise ImportError(
            "orfdata is required for GDC IARC viral subsets and RefSeq viral ORFs; "
            "install with `pip install orfdata`."
        ) from exc
    return orfdata


def orfdata_viral_subsets() -> tuple[dict[str, object], ...]:
    """Return IARC cancer-viral subsets exposed by optional ``orfdata``.

    Each entry contains ``name``, ``title``, ``description``,
    ``abbreviations`` (tuple of virus abbreviations), and ``virus_count``.
    Works completely offline once ``orfdata`` is installed; raises
    :class:`ImportError` with installation guidance when ``orfdata`` is absent.
    """
    orfdata = _require_orfdata()
    subsets: list[dict[str, object]] = []
    for entry in orfdata.viral_subsets():
        abbrevs = tuple(str(a).strip() for a in entry.get("abbreviations", ()))
        subsets.append(
            {
                "name": str(entry.get("name", "")),
                "title": str(entry.get("title", "")),
                "description": str(entry.get("description", "")),
                "abbreviations": abbrevs,
                "virus_count": len(abbrevs),
            }
        )
    return tuple(subsets)


def orfdata_viral_peptides(
    *,
    reference_id: str | None = None,
    subsets: str | list[str] | tuple[str, ...] | None = None,
    lengths: tuple[int, ...] = DEFAULT_PEPTIDE_LENGTHS,
    installed: bool = False,
    download_if_missing: bool = False,
    host_exclusive: bool = False,
    ensembl_release: int = 112,
    host_species: str = "human",
    host_proteins: dict[str, str] | list[str] | None = None,
) -> pd.DataFrame:
    """Generate viral peptides from ``orfdata`` RefSeq viral ORFs or GDC viral subsets.

    Provide exactly one of ``reference_id`` (uses :func:`orfdata.viral_orfs`,
    including bundled spliced CDS annotations such as HPV16 ``NC_001526.4``
    ``E1^E4`` and ``E6*``) or ``subsets`` (uses
    :func:`orfdata.select_viral_references` over the IARC GDC viral panel, e.g.
    ``"cancer-causing"``, ``"cancer-associated"``, or
    ``"cancer-integrated-driver"``).
    """
    if (reference_id is None) == (subsets is None):
        raise ValueError("Provide exactly one of reference_id or subsets.")

    orfdata = _require_orfdata()

    if reference_id is not None:
        orfs = orfdata.viral_orfs(reference_id=reference_id, installed=installed)
        proteins: dict[str, str] = {}
        gene_by_prot: dict[str, str] = {}
        for orf in orfs:
            seq = str(getattr(orf, "protein_sequence", "") or "").strip().upper()
            if not seq:
                continue
            prot_id = str(getattr(orf, "annotation_id", None) or getattr(orf, "orf_id", ""))
            proteins[prot_id] = seq
            gene_by_prot[prot_id] = str(getattr(orf, "gene_name", "") or "")
        if host_exclusive:
            df = host_exclusive_viral_peptides(
                proteins=proteins,
                lengths=lengths,
                ensembl_release=ensembl_release,
                host_species=host_species,
                host_proteins=host_proteins,
            )
        else:
            df = viral_peptides(proteins=proteins, lengths=lengths)
        if not df.empty:
            df["virus"] = str(reference_id)
            df["gene_name"] = df["protein_id"].map(gene_by_prot).fillna("")
        else:
            df["gene_name"] = pd.Series(dtype=str)
        return df

    subset_list = [subsets] if isinstance(subsets, str) else list(subsets)  # type: ignore[arg-type]
    refs = orfdata.select_viral_references(
        subsets=subset_list,
        download_if_missing=download_if_missing,
    )
    if host_exclusive:
        df = host_exclusive_viral_peptides(
            fasta_path=refs.protein_fasta,
            lengths=lengths,
            ensembl_release=ensembl_release,
            host_species=host_species,
            host_proteins=host_proteins,
        )
    else:
        df = viral_peptides(fasta_path=refs.protein_fasta, lengths=lengths)
    if not df.empty:
        df["virus"] = ",".join(subset_list)
    return df
