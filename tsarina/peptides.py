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

Requires the ``pyensembl`` package (install with ``pip install tsarina[peptides]``).

Typical usage::

    from tsarina.peptides import cta_peptides, cta_exclusive_peptides

    # All 8-11mer peptides from expressed CTAs
    df = cta_peptides()

    # Only peptides found exclusively in CTA proteins (not in any non-CTA protein)
    exclusive = cta_exclusive_peptides()
"""

from __future__ import annotations

import sys
from collections.abc import Callable, Iterable, Iterator
from functools import lru_cache
from typing import TextIO, Union

import pandas as pd

AA20 = set("ACDEFGHIKLMNPQRSTVWY")

DEFAULT_PEPTIDE_LENGTHS = (8, 9, 10, 11)


def _report_progress(on_progress: Callable[[str], None] | None, message: str) -> None:
    if on_progress is not None:
        on_progress(message)


def _progress_iter(
    values: Iterable,
    *,
    desc: str,
    progress_bar: bool,
    progress_file: TextIO | None,
) -> Iterator:
    if not progress_bar:
        yield from values
        return

    from tqdm import tqdm

    yield from tqdm(values, desc=desc, unit="gene", file=progress_file or sys.stderr, leave=False)


def _split_semicolon_values(value: object) -> list[str]:
    if not isinstance(value, str):
        return []
    return [part.strip() for part in value.split(";") if part.strip()]


def _cta_gene_ids_for_names(gene_names: Iterable[str] | None) -> list[str]:
    from .gene_sets import CTA_gene_ids

    all_cta_ids = CTA_gene_ids()
    if gene_names is None:
        return sorted(all_cta_ids)

    wanted = {str(name).strip() for name in gene_names if str(name).strip()}
    if not wanted:
        return []

    from .loader import cta_dataframe

    df = cta_dataframe()
    if "Symbol" not in df.columns or "Ensembl_Gene_ID" not in df.columns:
        return []

    selected: set[str] = set()
    for symbol_cell, id_cell in df[["Symbol", "Ensembl_Gene_ID"]].itertuples(
        index=False,
        name=None,
    ):
        if wanted.isdisjoint(_split_semicolon_values(symbol_cell)):
            continue
        selected.update(
            gene_id for gene_id in _split_semicolon_values(id_cell) if gene_id in all_cta_ids
        )
    return sorted(selected)


@lru_cache(maxsize=4)
def _non_cta_gene_ids(ensembl_release: int) -> frozenset[str]:
    """Cached frozenset of non-CTA Ensembl gene IDs for a release.

    Pre-building the frozenset means downstream
    :func:`hitlist.proteome.proteome_kmer_set` calls reuse a single
    hashable instance rather than rebuilding from the mutable
    ``CTA_partition_gene_ids(...).non_cta`` set on every call.  The
    CTA / non-CTA partition is tsarina-specific so it stays here; the
    proteome walk itself is delegated to hitlist.
    """
    from .partition import CTA_partition_gene_ids

    return frozenset(CTA_partition_gene_ids(ensembl_release).non_cta)


@lru_cache(maxsize=4)
def _cta_gene_ids(ensembl_release: int) -> frozenset[str]:
    """Cached frozenset of CTA Ensembl gene IDs for a release."""
    from .partition import CTA_partition_gene_ids

    return frozenset(CTA_partition_gene_ids(ensembl_release).cta)


def cta_peptides(
    ensembl_release: int = 112,
    lengths: tuple[int, ...] = DEFAULT_PEPTIDE_LENGTHS,
    flank_length: int = 15,
    gene_names: Iterable[str] | None = None,
    on_progress: Callable[[str], None] | None = None,
    progress_bar: bool = False,
    progress_file: TextIO | None = None,
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
    gene_names
        Optional CTA gene symbols to enumerate. Defaults to all expressed CTAs.
    on_progress
        Optional progress callback.
    progress_bar
        If True, show a tqdm bar while walking CTA genes.
    progress_file
        File handle for tqdm output. Defaults to ``sys.stderr``.

    Returns
    -------
    pd.DataFrame
        Columns: ``gene_name``, ``gene_id``, ``transcript_id``,
        ``peptide``, ``length``, ``start`` (1-based), ``end``,
        ``n_flank``, ``c_flank``.
    """
    from pyensembl import EnsemblRelease

    ensembl = EnsemblRelease(ensembl_release)
    cta_ids = _cta_gene_ids_for_names(gene_names)
    _report_progress(
        on_progress,
        f"Generating CTA peptide candidates from {len(cta_ids)} CTA genes "
        f"(Ensembl {ensembl_release})...",
    )

    rows: list[dict] = []
    for gene_id in _progress_iter(
        cta_ids,
        desc="CTA genes",
        progress_bar=progress_bar,
        progress_file=progress_file,
    ):
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
            except (ValueError, KeyError, TypeError):
                # pyensembl raises these for transcripts lacking a usable
                # translation; let unexpected errors surface rather than
                # silently dropping the transcript.
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

    out = pd.DataFrame(rows)
    _report_progress(
        on_progress,
        f"Generated {len(out)} CTA peptide rows "
        f"({out['peptide'].nunique() if not out.empty else 0} unique peptides).",
    )
    return out


def _non_cta_overlapping_peptides(
    candidate_peptides: set[str],
    *,
    ensembl_release: int,
    lengths: tuple[int, ...],
    on_progress: Callable[[str], None] | None = None,
    progress_bar: bool = False,
    progress_file: TextIO | None = None,
) -> set[str]:
    """Find candidate peptides that appear in non-CTA protein-coding transcripts.

    This intentionally streams non-CTA proteins and checks membership in the
    CTA candidate set. It avoids building hitlist's full human k-mer posting
    index when tsarina only needs to subtract a much smaller set of CTA
    candidate peptides.
    """
    if not candidate_peptides:
        return set()

    from pyensembl import EnsemblRelease

    non_cta_gene_ids = sorted(_non_cta_gene_ids(ensembl_release))
    _report_progress(
        on_progress,
        f"Scanning {len(non_cta_gene_ids)} non-CTA genes for overlapping "
        f"{','.join(str(length) for length in lengths)}-mers...",
    )
    ensembl = EnsemblRelease(ensembl_release)
    seen: set[str] = set()
    for gene_id in _progress_iter(
        non_cta_gene_ids,
        desc="Non-CTA genes",
        progress_bar=progress_bar,
        progress_file=progress_file,
    ):
        try:
            gene = ensembl.gene_by_id(gene_id)
        except ValueError:
            continue
        for transcript in gene.transcripts:
            if getattr(transcript, "biotype", "") != "protein_coding":
                continue
            try:
                protein = transcript.protein_sequence
            except (ValueError, KeyError, TypeError):
                # pyensembl raises these for transcripts lacking a usable
                # translation. Narrowed deliberately: silently swallowing every
                # error here could skip a non-CTA protein and let a peptide pass
                # as CTA-exclusive when it is not.
                continue
            if not protein:
                continue
            protein_len = len(protein)
            for length in lengths:
                end = protein_len - length + 1
                for i in range(end):
                    peptide = protein[i : i + length]
                    if peptide in candidate_peptides:
                        seen.add(peptide)
        if len(seen) == len(candidate_peptides):
            break

    _report_progress(
        on_progress,
        f"Found {len(seen)} CTA candidate peptides that also occur in non-CTA proteins.",
    )
    return seen


def cta_exclusive_peptides(
    ensembl_release: int = 112,
    lengths: tuple[int, ...] = DEFAULT_PEPTIDE_LENGTHS,
    gene_names: Iterable[str] | None = None,
    on_progress: Callable[[str], None] | None = None,
    progress_bar: bool = False,
    progress_file: TextIO | None = None,
) -> pd.DataFrame:
    """Return CTA peptides that do NOT appear in any non-CTA protein.

    This filters the output of :func:`cta_peptides` to only include
    peptides whose sequence is exclusive to CTA proteins -- i.e., the
    peptide is not found as a substring of any non-CTA protein-coding
    transcript sequence.

    Parameters
    ----------
    ensembl_release
        Ensembl release (default 112).
    lengths
        Peptide lengths to generate (default 8, 9, 10, 11).
    gene_names
        Optional CTA gene symbols to enumerate. Defaults to all expressed CTAs.
    on_progress
        Optional progress callback.
    progress_bar
        If True, show tqdm bars while walking CTA and non-CTA genes.
    progress_file
        File handle for tqdm output. Defaults to ``sys.stderr``.

    Returns
    -------
    pd.DataFrame
        Same columns as :func:`cta_peptides`, filtered to exclusive peptides.

    Notes
    -----
    The non-CTA background scan is candidate-driven: tsarina first generates
    CTA k-mers, then streams non-CTA protein-coding transcripts and records
    only overlaps with those CTA candidates. This avoids building a full human
    k-mer posting index when the caller only needs to subtract CTA overlaps.
    """
    cta_df = cta_peptides(
        ensembl_release=ensembl_release,
        lengths=lengths,
        gene_names=gene_names,
        on_progress=on_progress,
        progress_bar=progress_bar,
        progress_file=progress_file,
    )
    if cta_df.empty:
        return cta_df

    candidate_peptides = set(cta_df["peptide"])
    noncta_peptides = _non_cta_overlapping_peptides(
        candidate_peptides,
        ensembl_release=ensembl_release,
        lengths=tuple(lengths),
        on_progress=on_progress,
        progress_bar=progress_bar,
        progress_file=progress_file,
    )
    mask = ~cta_df["peptide"].isin(noncta_peptides)
    out = cta_df[mask].reset_index(drop=True)
    _report_progress(
        on_progress,
        f"CTA exclusivity filter kept {len(out)} peptide rows "
        f"({out['peptide'].nunique() if not out.empty else 0} unique peptides).",
    )
    return out


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
