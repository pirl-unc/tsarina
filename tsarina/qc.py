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

"""Quality-control checks for the curated CTA universe.

The headline check is :func:`find_fragment_gene_models`, which flags candidate
CTA Ensembl gene IDs whose annotated protein is a degenerate **fragment** rather
than a real antigen.  This is the GAGE12B failure mode (tsarina#108): the ID
``ENSG00000236737`` is a 117-bp, single-exon locus with a 7-residue ORF (the
shared GAGE C-terminus), so its HPA gene-level nTPM is quantification noise on a
malformed model rather than evidence of a real cancer-testis antigen.

Detection uses a conservative **absolute** protein-length floor rather than the
"< 50% of paralog-family median" heuristic first proposed in tsarina#108.  Real
CTAs can be genuinely small -- the shortest in the bundled table is PRM1
(protamine 1) at 51 aa, with several SPANX/XAGE members at 72-81 aa -- so a
family-relative rule risks flagging legitimate small antigens, and grouping
paralog "families" from symbols alone is unreliable.  A floor well below the
smallest real CTA (51 aa) and well above a true fragment (GAGE12B, 7 aa) cleanly
separates the two without that ambiguity.

The length lookup is injectable (``length_fn``) so the flagging logic is unit-
testable without pyensembl reference data, which CI does not download.
"""

from __future__ import annotations

from typing import Callable

import pandas as pd

#: Minimum plausible protein length (aa) for a real CTA.  Below this a candidate
#: is treated as a fragment / mis-annotated model.  Chosen conservatively: the
#: smallest legitimate CTA in the bundled table is PRM1 at 51 aa, and the
#: canonical fragment (GAGE12B, tsarina#108) is 7 aa, so this floor flags gross
#: fragments without false-positiving on small-but-real antigens.
MIN_CTA_PROTEIN_AA = 30


def gene_max_protein_length(
    gene_id: str,
    *,
    ensembl_release: int = 112,
    ensembl=None,
) -> int | None:
    """Length (aa) of the longest protein-coding transcript for *gene_id*.

    Returns ``None`` if the gene is unknown or has no translatable
    protein-coding transcript.  Pass a pre-built ``ensembl``
    (:class:`pyensembl.EnsemblRelease`) to amortize the lookup across many
    genes.
    """
    if ensembl is None:
        from pyensembl import EnsemblRelease

        ensembl = EnsemblRelease(ensembl_release)

    unversioned = str(gene_id).split(".")[0].strip()
    try:
        gene = ensembl.gene_by_id(unversioned)
    except ValueError:
        return None

    best = 0
    for transcript in gene.transcripts:
        if transcript.biotype != "protein_coding":
            continue
        try:
            protein = transcript.protein_sequence
        except (ValueError, KeyError, TypeError):
            # pyensembl raises these for transcripts lacking a usable
            # translation; skip them rather than letting one bad transcript
            # mask the gene's real protein-coding transcripts.
            continue
        if protein and len(protein) > best:
            best = len(protein)
    return best or None


def find_fragment_gene_models(
    df: pd.DataFrame | None = None,
    *,
    min_protein_aa: int = MIN_CTA_PROTEIN_AA,
    length_fn: Callable[[str], int | None] | None = None,
    ensembl_release: int = 112,
) -> list[dict]:
    """Flag CTA-table gene IDs whose annotated protein is a fragment.

    Parameters
    ----------
    df
        CTA evidence table (defaults to the bundled table via
        :func:`tsarina.loader.cta_dataframe`).  Must have ``Symbol`` and
        ``Ensembl_Gene_ID`` columns; the latter may hold ``;``-separated IDs.
    min_protein_aa
        Protein-length floor; IDs strictly below it are flagged ``"fragment"``.
    length_fn
        ``gene_id -> length_aa | None`` lookup.  Injectable so the flagging
        logic is testable without pyensembl data.  Defaults to a pyensembl-
        backed lookup at ``ensembl_release``.
    ensembl_release
        Ensembl release for the default ``length_fn``.

    Returns
    -------
    list[dict]
        One entry per flagged ID: ``{"symbol", "gene_id", "protein_length",
        "reason"}``.  ``reason`` is ``"fragment"`` (length below the floor) or
        ``"no_protein"`` (no protein-coding transcript at all).  Empty list when
        every ID resolves to a plausible protein.
    """
    if df is None:
        from .loader import cta_dataframe

        df = cta_dataframe()

    if length_fn is None:
        from pyensembl import EnsemblRelease

        ensembl = EnsemblRelease(ensembl_release)

        def length_fn(gene_id: str) -> int | None:
            return gene_max_protein_length(gene_id, ensembl=ensembl)

    flagged: list[dict] = []
    for symbol, id_cell in df[["Symbol", "Ensembl_Gene_ID"]].itertuples(index=False):
        for raw_id in str(id_cell).split(";"):
            gene_id = raw_id.split(".")[0].strip()
            if not gene_id:
                continue
            length = length_fn(gene_id)
            if length is None:
                flagged.append(
                    {
                        "symbol": symbol,
                        "gene_id": gene_id,
                        "protein_length": None,
                        "reason": "no_protein",
                    }
                )
            elif length < min_protein_aa:
                flagged.append(
                    {
                        "symbol": symbol,
                        "gene_id": gene_id,
                        "protein_length": length,
                        "reason": "fragment",
                    }
                )
    return flagged
