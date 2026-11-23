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

"""Curated cancer fusions and junction-spanning peptide enumeration.

Wraps OncoRef's curated cancer fusion catalog (:func:`oncoref.cancer_fusions_df`)
for discovery and provides breakpoint-spanning k-mer enumeration with
reference-proteome self-exclusivity screening so fusion neoantigens can be
fed directly into MHC scoring and :func:`tsarina.personalized_targets`.

Typical usage::

    from tsarina import cancer_fusions, fusion_junction_peptides

    # Discover curated fusions in OncoRef
    ewing = cancer_fusions(cancer_code="SARC_EWS")

    # Enumerate human-exclusive junction peptides from an explicit breakpoint
    peps = fusion_junction_peptides(
        [
            {
                "label": "EWSR1-FLI1",
                "gene_5prime": "EWSR1",
                "gene_3prime": "FLI1",
                "n_terminal_seq": "SYGQQSSYPGYGQQ",
                "c_terminal_seq": "NPSYDSVRRGAWGN",
            }
        ],
        lengths=(9,),
    )
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from .peptides import AA20, _reference_overlapping_peptides
from .viral import DEFAULT_PEPTIDE_LENGTHS, _proteome_kmers_from_sequences

FUSION_PEPTIDE_COLUMNS: tuple[str, ...] = (
    "label",
    "gene_5prime",
    "gene_3prime",
    "peptide",
    "length",
    "n_terminal_residues",
    "c_terminal_residues",
    "junction_offset",
    "in_reference_proteome",
)


def cancer_fusions(
    *,
    cancer_code: str | None = None,
    gene: str | None = None,
    fusion_family: str | None = None,
    defining_only: bool = False,
    pathognomonic_only: bool = False,
    rnaseq_detectable_only: bool = False,
    mechanism: str | None = None,
) -> pd.DataFrame:
    """Return curated cancer fusions from OncoRef with optional filters.

    Parameters
    ----------
    cancer_code
        Optional case-insensitive cancer code filter (e.g. ``"SARC_EWS"``,
        ``"LUAD"``).
    gene
        Optional case-insensitive gene symbol matched against either
        ``gene_5prime`` or ``gene_3prime`` (e.g. ``"ALK"``, ``"EWSR1"``).
    fusion_family
        Optional case-insensitive ``fusion_family`` filter (e.g.
        ``"EWSR1-ETS"``, ``"SS18-SSX"``).
    defining_only
        When ``True``, keep only fusions with ``is_defining == True``.
    pathognomonic_only
        When ``True``, keep only fusions with ``pathognomonic == True``.
    rnaseq_detectable_only
        When ``True``, keep only fusions with ``rnaseq_detectable == "yes"``.
    mechanism
        Optional case-insensitive mechanism filter (e.g.
        ``"fusion_transcript"``).

    Returns
    -------
    pd.DataFrame
        Filtered copy of :func:`oncoref.cancer_fusions_df`.
    """
    import oncoref

    df = oncoref.cancer_fusions_df().copy()
    if df.empty:
        return df

    if cancer_code is not None:
        target = cancer_code.strip().lower()
        df = df[df["cancer_code"].astype(str).str.lower() == target]

    if gene is not None:
        target = gene.strip().lower()
        df = df[
            (df["gene_5prime"].astype(str).str.lower() == target)
            | (df["gene_3prime"].astype(str).str.lower() == target)
        ]

    if fusion_family is not None:
        target = fusion_family.strip().lower()
        df = df[df["fusion_family"].astype(str).str.lower() == target]

    if defining_only and "is_defining" in df.columns:
        df = df[df["is_defining"].astype(bool)]

    if pathognomonic_only and "pathognomonic" in df.columns:
        df = df[df["pathognomonic"].astype(bool)]

    if rnaseq_detectable_only and "rnaseq_detectable" in df.columns:
        df = df[df["rnaseq_detectable"].astype(str).str.lower() == "yes"]

    if mechanism is not None and "mechanism" in df.columns:
        target = mechanism.strip().lower()
        df = df[df["mechanism"].astype(str).str.lower() == target]

    return df.reset_index(drop=True)


def _normalize_fusion_records(
    fusions: Mapping[str, object] | Sequence[Mapping[str, object]] | pd.DataFrame,
) -> list[dict[str, object]]:
    if isinstance(fusions, pd.DataFrame):
        return [dict(row) for row in fusions.to_dict(orient="records")]
    if isinstance(fusions, Mapping):
        return [dict(fusions)]
    return [dict(item) for item in fusions]


def _extract_breakpoint(record: Mapping[str, object]) -> tuple[str, int]:
    has_split = "n_terminal_seq" in record or "c_terminal_seq" in record
    has_joined = "sequence" in record or "junction_pos" in record
    if has_split and has_joined:
        raise ValueError(
            "Provide either (n_terminal_seq, c_terminal_seq) or (sequence, junction_pos), not both."
        )
    if has_split:
        n_seq = str(record.get("n_terminal_seq") or "").strip().upper()
        c_seq = str(record.get("c_terminal_seq") or "").strip().upper()
        if not n_seq or not c_seq:
            raise ValueError("Both n_terminal_seq and c_terminal_seq must be non-empty strings.")
        return n_seq + c_seq, len(n_seq)
    if has_joined:
        seq = str(record.get("sequence") or "").strip().upper()
        pos_raw = record.get("junction_pos")
        if not seq or pos_raw is None:
            raise ValueError("Both sequence and junction_pos are required.")
        pos = int(pos_raw)
        if pos <= 0 or pos >= len(seq):
            raise ValueError(
                f"junction_pos must satisfy 1 <= junction_pos < len(sequence); got {pos}."
            )
        return seq, pos
    raise ValueError(
        "Each fusion record must provide (n_terminal_seq, c_terminal_seq) "
        "or (sequence, junction_pos)."
    )


def fusion_junction_peptides(
    fusions: Mapping[str, object] | Sequence[Mapping[str, object]] | pd.DataFrame,
    *,
    lengths: tuple[int, ...] = DEFAULT_PEPTIDE_LENGTHS,
    ensembl_release: int = 112,
    require_host_exclusive: bool = True,
    host_species: str = "human",
    host_proteins: dict[str, str] | list[str] | tuple[str, ...] | None = None,
) -> pd.DataFrame:
    """Enumerate junction-spanning k-mer peptides for explicit fusion sequences.

    Every emitted peptide contains at least one residue from the 5' partner
    (N-terminal side of the breakpoint) and at least one residue from the 3'
    partner (C-terminal side of the breakpoint), restricted to the 20 standard
    amino acids.

    By default (``require_host_exclusive=True``), peptides that match any
    protein in the reference host proteome are excluded so only non-self
    junction candidates are returned.

    Parameters
    ----------
    fusions
        Single fusion dict, sequence of fusion dicts, or DataFrame. Each entry
        must specify either ``(n_terminal_seq, c_terminal_seq)`` or
        ``(sequence, junction_pos)`` (0-based index of the first 3' residue),
        plus optional ``label``, ``gene_5prime``, and ``gene_3prime``.
    lengths
        Peptide lengths to enumerate (default ``(8, 9, 10, 11)``).
    ensembl_release
        Ensembl release used for reference host-proteome screening when
        ``host_proteins`` is not supplied (default ``112``).
    require_host_exclusive
        When ``True`` (default), drop junction peptides that occur anywhere in
        the reference host proteome. When ``False``, retain all rows and flag
        matches via ``in_reference_proteome``.
    host_species
        Host species name for reference proteome screening (default ``"human"``).
    host_proteins
        Optional explicit host protein sequences (dict or sequence of AA strings)
        used instead of Ensembl for offline screening or custom backgrounds.

    Returns
    -------
    pd.DataFrame
        Columns: :data:`FUSION_PEPTIDE_COLUMNS`.
    """
    records = _normalize_fusion_records(fusions)
    if not records:
        return pd.DataFrame(columns=list(FUSION_PEPTIDE_COLUMNS))

    lengths_t = tuple(int(k) for k in lengths)
    if any(k < 2 for k in lengths_t):
        raise ValueError("Peptide lengths for junction enumeration must be >= 2.")

    rows: list[dict[str, object]] = []
    for rec in records:
        full_seq, junction_pos = _extract_breakpoint(rec)
        g5 = str(rec.get("gene_5prime") or "").strip()
        g3 = str(rec.get("gene_3prime") or "").strip()
        default_label = f"{g5}-{g3}" if g5 and g3 else (g5 or g3 or "fusion")
        label = str(rec.get("label") or default_label).strip()

        seq_len = len(full_seq)
        for k in lengths_t:
            start_min = max(0, junction_pos - k + 1)
            start_max = min(junction_pos, seq_len - k + 1)
            for start in range(start_min, start_max):
                pep = full_seq[start : start + k]
                if not set(pep).issubset(AA20):
                    continue
                n_res = junction_pos - start
                c_res = start + k - junction_pos
                rows.append(
                    {
                        "label": label,
                        "gene_5prime": g5,
                        "gene_3prime": g3,
                        "peptide": pep,
                        "length": k,
                        "n_terminal_residues": n_res,
                        "c_terminal_residues": c_res,
                        "junction_offset": n_res + 1,
                        "in_reference_proteome": False,
                    }
                )

    if not rows:
        return pd.DataFrame(columns=list(FUSION_PEPTIDE_COLUMNS))

    out = pd.DataFrame(rows)
    candidate_peptides = set(out["peptide"].unique())
    used_lengths = tuple(sorted({int(k) for k in out["length"].unique()}))

    if host_proteins is not None:
        host_kmers = _proteome_kmers_from_sequences(host_proteins, used_lengths)
        overlapping = candidate_peptides & host_kmers
    else:
        from .species import resolve_species_context

        ctx = resolve_species_context(host_species, ensembl_release=ensembl_release)
        if ctx.is_human_host:
            overlapping = _reference_overlapping_peptides(
                candidate_peptides,
                ensembl_release=ctx.ensembl_release,
                lengths=used_lengths,
            )
        else:
            from hitlist.proteome import proteome_kmer_set

            host_kmers = proteome_kmer_set(
                release=ctx.ensembl_release,
                lengths=used_lengths,
                gene_ids=None,
                species=ctx.pyensembl_species,
            )
            overlapping = candidate_peptides & host_kmers

    out["in_reference_proteome"] = out["peptide"].isin(overlapping)
    if require_host_exclusive:
        out = out[~out["in_reference_proteome"]].copy()

    return out[list(FUSION_PEPTIDE_COLUMNS)].reset_index(drop=True)
