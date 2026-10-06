"""Native interval subtraction and ligand-preserving vaccine segments.

Coordinates throughout the Python API and CSVs are zero-based, half-open.
"""

from __future__ import annotations

from dataclasses import dataclass

from .peptides import AA20


def shared_kmers(proteins, cta_ids, selected_sequences, k=8):
    """Only materialize relevant k-mers, scanning every non-CTA isoform.

    Identical-sequence background occurrences remain disqualifying even when
    the same full sequence is also encoded by a CTA gene.
    """
    wanted = {s[i : i + k] for s in selected_sequences for i in range(len(s) - k + 1)}
    remaining = wanted.copy()
    for protein in proteins:
        if protein.gene_id in cta_ids:
            continue
        seq = protein.sequence
        for i in range(len(seq) - k + 1):
            remaining.discard(seq[i : i + k])
        if not remaining:
            break
    return wanted - remaining


def specific_intervals(sequence, forbidden_kmers, k=8):
    """Complement of the union of residues in shared k-mers or ambiguous AA."""
    blocked = [aa not in AA20 for aa in sequence]
    for i in range(len(sequence) - k + 1):
        if sequence[i : i + k] in forbidden_kmers:
            blocked[i : i + k] = [True] * k
    intervals, start = [], None
    for i, excluded in enumerate([*blocked, True]):
        if not excluded and start is None:
            start = i
        elif excluded and start is not None:
            intervals.append((start, i))
            start = None
    return intervals


def peptide_occurrences(sequence, peptide):
    """All exact occurrences, including overlapping repeats."""
    start = sequence.find(peptide)
    while start >= 0:
        yield start, start + len(peptide)
        start = sequence.find(peptide, start + 1)


@dataclass(frozen=True)
class Segment:
    segment_id: str
    proteoform_key: str
    name: str
    rank: int
    score: float
    protein_sequence: str
    specific_start: int
    specific_end: int
    ligand_start: int
    ligand_end: int
    alleles: tuple[str, ...]
    peptides: tuple[str, ...]

    def bounds(self, n_padding, c_padding):
        return (
            max(self.specific_start, self.ligand_start - n_padding),
            min(self.specific_end, self.ligand_end + c_padding),
        )

    def sequence(self, n_padding, c_padding):
        start, end = self.bounds(n_padding, c_padding)
        return self.protein_sequence[start:end]


def supported_segments(selected, intervals, support):
    """Retain native pieces with qualifying MS support; preserve all occurrences."""
    segments, ligands = [], []
    for protein in selected.itertuples(index=False):
        for i, (start, end) in enumerate(intervals[protein.proteoform_key], 1):
            matches = []
            for hit in support.itertuples(index=False):
                for a, b in peptide_occurrences(protein.sequence[start:end], hit.peptide):
                    matches.append((start + a, start + b, hit))
            if not matches:
                continue
            sid = f"{protein.proteoform_key}:{start}-{end}"
            segments.append(
                Segment(
                    sid,
                    protein.proteoform_key,
                    protein.name,
                    int(protein.rank),
                    float(protein.mortality_weighted_score),
                    protein.sequence,
                    start,
                    end,
                    min(a for a, _, _ in matches),
                    max(b for _, b, _ in matches),
                    tuple(sorted({h.allele for _, _, h in matches})),
                    tuple(sorted({h.peptide for _, _, h in matches})),
                )
            )
            for a, b, hit in matches:
                ligands.append(
                    {
                        "proteoform_key": protein.proteoform_key,
                        "name": protein.name,
                        "segment_id": sid,
                        "specific_piece": i,
                        "start": a,
                        "end": b,
                        **hit._asdict(),
                    }
                )
    return segments, ligands


def assemble_layers(placements):
    """Assemble (Segment, n padding, c padding, incoming linker) placements."""
    layers, parts, offset = [], [], 0

    def emit(kind, seq, **meta):
        nonlocal offset
        if not seq:
            return
        layers.append(
            {"kind": kind, "start_aa": offset, "end_aa": offset + len(seq), "sequence": seq, **meta}
        )
        parts.append(seq)
        offset += len(seq)

    if placements and not placements[0][0].sequence(*placements[0][1:3]).startswith("M"):
        emit("start_methionine", "M")
    for i, (segment, n_pad, c_pad, linker) in enumerate(placements):
        if i:
            emit("linker", linker)
        start, end = segment.bounds(n_pad, c_pad)
        emit(
            "cta_segment",
            segment.protein_sequence[start:end],
            segment_id=segment.segment_id,
            proteoform_key=segment.proteoform_key,
            name=segment.name,
            native_start=start,
            native_end=end,
            n_padding=segment.ligand_start - start,
            c_padding=end - segment.ligand_end,
        )
    return "".join(parts), layers


def junction_windows(sequence, layers, lengths=(8, 9, 10, 11)):
    """Every boundary-spanning window, including linker-only/multi-boundary ones."""
    boundaries = [layer["start_aa"] for layer in layers[1:]]
    rows = []
    # Windows wholly inside a linker must also be checked: they are synthetic.
    for k in lengths:
        for start in range(len(sequence) - k + 1):
            end = start + k
            crossed = [b for b in boundaries if start < b < end]
            in_linker = any(
                layer["kind"] == "linker" and layer["start_aa"] <= start and end <= layer["end_aa"]
                for layer in layers
            )
            if crossed or in_linker:
                rows.append(
                    {
                        "peptide": sequence[start:end],
                        "start": start,
                        "end": end,
                        "boundaries": ";".join(map(str, crossed)),
                        "inside_linker": in_linker,
                    }
                )
    return rows
