"""Bounded construct search with actual sequence-level binding/cleavage audits."""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass

from .peptides import AA20
from .vaccine_sequences import assemble_layers, junction_windows


class NoFeasibleConstruct(ValueError):
    """Supported native pieces cannot satisfy the requested construct limits."""


@dataclass(frozen=True)
class VaccineConfig:
    top_k: int = 10
    panel: str = "global54_abc"
    alleles: tuple[str, ...] | None = None
    definition: str = "strict"
    ensembl_release: int = 112
    shared_k: int = 8
    min_padding: int = 0
    max_padding: int = 10
    padding_step: int = 5
    lengths: tuple[int, ...] = (8, 9, 10, 11)
    predictor: str = "mhcflurry"
    junction_affinity_nm: float = 1000
    linkers: tuple[str, ...] = ("", "AAY")
    beam_width: int = 4
    optimization_rounds: int = 3
    max_length_aa: int | None = None
    max_length_nt: int | None = None
    vaccine_type: str = "rna"
    include_utrs: bool = False
    utr_5p: str = "HBB"
    utr_3p: str = "HBB_FI"
    poly_a_length: int = 0
    require_clean_junctions: bool = False
    codon_species: str = "h_sapiens"

    def validate(self):
        for name in ("top_k", "shared_k", "padding_step", "beam_width"):
            if type(getattr(self, name)) is not int or getattr(self, name) < 1:
                raise ValueError(f"{name} must be a positive integer")
        if self.shared_k != 8:
            raise ValueError("CTA-specific subtraction currently requires shared_k=8")
        for name in ("min_padding", "max_padding", "poly_a_length", "optimization_rounds"):
            if type(getattr(self, name)) is not int or getattr(self, name) < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
        if self.min_padding > self.max_padding:
            raise ValueError("min_padding cannot exceed max_padding")
        if not self.lengths or any(type(k) is not int or not 8 <= k <= 15 for k in self.lengths):
            raise ValueError("Class-I ligand lengths must be in 8..15")
        if self.definition not in {"strict", "loose"} or self.vaccine_type not in {"dna", "rna"}:
            raise ValueError("Invalid CTA definition or vaccine type")
        if not math.isfinite(self.junction_affinity_nm) or self.junction_affinity_nm <= 0:
            raise ValueError("junction_affinity_nm must be finite and positive")
        for name in ("max_length_aa", "max_length_nt"):
            value = getattr(self, name)
            if value is not None and (type(value) is not int or value < 1):
                raise ValueError(f"{name} must be positive")
        if not self.linkers or any(set(linker) - AA20 for linker in self.linkers):
            raise ValueError("Linkers must be strings containing canonical amino acids")
        if "" not in self.linkers:
            raise ValueError("Linker candidates must include a direct join (empty string)")


def resolve_utr(value, which):
    from .vaccine_elements import UTRS_3P, UTRS_5P

    table = UTRS_5P if which == "5p" else UTRS_3P
    if value.lower() == "none":
        return ""
    seq = table.get(value, value).upper().replace("U", "T")
    if not seq or set(seq) - set("ACGT"):
        raise ValueError(f"Unknown {which} UTR {value!r}; use a Vaxrank name or A/C/G/T/U sequence")
    return seq


def nucleotide_elements(config):
    return (
        resolve_utr(config.utr_5p, "5p") if config.include_utrs else "",
        resolve_utr(config.utr_3p, "3p") if config.include_utrs else "",
        "A" * config.poly_a_length,
    )


class PredictionAudit:
    """Cache complete peptide/allele and local cleavage predictions across search."""

    def __init__(self, alleles, config, affinity_fn=None, cleavage_fn=None):
        self.alleles, self.config = alleles, config
        if affinity_fn is None:
            from .scoring import score_affinity

            def affinity_fn(peptides, alleles):
                return score_affinity(peptides, alleles, predictor=config.predictor)

        self.affinity_fn = affinity_fn
        self.affinities, self.profiles = {}, {}
        self.cleavage_fn = cleavage_fn
        self.cleavage_model = {"name": "injected"} if cleavage_fn else {}

    def prepare(self, assemblies):
        peptides, contexts = set(), set()
        for sequence, layers in assemblies:
            peptides.update(
                row["peptide"] for row in junction_windows(sequence, layers, self.config.lengths)
            )
            for layer in layers[1:]:
                b = layer["start_aa"]
                # Pepsickle uses 8 preceding residues + the cleavage residue
                # itself + 8 following residues (17, not 16, total).
                contexts.add(sequence[max(0, b - 9) : b + 8])
        missing = sorted(peptides - self.affinities.keys())
        if missing:
            frame = self.affinity_fn(missing, self.alleles)
            if frame.duplicated(["peptide", "allele"]).any():
                raise ValueError("Duplicate junction affinity predictions")
            lookup = frame.set_index(["peptide", "allele"])
            for peptide in missing:
                scores = {}
                for allele in self.alleles:
                    if (peptide, allele) not in lookup.index:
                        raise ValueError(f"Missing junction affinity for {peptide}/{allele}")
                    value = float(lookup.loc[(peptide, allele), "affinity_nm"])
                    if not math.isfinite(value) or value <= 0:
                        raise ValueError(f"Invalid junction affinity for {peptide}/{allele}")
                    scores[allele] = value
                self.affinities[peptide] = scores
        missing_contexts = sorted(contexts - self.profiles.keys())
        if missing_contexts:
            if self.cleavage_fn is None:
                from mhctools import Pepsickle

                predictor = Pepsickle(human_only=True, isolate_subprocess=True)
                self.cleavage_model = asdict(predictor.cleavage_model())
                self.cleavage_fn = predictor.cleavage_probs_many
            profiles = self.cleavage_fn(missing_contexts)
            for context in missing_contexts:
                profile = profiles[context]
                if len(profile) != len(context) or any(
                    not math.isfinite(float(x)) or not 0 <= x <= 1 for x in profile
                ):
                    raise ValueError("Invalid proteasomal cleavage predictions")
                self.profiles[context] = list(map(float, profile))

    def assess(self, sequence, layers):
        junctions, cleavage = [], []
        for window in junction_windows(sequence, layers, self.config.lengths):
            for allele, affinity in self.affinities[window["peptide"]].items():
                junctions.append(
                    {
                        **window,
                        "allele": allele,
                        "affinity_nm": affinity,
                        "below_threshold": affinity < self.config.junction_affinity_nm,
                    }
                )
        for layer in layers[1:]:
            b = layer["start_aa"]
            start = max(0, b - 9)
            context = sequence[start : b + 8]
            cleavage.append(
                {
                    "bond": b,
                    "context": context,
                    "right_kind": layer["kind"],
                    "cleavage_probability": self.profiles[context][b - start - 1],
                }
            )
        bad = [r["affinity_nm"] for r in junctions if r["below_threshold"]]
        # Failures first; distance from threshold differentiates two equally
        # burdened joins. Cleavage is a secondary predictive design criterion.
        key = (
            len(bad),
            sum(math.log(self.config.junction_affinity_nm / a) for a in bad),
            -sum(r["cleavage_probability"] for r in cleavage) / max(1, len(cleavage)),
            sum(len(layer["sequence"]) for layer in layers if layer["kind"] == "linker"),
            -sum(len(layer["sequence"]) for layer in layers if layer["kind"] == "cta_segment"),
        )
        return key, junctions, cleavage


def optimize_construct(
    segments, alleles, config, affinity_fn=None, cleavage_fn=None, on_progress=None
):
    """Bounded beam of complete constructs, preserving every retained ligand.

    Each round explores adjacent swaps, reversal/rotations, independent edge
    padding changes, and linkers. Full-sequence windows avoid the false
    independence assumption for adjacent junctions near short segments.
    This deterministic heuristic does not claim a global optimum.
    """
    utr5, utr3, poly_a = nucleotide_elements(config)
    overhead_nt = len(utr5) + len(utr3) + len(poly_a) + 3  # stop codon
    audit = PredictionAudit(alleles, config, affinity_fn, cleavage_fn)

    def feasible(placements):
        sequence, _ = assemble_layers(placements)
        return (config.max_length_aa is None or len(sequence) <= config.max_length_aa) and (
            config.max_length_nt is None or 3 * len(sequence) + overhead_nt <= config.max_length_nt
        )

    retained, excluded = [], []
    for segment in sorted(segments, key=lambda s: (s.rank, -len(s.alleles), s.segment_id)):
        placement = (segment, config.min_padding, config.min_padding, "")
        if feasible([*retained, placement]):
            retained.append(placement)
        else:
            excluded.append({"segment_id": segment.segment_id, "reason": "construct_length_limit"})
    if not retained:
        raise NoFeasibleConstruct("No MS-supported segment fits the construct length constraints")
    # Start with maximal context where it fits; it is then free to shrink.
    for i, (segment, _, _, linker) in enumerate(retained):
        trial = retained.copy()
        trial[i] = (segment, config.max_padding, config.max_padding, linker)
        if feasible(trial):
            retained = trial
    padding = sorted(
        {
            *range(config.min_padding, config.max_padding + 1, config.padding_step),
            config.max_padding,
        }
    )
    beam = [retained]
    seen = set()
    history = []

    def signature(state):
        # Padding beyond a native interval end produces the same construct.
        # Count actual boundaries so these aliases cannot fill the beam.
        return tuple(
            (s.segment_id, *s.bounds(n, c), link if i else "")
            for i, (s, n, c, link) in enumerate(state)
        )

    def candidates(state):
        yield state
        if len(state) > 1:
            yield list(reversed(state))
            yield [*state[1:], state[0]]
            yield [state[-1], *state[:-1]]
        for i, (segment, n, c, link) in enumerate(state):
            if i:
                trial = state.copy()
                trial[i - 1], trial[i] = trial[i], trial[i - 1]
                yield trial
                for linker in config.linkers:
                    trial = state.copy()
                    trial[i] = (segment, n, c, linker)
                    yield trial
                    if len(linker) > len(link) and not feasible(trial):
                        # At the length cap, a linker and the padding needed
                        # to fit it must be proposed together. Requiring a
                        # worse direct join first can trap a narrow beam.
                        left, left_n, left_c, left_link = state[i - 1]
                        for cpad in padding:
                            for npad in padding:
                                if cpad > left_c or npad > n:
                                    continue
                                trial = state.copy()
                                trial[i - 1] = (left, left_n, cpad, left_link)
                                trial[i] = (segment, npad, c, linker)
                                yield trial
            for npad in padding:
                for cpad in padding:
                    trial = state.copy()
                    trial[i] = (segment, npad, cpad, link)
                    yield trial

    initial_key = None
    for round_index in range(config.optimization_rounds + 1):
        states = {}
        for state in beam:
            for trial in [state] if round_index == 0 else candidates(state):
                # Incoming linker belongs to the boundary, not the source protein.
                trial = [(s, n, c, link if i else "") for i, (s, n, c, link) in enumerate(trial)]
                sig = signature(trial)
                if feasible(trial) and (sig not in seen or trial in beam):
                    states[sig] = trial
        assemblies = [assemble_layers(state) for state in states.values()]
        if on_progress:
            on_progress(
                f"Construct search round {round_index}: auditing {len(assemblies)} candidates"
            )
        audit.prepare(assemblies)
        scored = []
        for (sig, state), (sequence, layers) in zip(states.items(), assemblies):
            key, _, _ = audit.assess(sequence, layers)
            scored.append((key, sig, state))
            seen.add(sig)
        scored.sort(key=lambda row: (row[0], row[1]))
        if initial_key is None:
            initial_key = scored[0][0]
        history.append({"round": round_index, "candidates": len(scored), "best_key": scored[0][0]})
        beam = [state for _, _, state in scored[: config.beam_width]]
    chosen = beam[0]
    protein, layers = assemble_layers(chosen)
    # Reassess the complete final product, not just the pairwise joins.
    audit.prepare([(protein, layers)])
    key, junctions, cleavage = audit.assess(protein, layers)
    return {
        "protein": protein,
        "layers": layers,
        "junctions": junctions,
        "cleavage": cleavage,
        "excluded_segments": excluded,
        "search_history": history,
        "initial_objective": initial_key,
        "final_objective": key,
        "clean_junctions": key[0] == 0,
        "cleavage_model": audit.cleavage_model,
        "predicted_unique_peptides": len(audit.affinities),
        "utr5": utr5,
        "utr3": utr3,
        "poly_a": poly_a,
    }


def encode_construct(design, config):
    """Use Vaxrank's human codon optimization; validate exact translation."""
    from Bio.Seq import Seq

    from .vaccine_elements import codon_optimize

    dna = codon_optimize(design["protein"], species=config.codon_species).upper() + "TAA"
    if not dna.startswith("ATG") or str(Seq(dna).translate()) != design["protein"] + "*":
        raise ValueError("Codon-optimized construct failed exact translation validation")
    full = design["utr5"] + dna + design["utr3"] + design["poly_a"]
    if config.max_length_nt is not None and len(full) > config.max_length_nt:
        raise ValueError("Encoded construct exceeds total nucleotide limit")
    if config.vaccine_type == "rna":
        dna, full = dna.replace("T", "U"), full.replace("T", "U")
    return dna, full
