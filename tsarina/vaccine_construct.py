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
    selection_mode: str = "ranked"
    exclude_gene_patterns: tuple[str, ...] = ()
    allow_genes: tuple[str, ...] = ()
    ms_support_mode: str = "presentation"
    ms_affinity_nm: float = 1000
    allow_untyped_ms: bool = False
    normal_ms_policy: str = "audit"
    normal_ms_atlas_dir: str | None = None
    normal_ms_min_donors: int = 1
    species: str = "human"
    canine_cohort: str | None = None
    allow_exploratory_dla: bool = False

    @classmethod
    def for_canine(cls, cohort, **kwargs):
        """Explicit dog defaults without changing the human constructor contract."""
        values = {
            "species": "canine",
            "panel": "bundle",
            "canine_cohort": cohort,
            "ms_support_mode": "sample_affinity",
            "normal_ms_policy": "exclude",
            "codon_species": "generic",
            "predictor": "frozen",
        }
        values.update(kwargs)
        return cls(**values)

    def validate(self):
        if self.species not in {"human", "canine"}:
            raise ValueError("species must be human or canine")
        if self.species == "canine":
            if not self.canine_cohort:
                raise ValueError("Canine design requires an explicit canine_cohort")
            if self.panel != "bundle":
                raise ValueError("Canine design requires the bundle DLA panel or explicit alleles")
            if self.normal_ms_atlas_dir:
                raise ValueError("Human Atlas inputs cannot supply canine normal evidence")
            if self.ms_support_mode != "sample_affinity":
                raise ValueError(
                    "Canine MS support requires sample_affinity, not human percentiles"
                )
            if self.normal_ms_policy != "exclude":
                raise ValueError("Canine mode requires verified healthy-primary MS exclusion")
            if self.predictor != "frozen":
                raise ValueError(
                    "Canine design uses frozen affinity or explicit callbacks, not an implicit human predictor"
                )
            if self.include_utrs and any(
                value.upper() in {"HBB", "HBB_FI"} for value in (self.utr_5p, self.utr_3p)
            ):
                raise ValueError("Canine UTRs require explicit sequences or none")
            if self.require_clean_junctions:
                raise ValueError(
                    "Canine final assessment is unassessed; clean-junction certification unavailable"
                )
        if self.normal_ms_policy not in {"audit", "exclude"}:
            raise ValueError("normal_ms_policy must be audit or exclude")
        if (
            self.normal_ms_policy == "exclude"
            and not self.normal_ms_atlas_dir
            and self.species == "human"
        ):
            raise ValueError(
                "Normal-MS exclusion requires normal_ms_atlas_dir with verified Atlas tables"
            )
        if type(self.normal_ms_min_donors) is not int or self.normal_ms_min_donors < 1:
            raise ValueError("normal_ms_min_donors must be a positive integer")
        if self.selection_mode not in {"ranked", "supported", "budget"}:
            raise ValueError("selection_mode must be ranked, supported or budget")
        if (
            self.selection_mode == "budget"
            and self.max_length_aa is None
            and self.max_length_nt is None
        ):
            raise ValueError("Budget selection requires max_length_aa or max_length_nt")
        if self.ms_support_mode not in {"presentation", "sample_affinity"}:
            raise ValueError("ms_support_mode must be presentation or sample_affinity")
        if not math.isfinite(self.ms_affinity_nm) or self.ms_affinity_nm <= 0:
            raise ValueError("ms_affinity_nm must be finite and positive")
        for name in ("exclude_gene_patterns", "allow_genes"):
            if any(not isinstance(s, str) or not s.strip() for s in getattr(self, name)):
                raise ValueError(f"{name} must contain nonempty gene symbols/patterns")
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
                if self.config.species == "canine":
                    self.cleavage_model = {"name": "unassessed", "species": "canine"}
                    self.profiles.update({c: [None] * len(c) for c in missing_contexts})
                    return
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
        assessed_cleavage = [
            r["cleavage_probability"] for r in cleavage if r["cleavage_probability"] is not None
        ]
        # Failures first; distance from threshold differentiates two equally
        # burdened joins. Cleavage is a secondary predictive design criterion.
        key = (
            len(bad),
            sum(math.log(self.config.junction_affinity_nm / a) for a in bad),
            -sum(assessed_cleavage) / max(1, len(assessed_cleavage)),
            sum(len(layer["sequence"]) for layer in layers if layer["kind"] == "linker"),
            -sum(len(layer["sequence"]) for layer in layers if layer["kind"] == "cta_segment"),
        )
        return key, junctions, cleavage


def construct_fits(placements, config):
    """Check actual initiation and all nucleotide elements at the length caps."""
    utr5, utr3, poly_a = nucleotide_elements(config)
    sequence, _ = assemble_layers(placements)
    return (config.max_length_aa is None or len(sequence) <= config.max_length_aa) and (
        config.max_length_nt is None
        or 3 * len(sequence) + len(utr5) + len(utr3) + len(poly_a) + 3 <= config.max_length_nt
    )


def reserve_target_segments(segments, config):
    """Reserve one shortest whole ligand-bearing segment per target in rank order."""
    groups = {}
    for segment in sorted(segments, key=lambda s: (s.rank, s.proteoform_key)):
        groups.setdefault(segment.proteoform_key, []).append(segment)
    retained, rejected = [], []
    for key, group in groups.items():
        shortest = min(
            group,
            key=lambda s: (
                len(s.sequence(config.min_padding, config.min_padding)),
                -len(s.alleles),
                s.segment_id,
            ),
        )
        placement = (shortest, config.min_padding, config.min_padding, "")
        if construct_fits([*retained, placement], config):
            retained.append(placement)
        else:
            rejected.append(key)
        if len(retained) == config.top_k:
            break
    return retained, rejected


def optimize_construct(
    segments,
    alleles,
    config,
    affinity_fn=None,
    cleavage_fn=None,
    on_progress=None,
    background_fn=None,
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
    if config.selection_mode == "supported":
        retained, rejected = reserve_target_segments(segments, config)
        if rejected or len(retained) != len({s.proteoform_key for s in segments}):
            raise NoFeasibleConstruct("Cannot preserve every selected supported proteoform")
    reserved_ids = {placement[0].segment_id for placement in retained}
    for segment in sorted(segments, key=lambda s: (s.rank, -len(s.alleles), s.segment_id)):
        if segment.segment_id in reserved_ids:
            continue
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
        forbidden = (
            background_fn([sequence for sequence, _ in assemblies]) if background_fn else set()
        )
        if on_progress:
            on_progress(
                f"Construct search round {round_index}: auditing {len(assemblies)} candidates"
            )
        audit.prepare(assemblies)
        scored = []
        for (sig, state), (sequence, layers) in zip(states.items(), assemblies):
            key, _, _ = audit.assess(sequence, layers)
            if background_fn:
                key = (
                    sum(sequence[i : i + 8] in forbidden for i in range(len(sequence) - 7)),
                    *key,
                )
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
    result = {
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
    if background_fn:
        result["initial_background_overlap_windows"] = initial_key[0]
        result["initial_objective"] = initial_key[1:]
        for row in result["search_history"]:
            row["background_overlap_windows"] = row["best_key"][0]
            row["best_key"] = row["best_key"][1:]
        forbidden = background_fn([protein])
        result["background_overlaps"] = [
            {"start": i, "end": i + 8, "peptide": protein[i : i + 8]}
            for i in range(len(protein) - 7)
            if protein[i : i + 8] in forbidden
        ]
        result["background_safe"] = not result["background_overlaps"]
    return result


def encode_construct(design, config):
    """Apply the explicit coding policy and validate exact translation."""
    from Bio.Seq import Seq

    from .vaccine_elements import codon_optimize

    if config.codon_species == "generic":
        from Bio.Data.CodonTable import unambiguous_dna_by_id

        codons = {}
        for codon, aa in sorted(unambiguous_dna_by_id[1].forward_table.items()):
            codons.setdefault(aa, codon)
        dna = "".join(codons[aa] for aa in design["protein"]) + "TAA"
    else:
        dna = codon_optimize(design["protein"], species=config.codon_species).upper() + "TAA"
    if not dna.startswith("ATG") or str(Seq(dna).translate()) != design["protein"] + "*":
        raise ValueError("Codon-optimized construct failed exact translation validation")
    full = design["utr5"] + dna + design["utr3"] + design["poly_a"]
    if config.max_length_nt is not None and len(full) > config.max_length_nt:
        raise ValueError("Encoded construct exceeds total nucleotide limit")
    if config.vaccine_type == "rna":
        dna, full = dna.replace("T", "U"), full.replace("T", "U")
    return dna, full
