"""Minimal Vaxrank-compatible untranslated elements and coding design.

HBB sequences and nomenclature match openvax/vaxrank's Apache-2.0
``mrna_library.py``. Primary source: RefSeq NM_000518.5, positions 1-50
and 495-628. Tandem 3' HBB: Sahin et al., Nature (2017),
https://doi.org/10.1038/nature23003; Holtkamp et al., Blood (2006),
https://doi.org/10.1182/blood-2006-04-015024.

Keep this small downstream adapter independent of Vaxrank's complete variant
pipeline, which pins a conflicting OncoRef data release. Do not add unsourced
UTRs or make stronger claims than the reference construct literature.
"""

UTR_5P_HBB = "ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACC"
UTR_3P_HBB = (
    "GCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTAC"
    "TAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATT"
    "TATTTTCATTGCAA"
)
UTRS_5P = {"HBB": UTR_5P_HBB}
UTRS_3P = {"HBB": UTR_3P_HBB, "HBB_FI": UTR_3P_HBB * 2}


def codon_optimize(amino_acids, species="h_sapiens"):
    """Vaxrank's DnaChisel use_best_codon method with exact translation constraint."""
    from dnachisel import (
        CodonOptimize,
        DnaOptimizationProblem,
        EnforceTranslation,
        reverse_translate,
    )

    problem = DnaOptimizationProblem(
        sequence=reverse_translate(amino_acids),
        constraints=[EnforceTranslation()],
        objectives=[CodonOptimize(species=species, method="use_best_codon")],
        logger=None,
    )
    problem.resolve_constraints()
    problem.optimize()
    return problem.sequence
