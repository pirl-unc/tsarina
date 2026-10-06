"""Versioned scientific inputs for mortality-prioritized CTA vaccine design."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from hashlib import sha256
from importlib.metadata import version

import pandas as pd

# Broad histology cohorts, following OncoRef's mortality coverage report. No
# narrow subtype is used to stand in for an entire mortality category.
DEFAULT_CANCER_COHORTS = {
    "lung": ["LUAD", "LUSC"],
    "colorectal": ["COAD", "READ"],
    "liver": ["LIHC"],
    "breast": ["BRCA"],
    "stomach": ["STAD"],
    "pancreas": ["PAAD"],
    "esophagus": ["ESCA"],
    "head_and_neck": ["HNSC"],
    "prostate": ["PRAD"],
    "cervix": ["CESC"],
    "bladder": ["BLCA"],
    "kidney": ["KICH", "KIRC", "KIRP"],
    "melanoma": ["SKCM"],
    "ovary": ["OV"],
    "uterus_endometrium": ["UCEC"],
    "thyroid": ["THCA"],
    "brain_cns": ["GBM", "LGG"],
    "non_hodgkin_lymphoma": ["DLBC"],
    "leukemia_AML": ["LAML"],
    "mesothelioma": ["MESO"],
    "testicular_germ_cell": ["TGCT"],
    "adrenal": ["ACC"],
}


@dataclass(frozen=True)
class Protein:
    """One translated occurrence; all isoforms participate in the background."""

    gene_id: str
    gene_name: str
    protein_id: str
    sequence: str
    transcript_id: str = ""
    canonical: bool = True
    contig: str = ""


@dataclass
class VaccineInputs:
    """Portable inputs; prevalence must already be summed/ranked as proteoforms.

    ``prevalence`` columns: proteoform_key, cancer_code, prevalence_p95,
    n_samples. ``gene_keys`` maps genes to those exact expression identities.
    Supplying these inputs permits offline runs without substituting gene-level
    prevalence arithmetic. MS observations may be supplied or read live.
    """

    proteins: list[Protein]
    cta_gene_ids: set[str]
    gene_keys: dict[str, str]
    prevalence: pd.DataFrame
    burden: pd.DataFrame
    cohorts: dict[str, list[str]]
    ms_hits: pd.DataFrame | None = None
    provenance: dict = field(default_factory=dict)
    source_tables: dict[str, pd.DataFrame] = field(default_factory=dict)


def require_current_hitlist():
    """Check imported code, not possibly shadowing distribution metadata."""
    from hitlist.version import __version__
    from packaging.version import Version

    if Version(__version__) < Version("1.64.7"):
        raise ImportError(
            f"Vaccine design requires Hitlist >=1.64.7; imported {__version__}. Install tsarina[vaccine] in a compatible environment."
        )


def resolve_background_cta_ids(proteins, cta_ids):
    """Resolve same-gene alternate haplotypes without admitting unrelated loci.

    OncoRef curates primary gene IDs. Ensembl's HSCHR alternate haplotypes
    can annotate the same gene with another ID (PRAME, MAGEA3, MAGEA6).
    Those are CTA gene occurrences, not independent somatic sources. Only
    an exact symbol match to a curated non-HSCHR occurrence is inherited;
    ordinary non-CTA loci, even with an identical sequence, remain background.
    This resolution never changes expression keys or adds prevalence values.
    """
    names = {}
    for protein in proteins:
        if protein.gene_id in cta_ids and not protein.contig.startswith("HSCHR"):
            names.setdefault(protein.gene_name, set()).add(protein.gene_id)
    aliases = {}
    for protein in proteins:
        if (
            protein.gene_id not in cta_ids
            and protein.contig.startswith("HSCHR")
            and protein.gene_name in names
        ):
            aliases[protein.gene_id] = {
                "gene_id": protein.gene_id,
                "gene_name": protein.gene_name,
                "contig": protein.contig,
                "cta_reference_gene_ids": ";".join(sorted(names[protein.gene_name])),
            }
    table = pd.DataFrame(
        [aliases[key] for key in sorted(aliases)],
        columns=["gene_id", "gene_name", "contig", "cta_reference_gene_ids"],
    )
    return set(cta_ids) | aliases.keys(), table


def load_vaccine_inputs(
    definition="strict", ensembl_release=112, cohorts=None, auto_fetch=False, on_progress=None
):
    """Read OncoRef references and every translated Ensembl coding isoform.

    Ensembl assets must be installed/indexed (``pyensembl install --release``).
    Expression downloads are explicit via ``auto_fetch``. Never change the
    active Ensembl release or repair scientific reference data silently.
    """
    import oncoref
    from oncoref import cta, expression, proteoforms
    from pyensembl import EnsemblRelease

    require_current_hitlist()

    if definition not in {"strict", "loose"}:
        raise ValueError("CTA definition must be strict or loose")
    if definition == "loose" and not hasattr(cta, "cta_extended_gene_ids"):
        raise ImportError(
            "Loose CTAs require OncoRef >=1.8.207; install tsarina[vaccine] in a compatible environment"
        )
    ids = cta.cta_gene_ids() if definition == "strict" else cta.cta_extended_gene_ids()
    cohorts = DEFAULT_CANCER_COHORTS if cohorts is None else cohorts
    _validate_cohorts(cohorts)
    # The genome registry includes identical proteins omitted by the focused
    # CTA registry (e.g. RBMY1F/J); expression is collapsed BEFORE percentiling.
    keys = {g: proteoforms.proteoform_key(g, scope="genome") for g in ids}
    frames = []
    for code in dict.fromkeys(c for cs in cohorts.values() for c in cs):
        if on_progress:
            on_progress(f"Reading OncoRef p95 proteoform prevalence: {code}")
        frame = expression.proteoform_within_sample_top_fraction(
            code, threshold=0.95, scope="genome", auto_fetch=auto_fetch
        ).copy()
        frame["cancer_code"] = code
        frame = frame.rename(columns={"frac_samples_top5pct": "prevalence_p95"})
        # Keep the full denominator and provenance frame, including non-CTAs.
        frames.append(frame)
    if on_progress:
        on_progress(f"Loading translated Ensembl r{ensembl_release} proteins")
    proteins = []
    from .gene_sets import is_coding_gene, is_coding_transcript

    # Only sequences are needed here, not Hitlist's large inverted k-mer index.
    # Keep every translated transcript and its source occurrence, including
    # repeated ENSPs, so shared IDs cannot hide a non-CTA source gene.
    for gene in EnsemblRelease(ensembl_release).genes():
        if not is_coding_gene(gene):
            continue
        translated = [
            (t, t.protein_sequence)
            for t in gene.transcripts
            if is_coding_transcript(t) and t.protein_sequence
        ]
        if not translated:
            continue
        best = max(len(seq) for _, seq in translated)
        for transcript, seq in translated:
            proteins.append(
                Protein(
                    gene.id,
                    gene.name,
                    transcript.protein_id or transcript.id,
                    seq,
                    transcript.id,
                    len(seq) == best,
                    getattr(gene, "contig", ""),
                )
            )
    burden = oncoref.cancer_burden_df()
    tables = {
        "cta_reference": cta.cta_df() if definition == "strict" else cta.cta_extended_df(),
        "proteoform_reference": proteoforms.proteoform_groups(scope="genome"),
        "burden_reference": burden,
        "expression_availability": expression.cancer_reference_expression_availability(),
    }
    return VaccineInputs(
        proteins,
        ids,
        keys,
        pd.concat(frames, ignore_index=True),
        burden,
        cohorts,
        provenance={
            "definition": definition,
            "oncoref_version": version("oncoref"),
            "hitlist_version": version("hitlist"),
            "ensembl_release": ensembl_release,
            "expression_scope": "genome",
            "expression_threshold": 0.95,
            "protein_reference_sha256": sha256(
                "\n".join(
                    sorted(
                        f"{p.protein_id}:{p.transcript_id}:{p.gene_id}:{p.contig}:{p.sequence}"
                        for p in proteins
                    )
                ).encode()
            ).hexdigest(),
        },
        source_tables=tables,
    )


def _validate_cohorts(cohorts):
    if not cohorts or any(not cs or isinstance(cs, str) for cs in cohorts.values()):
        raise ValueError("Cancer cohort mapping must contain nonempty lists")
    flat = [c for cs in cohorts.values() for c in cs]
    if len(flat) != len(set(flat)):
        raise ValueError("A cohort cannot be counted twice across mortality categories")


def rank_proteoforms(inputs: VaccineInputs):
    """Collapse byte-identical canonical sequences; score each category once.

    Return ranking, complete per-cohort contributions, category summaries.
    Missing measurements yield an observed partial score with explicit
    missingness, rather than asserted zero expression. Registry disagreements fail.
    """
    _validate_cohorts(inputs.cohorts)
    if not inputs.cta_gene_ids:
        raise ValueError("CTA input set is empty")
    if set(inputs.cta_gene_ids) - set(inputs.gene_keys):
        raise ValueError("CTA genes are missing their proteoform expression keys")
    burden = inputs.burden.set_index("burden_category", verify_integrity=True)
    missing = set(inputs.cohorts) - set(burden.index)
    if missing:
        raise ValueError(f"Unknown mortality categories: {sorted(missing)}")
    weights = pd.to_numeric(burden.world_mortality_pct, errors="raise")
    # Published shares are rounded (OncoRef 1.8.207 sums to 100.2%). Preserve
    # their source values, without imposing spurious exact normalization.
    if weights.isna().any() or (weights < 0).any() or (weights > 100).any():
        raise ValueError("Invalid world mortality shares")
    best, occurrences = {}, {}
    for p in inputs.proteins:
        if p.gene_id not in inputs.cta_gene_ids:
            continue
        occurrences.setdefault(p.gene_id, []).append(p)
        previous = best.get(p.gene_id)
        if previous is None or (-len(p.sequence), p.protein_id) < (
            -len(previous.sequence),
            previous.protein_id,
        ):
            best[p.gene_id] = p
    grouped = {}
    for gene in sorted(inputs.cta_gene_ids):
        p = best.get(gene, Protein(gene, gene, "", ""))
        grouped.setdefault(p.sequence or f"missing:{gene}", []).append(p)
    prevalence = inputs.prevalence.copy()
    if prevalence.duplicated(["proteoform_key", "cancer_code"]).any():
        raise ValueError("Duplicate proteoform/cohort prevalence measurements")
    values = pd.to_numeric(prevalence.prevalence_p95, errors="raise")
    ns = pd.to_numeric(prevalence.n_samples, errors="raise")
    if ((values.dropna() < 0) | (values.dropna() > 1)).any() or (
        ns.isna().any() or (ns <= 0).any() or (ns % 1 != 0).any()
    ):
        raise ValueError("Prevalence must be a fraction with positive integer sample denominators")
    totals = {}
    for code, sub in prevalence.groupby("cancer_code"):
        counts = sub.n_samples.unique()
        if len(counts) != 1:
            raise ValueError(f"Inconsistent sample denominators for {code}")
        totals[code] = int(counts[0])
    lookup = prevalence.set_index(["proteoform_key", "cancer_code"])
    ranked, detailed, summaries = [], [], []
    for members in grouped.values():
        keys = {inputs.gene_keys[p.gene_id] for p in members}
        if len(keys) != 1:
            raise ValueError(
                "Identical protein sequences disagree with OncoRef expression grouping: "
                + "/".join(p.gene_name for p in members)
                + ". Update the upstream proteoform registry; prevalences cannot be added."
            )
        key = keys.pop()
        seq = members[0].sequence
        sources = [
            occurrence
            for member in members
            for occurrence in occurrences.get(member.gene_id, [member])
            if occurrence.sequence == seq
        ]
        name = "/".join(sorted({p.gene_name for p in members}))
        score, complete = 0.0, True
        for category, codes in inputs.cohorts.items():
            reference = burden.loc[category].to_dict()
            total = sum(totals.get(c, 0) for c in codes)
            weighted, measured = 0.0, 0
            for code in codes:
                loc = (key, code)
                value = lookup.loc[loc, "prevalence_p95"] if loc in lookup.index else float("nan")
                n = totals.get(code, 0)
                if pd.notna(value):
                    weighted += float(value) * n
                    measured += n
                detailed.append(
                    {
                        "proteoform_key": key,
                        "name": name,
                        "burden_category": category,
                        "cancer_code": code,
                        "prevalence_p95": value,
                        "n_samples": n or None,
                        "measurement_status": "observed" if pd.notna(value) else "missing",
                        **reference,
                    }
                )
            fully_measured = measured == total and total > 0 and all(c in totals for c in codes)
            # No denominator for an absent cohort: observed score is not claimed
            # to be a category-wide lower bound. Its missingness remains explicit.
            fraction = weighted / total if total else float("nan")
            contribution = (
                float(reference["world_mortality_pct"])
                / 100
                * (fraction if pd.notna(fraction) else 0)
            )
            score += contribution
            complete &= fully_measured
            summaries.append(
                {
                    "proteoform_key": key,
                    "name": name,
                    "burden_category": category,
                    "cancer_codes": ";".join(codes),
                    "prevalence_p95": fraction,
                    "n_samples": total,
                    "n_measured_samples": measured,
                    "complete_measurement": fully_measured,
                    "score_contribution": contribution,
                    **reference,
                }
            )
        ranked.append(
            {
                "proteoform_key": key,
                "name": name,
                "gene_ids": ";".join(p.gene_id for p in members),
                "protein_ids": ";".join(sorted({p.protein_id for p in sources})),
                "transcript_ids": ";".join(sorted({p.transcript_id for p in sources})),
                "sequence": seq,
                "sequence_sha256": sha256(seq.encode()).hexdigest(),
                "length_aa": len(seq),
                "mortality_weighted_score": score,
                "complete_measurement": complete,
            }
        )
    ranking = (
        pd.DataFrame(ranked)
        .sort_values(["mortality_weighted_score", "name"], ascending=[False, True])
        .reset_index(drop=True)
    )
    if ranking.proteoform_key.duplicated().any():
        raise ValueError(
            "OncoRef groups non-identical reference sequences; reconcile Ensembl release"
        )
    ranking.insert(0, "rank", range(1, len(ranking) + 1))
    return ranking, pd.DataFrame(detailed), pd.DataFrame(summaries)


def filter_ms_modality(hits):
    """Require positive MS modality, retaining non-MS/unknown records for audit.

    Hitlist #644: nonbinding does not imply MS. Curated supplement adapters
    publish MS-only rows; supplied observations may declare assay_modality when
    no structured method exists. Explicit non-MS methods always take precedence.
    """
    from .spanning import _is_truthy

    method = hits.get("assay_method", pd.Series("", index=hits.index)).fillna("").astype(str)
    method = method.str.strip().str.lower()
    modality = hits.get("assay_modality", pd.Series("", index=hits.index)).fillna("").astype(str)
    source = hits.get("source", pd.Series("", index=hits.index)).fillna("").astype(str)
    accepted = method.str.contains("mass spectrometry", regex=False) | (
        method.eq("") & (modality.eq("mass_spectrometry") | source.eq("supplement"))
    )
    if "is_binding_assay" in hits:
        accepted &= ~hits.is_binding_assay.map(_is_truthy)
    rejected = hits[~accepted].copy()
    rejected["rejection_reason"] = method[~accepted].map(
        lambda m: "explicit_non_ms_assay" if m else "unknown_assay_modality"
    )
    if "is_binding_assay" in rejected:
        rejected.loc[rejected.is_binding_assay.map(_is_truthy), "rejection_reason"] = (
            "binding_assay"
        )
    return hits[accepted].copy(), rejected


def sample_affinity_support(hits, scores, alleles, cutoff, allow_untyped):
    """Every binding panel allele in a typed sample, or permitted untyped MS.

    No best-of-haplotype or presentation-percentile gate is applied. The exact
    observed peptide remains the evidence unit; nested peptide inference is not
    used. Study-wide allele pools are not individual sample genotypes.
    """
    from .spanning import (
        _SAMPLE_NARROWED_PROVENANCES,
        _add_evidence,
        _exact_hla_alleles,
        _finalize_evidence_bucket,
        _is_truthy,
    )

    lookup = scores.set_index(["peptide", "allele"])
    stats, assignments = {}, []
    for _, row in hits.iterrows():
        restriction = _exact_hla_alleles(row.get("mhc_restriction", ""))
        if row.get("mhc_allele_provenance", "") in _SAMPLE_NARROWED_PROVENANCES:
            sample = _exact_hla_alleles(row.get("mhc_allele_set", ""))
        else:
            sample = restriction
        mono = _is_truthy(row.get("is_monoallelic", False)) and len(restriction) == 1
        if mono:
            sample = restriction
        tier = "monoallelic_ms" if mono else "sample_allele_ms"
        if not sample:
            if not allow_untyped:
                continue
            tier = "unrestricted_ms"
        for allele in alleles:
            if sample and allele not in sample:
                continue
            affinity = float(lookup.loc[(row.peptide, allele), "affinity_nm"])
            if affinity >= cutoff:
                continue
            _add_evidence(stats, (row.peptide, allele, tier), row)
            assignments.append(
                {
                    "peptide": row.peptide,
                    "allele": allele,
                    "evidence_tier": tier,
                    "sample_alleles": ";".join(sorted(sample)),
                    "affinity_nm": affinity,
                    **{k: row.get(k, "") for k in ("provenance_id", "assay_iri", "pmid", "source")},
                }
            )
    rows = []
    for peptide, allele in sorted({(p, a) for p, a, _ in stats}):
        for tier in ("monoallelic_ms", "sample_allele_ms", "unrestricted_ms"):
            if (peptide, allele, tier) in stats:
                score = lookup.loc[(peptide, allele)].to_dict()
                rows.append(
                    {
                        "peptide": peptide,
                        "allele": allele,
                        "length": len(peptide),
                        "evidence_tier": tier,
                        **_finalize_evidence_bucket(stats[peptide, allele, tier]),
                        **score,
                    }
                )
                break
    result = pd.DataFrame(rows)
    result.attrs["ms_assignments"] = pd.DataFrame(
        assignments,
        columns=[
            "peptide",
            "allele",
            "evidence_tier",
            "sample_alleles",
            "affinity_nm",
            "provenance_id",
            "assay_iri",
            "pmid",
            "source",
        ],
    )
    return result


def panel_ms_support(
    peptides,
    alleles,
    hits=None,
    predictor="mhcflurry",
    *,
    mode="presentation",
    affinity_nm=1000,
    allow_untyped=False,
):
    """All qualifying MS-supported pMHCs using the existing panel tier policy."""
    if mode not in {"presentation", "sample_affinity"}:
        raise ValueError("Invalid MS support mode")
    if not math.isfinite(affinity_nm) or affinity_nm <= 0:
        raise ValueError("MS affinity cutoff must be finite and positive")
    from .indexing import load_ms_evidence
    from .scoring import score_presentation
    from .spanning import (
        _build_evidence_stats,
        _candidate_rows,
        _score_alleles_for_panel,
        _score_lookup,
    )

    if not peptides:
        return pd.DataFrame()
    if hits is None:
        require_current_hitlist()
        hits = load_ms_evidence(peptides=set(peptides))
    if hits.empty:
        return pd.DataFrame()
    hits = hits[hits.peptide.isin(peptides)].copy()
    if hits.empty:
        return pd.DataFrame()
    for column, value in (("mhc_class", "I"), ("species", "Homo sapiens")):
        if column in hits and not hits[column].eq(value).all():
            raise ValueError("MS inputs must contain only human class-I observations")
    hits, rejected = filter_ms_modality(hits)
    attrs = {
        "queried_ms_observations": hits,
        "rejected_ms_observations": rejected,
        "ms_input_sha256": sha256(hits.to_csv(index=False).encode()).hexdigest(),
    }
    if hits.empty:
        result = pd.DataFrame()
        result.attrs.update(attrs)
        return result
    wanted = sorted(set(hits.peptide))
    score_alleles = _score_alleles_for_panel(alleles, hits)
    scores = score_presentation(wanted, score_alleles, predictor=predictor)
    if scores.duplicated(["peptide", "allele"]).any():
        raise ValueError("Duplicate MS predictions")
    actual = set(zip(scores.peptide, scores.allele))
    if any((p, a) not in actual for p in wanted for a in score_alleles):
        raise ValueError(
            "Incomplete MS presentation prediction coverage; allele inference cannot proceed"
        )
    if mode == "sample_affinity":
        affinities = pd.to_numeric(scores.affinity_nm, errors="raise")
        if not affinities.gt(0).all() or not affinities.map(lambda a: a < float("inf")).all():
            raise ValueError("Invalid MS affinity predictions; allele inference cannot proceed")
        result = sample_affinity_support(hits, scores, alleles, affinity_nm, allow_untyped)
        result.attrs.update(attrs)
        return result
    percentiles = pd.to_numeric(scores.presentation_percentile, errors="raise")
    if not percentiles.between(0, 100).all():
        raise ValueError("Invalid MS presentation percentiles; allele inference cannot proceed")
    meta = pd.DataFrame({"peptide": wanted, "gene_name": "candidate"})
    meta["length"] = meta.peptide.str.len()
    evidence = _build_evidence_stats(hits, alleles, _score_lookup(scores))
    result = _candidate_rows(
        scores,
        meta,
        alleles,
        evidence,
        {"monoallelic_ms": 2.0, "sample_allele_ms": 1.0, "unrestricted_ms": 0.5},
        False,
    )
    result = result.drop(columns=[c for c in result if c.startswith("_")]).reset_index(drop=True)
    result.attrs.update(attrs)
    return result
