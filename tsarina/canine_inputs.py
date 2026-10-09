"""Frozen canine evidence import; no human expression or mortality defaults.

This imports reviewed sequence-group abundance and restriction decisions. It
does not discover CTAs or allocate ambiguous RNA reads to transcripts.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from hashlib import sha256
from pathlib import Path

import pandas as pd

from .mhc import parse_mhc
from .peptides import AA20
from .vaccine_inputs import Protein, VaccineInputs


def canonical_hash(value):
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def sequence_hash(sequence):
    if not isinstance(sequence, str) or not sequence or set(sequence) - AA20:
        raise ValueError("Canine proteins/peptides require complete canonical AA sequences")
    return sha256(sequence.encode("ascii")).hexdigest()


def dla_alleles(names):
    """Canonical class-I DLA names; preserve full multi-copy genotypes."""
    result = []
    for name in names:
        parsed = parse_mhc(name, expect="allele", species="DLA")
        if (
            parsed is None
            or not str(parsed.mhc_class).startswith("I")
            or str(parsed.mhc_class).startswith("II")
        ):
            raise ValueError(f"Unresolved or non-class-I DLA allele: {name}")
        if len(parsed.allele_fields) < 2:
            raise ValueError(f"Unresolved DLA allele: {name}")
        result.append(parsed.to_string())
    return list(dict.fromkeys(result))


def _required(row, names):
    if any(not isinstance(row.get(n), str) or not row[n].strip() for n in names):
        raise ValueError(f"Nonempty provenance/identity fields required: {', '.join(names)}")


def _unique(rows, key):
    values = [r[key] for r in rows]
    if len(set(values)) != len(values):
        raise ValueError(f"Duplicate {key}")
    return dict(zip(values, rows))


def _sha(value):
    return isinstance(value, str) and len(value) == 64 and not set(value) - set("0123456789abcdef")


def validate_canine_bundle(bundle):
    """Validate the complete reference/evidence contract before ranking."""
    if bundle.get("schema") != "tsarina.canine.v1" or bundle.get("taxon") != 9615:
        raise ValueError("Expected tsarina.canine.v1 bundle with taxon 9615")
    reference = bundle["reference"]
    _required(reference, ["assembly_accession", "annotation_release", "source_version"])
    if reference.get("taxon") != 9615:
        raise ValueError("Mixed reference taxon")
    assets = reference["asset_hashes"]
    if not {"assembly", "annotation"} <= assets.keys() or any(not _sha(v) for v in assets.values()):
        raise ValueError("Reference requires assembly/annotation SHA256s")
    ref_key = canonical_hash(reference)
    if bundle.get("reference_key") != ref_key:
        raise ValueError("Incompatible reference key")
    if bundle.get("background_complete") is not True:
        raise ValueError(
            "Declare a complete translated reference background; incomplete panels cannot design"
        )
    sources = bundle["sources"]
    if not sources:
        raise ValueError("Source inventory is required")
    for source in sources.values():
        _required(source, ["url", "version", "license"])
        if not source["url"].startswith("https://"):
            raise ValueError("Source URLs must use https")
        if not _sha(source.get("sha256")):
            raise ValueError("Each source requires its asset SHA256")
    policy = bundle["restriction_policy"]
    _required(policy, ["name", "version", "normal_assessment_source"])
    if policy["normal_assessment_source"] not in sources:
        raise ValueError("Normal assessment requires an inventoried source")
    normal_threshold = policy["somatic_tpm_threshold"]
    if (
        not isinstance(normal_threshold, (int, float))
        or not math.isfinite(normal_threshold)
        or normal_threshold <= 0
    ):
        raise ValueError("Explicit finite positive somatic TPM threshold required")
    allowed = policy["allowed_tissues"]
    if not allowed or any(not isinstance(t, str) or not t.strip() for t in allowed):
        raise ValueError("Declare the allowed canine normal tissues")
    if policy["definition"] not in {"strict", "loose"}:
        raise ValueError("Restriction definition must be strict or loose")
    threshold = bundle["expression_policy"]["threshold"]
    if not isinstance(threshold, (int, float)) or not math.isfinite(threshold) or threshold <= 0:
        raise ValueError("Finite positive absolute expression threshold required")
    if bundle["expression_policy"]["unit"] != "TPM":
        raise ValueError("Canine absolute prevalence currently requires TPM bounds")
    cohorts = bundle["cohorts"]
    for cohort in cohorts.values():
        _required(cohort, ["histology", "description", "population", "quantification"])
    if not cohorts:
        raise ValueError("At least one explicit canine cohort is required")
    occurrences = _unique(bundle["occurrences"], "occurrence_id")
    samples = _unique(bundle["samples"], "sample_id")
    seq_ids = set()
    for kind in ("occurrences", "samples", "rna_bounds", "ms_hits"):
        for row in bundle[kind]:
            if row.get("reference_key") != ref_key or row.get("taxon") != 9615:
                raise ValueError(f"Mixed reference or taxon in {kind}")
            if row.get("source_id") not in sources:
                raise ValueError(f"Unknown source in {kind}")
    for row in occurrences.values():
        _required(row, ["gene_id", "transcript_id", "protein_id", "restriction_reason"])
        if row.get("complete") is not True:
            raise ValueError("Incomplete source proteins must be quarantined before import")
        if row["sequence_id"] != sequence_hash(row["sequence"]):
            raise ValueError("Incorrect exact protein sequence identity")
        seq_ids.add(row["sequence_id"])
        if row["restriction_status"] not in {"admitted", "rejected", "unknown"}:
            raise ValueError("Unknown restriction status")
        if row["restriction_status"] == "admitted" and row.get("normal_assessment") != "assessed":
            raise ValueError("Admission requires assessed normal evidence")
    for sample in samples.values():
        _required(
            sample,
            [
                "study",
                "specimen",
                "preparation",
                "health",
                "histology",
                "breed",
                "treatment",
                "pooling",
                "qc",
                "quantification",
                "unit",
            ],
        )
        if sample["unit"] != "TPM":
            raise ValueError("RNA sample unit mismatch; do not relabel counts as TPM")
        if sample.get("donor") is not None and (
            not isinstance(sample["donor"], str) or not sample["donor"].strip()
        ):
            raise ValueError("Donor must be a reconciled identity or null")
        if sample.get("cohort") and sample["cohort"] not in cohorts:
            raise ValueError("Unknown cohort")
        if sample.get("cohort") and sample["histology"] != cohorts[sample["cohort"]]["histology"]:
            raise ValueError("Cohort histology mismatch")
        if (
            sample.get("cohort")
            and sample["quantification"] != cohorts[sample["cohort"]]["quantification"]
        ):
            raise ValueError("Cohort RNA quantification mismatch")
    seen = set()
    for row in bundle["rna_bounds"]:
        key = (row["sample_id"], row["sequence_id"])
        if key in seen or key[0] not in samples or key[1] not in seq_ids:
            raise ValueError("Duplicate/unresolved sequence-group RNA bounds")
        seen.add(key)
        _required(row, ["allocation_provenance"])
        lo, hi = row["lower"], row["upper"]
        if (lo is None) != (hi is None) or (
            lo is not None
            and (
                not all(isinstance(v, (int, float)) and math.isfinite(v) for v in (lo, hi))
                or not 0 <= lo <= hi
            )
        ):
            raise ValueError("Invalid RNA bounds")
    _unique(bundle["ms_hits"], "observation_id")
    for row in bundle["ms_hits"]:
        sequence_hash(row["peptide"])
        _required(
            row,
            [
                "assay_method",
                "qualitative_measurement",
                "assay_context",
                "mhc_class",
                "host_species",
                "sample_kind",
                "sample_health",
                "sample_tissue",
                "restriction_kind",
                "sample_id",
                "ms_admission_reason",
                "assay_response",
            ],
        )
        if type(row.get("is_ms_observation")) is not bool:
            raise ValueError("Explicit reviewed positive-MS source decision required")
        if row.get("host_taxon") not in {9615, 9606}:
            raise ValueError("Declare the independent dog/human host taxon")
        if row.get("mhc_taxon") != 9615 or row.get("peptide_source_taxon") not in {9615, 9606}:
            raise ValueError("Declare independent peptide-source and DLA taxon facets")
        if row.get("verified_normal") is True:
            _required(row, ["normal_verification_source"])
        if row["restriction_kind"] not in {"monoallelic", "sample_genotype", "untyped"}:
            raise ValueError("MS restriction must be monoallelic, sample_genotype or untyped")
        alleles = dla_alleles(row["sample_alleles"])
        if row["restriction_kind"] == "monoallelic" and len(alleles) != 1:
            raise ValueError("Monoallelic MS requires one exact DLA allele")
        if row["restriction_kind"] == "sample_genotype" and not alleles:
            raise ValueError("Typed sample MS requires a genotype")
        if row["restriction_kind"] == "untyped" and alleles:
            raise ValueError("Untyped MS cannot also assert a sample genotype")
        if row["host_taxon"] != 9615 and not (
            row["restriction_kind"] == "monoallelic"
            and row["sample_kind"] == "transfected_cell_line"
        ):
            raise ValueError("Non-dog host evidence requires explicit DLA monoallelic transfection")
        if row["host_taxon"] == 9615 and row["peptide_source_taxon"] != 9615:
            raise ValueError("Endogenous dog MS requires canine peptide-source attribution")
    capability_names = []
    for cap in bundle["capabilities"]:
        allele = dla_alleles([cap["allele"]])[0]
        capability_names.append(allele)
        if cap["tier"] not in {"unsupported", "sequence_extrapolated", "empirically_supported"}:
            raise ValueError("Unknown DLA capability tier")
        _required(cap, ["reason"])
        if cap["tier"] != "unsupported":
            _required(cap, ["model", "version", "validation_source"])
            if not _sha(cap.get("weights_sha256")) or not _sha(cap.get("sequence_sha256")):
                raise ValueError("Executable DLA capability needs model/sequence hashes")
            if not cap["lengths"] or any(
                type(k) is not int or not 8 <= k <= 15 for k in cap["lengths"]
            ):
                raise ValueError("DLA capability requires explicit 8..15 lengths")
            if cap["tier"] == "empirically_supported" and (
                cap.get("validation_taxon") != 9615 or not _sha(cap.get("validation_sha256"))
            ):
                raise ValueError(
                    "Empirical canine capability requires a species-scoped validation artifact hash"
                )
    if len(set(capability_names)) != len(capability_names):
        raise ValueError("Duplicate canonical DLA capabilities")
    prediction_keys = set()
    if bundle.get("predictions") and bundle.get("prediction_capability_sha256") != canonical_hash(
        bundle["capabilities"]
    ):
        raise ValueError("Affinity archive capability fingerprint mismatch")
    for row in bundle.get("predictions", []):
        sequence_hash(row["peptide"])
        allele = dla_alleles([row["allele"]])[0]
        key = (row["peptide"], allele)
        if key in prediction_keys or allele not in capability_names:
            raise ValueError("Duplicate/unresolved affinity predictions")
        prediction_keys.add(key)
        if (
            not isinstance(row["affinity_nm"], (int, float))
            or not math.isfinite(row["affinity_nm"])
            or row["affinity_nm"] <= 0
        ):
            raise ValueError("Invalid affinity nM")
    population = bundle["genotype_population"]
    _required(population, ["name", "description"])
    if population["cohort"] not in cohorts:
        raise ValueError("Genotype population needs an explicit named cohort")
    mass = population["missing_mass"]
    if not isinstance(mass, (int, float)) or not math.isfinite(mass) or not 0 <= mass <= 1:
        raise ValueError("Invalid missing genotype mass")
    _unique(bundle["tumor_genotype_pairs"], "pair_id")
    observed = set()
    for pair in bundle["tumor_genotype_pairs"]:
        _required(
            pair, ["pair_id", "tumor_donor", "genotype_donor", "cohort", "breed", "source_id"]
        )
        if pair["cohort"] not in cohorts or pair["source_id"] not in sources:
            raise ValueError("Unknown pair cohort/source")
        if pair["pairing"] not in {"observed", "simulated_independent_cohorts"}:
            raise ValueError("Declare observed or simulated genotype pairing")
        if pair["pairing"] == "observed":
            key = (pair["cohort"], pair["tumor_donor"])
            if pair["tumor_donor"] != pair["genotype_donor"] or key in observed:
                raise ValueError("Observed pairs require unique independent matching dogs")
            observed.add(key)
        dla_alleles(pair["alleles"])
        if (
            not isinstance(pair["weight"], (int, float))
            or not math.isfinite(pair["weight"])
            or pair["weight"] <= 0
        ):
            raise ValueError("Positive finite genotype pair weight required")
    return bundle


def load_canine_vaccine_inputs(path):
    """Read a versioned offline bundle into the existing portable input route."""
    raw = Path(path).read_bytes()
    bundle = validate_canine_bundle(json.loads(raw))
    proteins = [
        Protein(
            r["gene_id"],
            r.get("gene_symbol") or r["gene_id"],
            r["protein_id"],
            r["sequence"],
            r["transcript_id"],
        )
        for r in bundle["occurrences"]
    ]
    return VaccineInputs(
        proteins,
        set(),
        {},
        pd.DataFrame(),
        pd.DataFrame(),
        {},
        provenance={
            "taxon": 9615,
            "bundle_sha256": sha256(raw).hexdigest(),
            "bundle_canonical_sha256": canonical_hash(bundle),
            "reference_key": bundle["reference_key"],
        },
        canine_evidence=bundle,
    )


def canine_prevalence(bundle, cohort):
    """Independent-dog bounds; repeated biopsies cannot inflate prevalence."""
    validate_canine_bundle(bundle)
    if cohort not in bundle["cohorts"]:
        raise ValueError(f"Unknown canine cohort: {cohort}")
    threshold = bundle["expression_policy"]["threshold"]
    normal_threshold = bundle["restriction_policy"]["somatic_tpm_threshold"]
    allowed = set(bundle["restriction_policy"]["allowed_tissues"])
    samples = {s["sample_id"]: s for s in bundle["samples"]}
    rna = {(r["sample_id"], r["sequence_id"]): r for r in bundle["rna_bounds"]}
    groups = defaultdict(list)
    for row in bundle["occurrences"]:
        groups[row["sequence_id"]].append(row)
    dogs, sample_audit = defaultdict(list), []
    for sample in bundle["samples"]:
        if sample.get("cohort") != cohort:
            continue
        reason = "included"
        if (
            sample["health"] != "tumor"
            or sample["specimen"] != "primary_tumor"
            or sample["preparation"] != "bulk"
            or sample["treatment"] != "untreated"
        ):
            reason = "not_untreated_primary_bulk_tumor"
        elif sample["qc"] != "pass" or sample["pooling"] != "none":
            reason = "qc_or_pooling_unresolved"
        elif not sample["donor"]:
            reason = "unresolved_donor"
        sample_audit.append({**sample, "inclusion_reason": reason})
        if reason == "included":
            dogs[sample["donor"]].append(sample)
    if not dogs:
        raise ValueError("No independent untreated primary dogs in the named cohort")
    rankings, per_dog = [], []
    for sid, occurrences in sorted(groups.items()):
        statuses = {r["restriction_status"] for r in occurrences}
        reasons = sorted({r["restriction_reason"] for r in occurrences})
        status = (
            "admitted"
            if statuses == {"admitted"}
            else "rejected"
            if "rejected" in statuses
            else "unknown"
        )
        # A positive or unresolved healthy-primary somatic bound prevents a
        # whole-protein reproductive-restriction claim, even for identical loci.
        for (sample_id, sequence_id), bound in rna.items():
            sample = samples[sample_id]
            if (
                sequence_id == sid
                and sample["health"] == "healthy"
                and sample["specimen"] == "primary_tissue"
                and sample.get("tissue") not in allowed
            ):
                if bound["lower"] is not None and bound["lower"] >= normal_threshold:
                    status = "rejected"
                    reasons.append(f"healthy_primary_RNA:{sample_id}")
                elif bound["upper"] is None or bound["upper"] >= normal_threshold:
                    if status != "rejected":
                        status = "unknown"
                    reasons.append(f"unresolved_normal_RNA:{sample_id}")
        confirmed = possible = evaluable = 0
        for donor, biopsies in sorted(dogs.items()):
            values = [rna.get((s["sample_id"], sid)) for s in biopsies]
            complete = all(v is not None and v["lower"] is not None for v in values)
            known = [v for v in values if v is not None and v["lower"] is not None]
            lower = bool(complete and all(v["lower"] >= threshold for v in known))
            upper = bool(not complete or any(v["upper"] >= threshold for v in known))
            confirmed += lower
            possible += upper
            evaluable += complete
            per_dog.append(
                {
                    "proteoform_key": sid,
                    "donor": donor,
                    "expression_lower": int(lower),
                    "expression_upper": int(upper),
                    "complete": complete,
                    "n_biopsies": len(biopsies),
                    "sample_ids": ";".join(s["sample_id"] for s in biopsies),
                }
            )
        rankings.append(
            {
                "proteoform_key": sid,
                "name": "/".join(
                    sorted({r.get("gene_symbol") or r["gene_id"] for r in occurrences})
                ),
                "gene_ids": ";".join(sorted({r["gene_id"] for r in occurrences})),
                "protein_ids": ";".join(sorted({r["protein_id"] for r in occurrences})),
                "transcript_ids": ";".join(sorted({r["transcript_id"] for r in occurrences})),
                "occurrence_ids": ";".join(sorted(r["occurrence_id"] for r in occurrences)),
                "sequence": occurrences[0]["sequence"],
                "length_aa": len(occurrences[0]["sequence"]),
                "restriction_status": status,
                "restriction_reason": ";".join(sorted(set(reasons))),
                "cohort": cohort,
                "positive_dogs": confirmed,
                "possible_positive_dogs": possible,
                "evaluable_dogs": evaluable,
                "total_identified_dogs": len(dogs),
                "prevalence_lower": confirmed / len(dogs),
                "prevalence_upper": possible / len(dogs),
                "threshold": threshold,
                "unit": "TPM",
            }
        )
    frame = (
        pd.DataFrame(rankings)
        .sort_values(["prevalence_lower", "proteoform_key"], ascending=[False, True])
        .reset_index(drop=True)
    )
    frame.insert(0, "rank", range(1, len(frame) + 1))
    return frame, pd.DataFrame(per_dog), pd.DataFrame(sample_audit)
