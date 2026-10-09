"""Expression-union bounds and HLA carrier proxies, never clinical coverage."""

from __future__ import annotations

import math

import pandas as pd

TIERS = ("monoallelic_ms", "sample_allele_ms", "unrestricted_ms")


def paired_target_coverage(targets, pairs, supported_alleles, missing_mass=0.0):
    """Joint expression/MHC evidence bounds on explicit tumor/genotype pairs.

    Each target supplies ``expression[donor] = (lower, upper)`` and a set of
    qualified ``alleles``. Weights describe the declared sampling frame, not
    universal population frequencies. Full genotypes can include duplicated
    loci. Unsupported alleles and missing individuals remain unknown mass.
    This shared estimator is independent of human HWE/carrier assumptions.
    """
    if not math.isfinite(missing_mass) or not 0 <= missing_mass <= 1:
        raise ValueError("Invalid missing population mass")
    if len({p["pair_id"] for p in pairs}) != len(pairs):
        raise ValueError("Duplicate tumor/genotype pair")
    if any(not math.isfinite(p["weight"]) or p["weight"] <= 0 for p in pairs):
        raise ValueError("Positive finite pair weights required")
    if any(
        p["pairing"] not in {"observed", "simulated_independent_cohorts"}
        or (p["pairing"] == "observed" and p["tumor_donor"] != p["genotype_donor"])
        for p in pairs
    ):
        raise ValueError("Explicit observed matching dogs or simulated pairing required")
    if not pairs and missing_mass != 1:
        raise ValueError("Empty pairs require fully missing population mass")
    observed = [p["tumor_donor"] for p in pairs if p["pairing"] == "observed"]
    if len(set(observed)) != len(observed):
        raise ValueError("Repeated observed dogs cannot inflate coverage")
    total = sum(p["weight"] for p in pairs)
    scale = (1 - missing_mass) / total if total else 0.0
    rows = []
    for pair in pairs:
        genotype = set(pair["alleles"])
        unsupported = genotype - set(supported_alleles)
        unknown_genotype = not genotype or bool(unsupported)
        lower, upper, rna_lower, rna_upper = set(), set(), set(), set()
        for key, target in targets.items():
            lo, hi = target["expression"].get(pair["tumor_donor"], (0, 1))
            if not 0 <= lo <= hi <= 1 or lo not in {0, 1} or hi not in {0, 1}:
                raise ValueError("Pair expression requires binary lower/upper evidence bounds")
            match = bool(genotype & set(target["alleles"]))
            if lo:
                rna_lower.add(key)
            if hi:
                rna_upper.add(key)
            if lo and match:
                lower.add(key)
            if hi and (match or unknown_genotype):
                upper.add(key)
        rows.append(
            {
                "pair_id": pair["pair_id"],
                "tumor_donor": pair["tumor_donor"],
                "genotype_donor": pair["genotype_donor"],
                "pairing": pair["pairing"],
                "weight": pair["weight"] * scale,
                "unsupported_alleles": ";".join(sorted(unsupported)),
                "genotype_unassessed": unknown_genotype,
                "qualified_targets_lower": len(lower),
                "qualified_targets_upper": len(upper),
                "expression_targets_lower": len(rna_lower),
                "expression_targets_upper": len(rna_upper),
            }
        )
    result = {
        "n_pairs": len(pairs),
        "missing_mass": missing_mass,
        "unsupported_genotype_mass": sum(r["weight"] for r in rows if r["genotype_unassessed"]),
        "pairing_modes": sorted({p["pairing"] for p in pairs}),
        "pairs": rows,
    }
    for label, field in (("joint", "qualified"), ("expression", "expression")):
        for count in (1, 2):
            for bound in ("lower", "upper"):
                value = sum(r["weight"] for r in rows if r[f"{field}_targets_{bound}"] >= count)
                if bound == "upper" and len(targets) >= count:
                    value += missing_mass
                result[f"{label}_{count}_{bound}"] = value
    return result


def carrier_reach(alleles, frequencies):
    """HWE within locus, linkage equilibrium across loci; missing is explicit.

    Frequencies are allele probabilities, not carrier frequencies. They are
    never renormalized to the queried panel. The known-allele estimate excludes
    missing values and is a lower estimate within this model, not a confidence
    bound on a real population. See Bui et al., BMC Bioinformatics 2006;7:153.
    """
    sums = dict.fromkeys("ABC", 0.0)
    missing = []
    for allele in sorted(set(alleles)):
        if not allele.startswith(("HLA-A*", "HLA-B*", "HLA-C*")):
            raise ValueError(f"Unsupported coverage allele: {allele}")
        value = frequencies.get(allele)
        if value is None or pd.isna(value):
            missing.append(allele)
            continue
        if not math.isfinite(float(value)) or not 0 <= float(value) <= 1:
            raise ValueError(f"Invalid allele frequency for {allele}")
        sums[allele[4]] += float(value)
    if any(value > 1 + 1e-12 for value in sums.values()):
        raise ValueError("Covered allele frequencies exceed one within an HLA locus")
    loci = {locus: 1 - (1 - value) ** 2 for locus, value in sums.items()}
    return {
        "reach": 1 - math.prod(1 - p for p in loci.values()),
        "loci": loci,
        "frequency_sums": sums,
        "missing_alleles": missing,
    }


def union_bounds(values):
    """Fréchet bounds from marginal prevalences; unknown marginals widen upper."""
    values = list(values)
    known = [float(v) for v in values if v is not None and not pd.isna(v)]
    if any(not math.isfinite(v) or not 0 <= v <= 1 for v in known):
        raise ValueError("Expression prevalence must be finite and in [0,1]")
    missing = len(values) - len(known)
    return max(known, default=0.0), 1.0 if missing else min(1.0, sum(known)), missing


def coverage_tables(proteins, ligands, cancer, layers, frequencies):
    """Cumulative rank-ordered protein and final-order segment evidence.

    Cancer-weighted bounds use *global* burden shares (not a renormalized
    subset); represented mortality/incidence percentages are separate columns.
    Segment curves describe prefixes of the final construct, not independently
    optimized shorter vaccines. Deduplication takes strongest restriction tier.
    """
    retained = ligands[ligands.assembled.eq(True)].copy()
    summaries, cancer_rows = [], []

    def append(axis, step, name, keys, hits, aa):
        tiers = {tier: set() for tier in TIERS}
        for row in hits.itertuples(index=False):
            tiers[row.evidence_tier].add(row.allele)
        row = {"axis": axis, "step": step, "name": name, "length_aa": aa}
        for label, chosen in (
            ("all", set().union(*tiers.values())),
            ("typed", tiers[TIERS[0]] | tiers[TIERS[1]]),
            ("measured", tiers[TIERS[0]]),
        ):
            estimate = carrier_reach(chosen, frequencies)
            row[f"hla_{label}"] = estimate["reach"]
            row[f"hla_{label}_missing"] = ";".join(estimate["missing_alleles"])
        row["alleles"] = hits.allele.nunique()
        row["pmhc"] = len(hits.drop_duplicates(["peptide", "allele"]))
        strongest = hits.assign(tier_order=hits.evidence_tier.map(dict(zip(TIERS, range(3)))))
        strongest = strongest.sort_values("tier_order").drop_duplicates("peptide")
        row["ms_peptides"] = len(strongest)
        for tier in TIERS:
            row[f"peptides_{tier}"] = int(strongest.evidence_tier.eq(tier).sum())
        for metric in ("mortality", "incidence"):
            row[f"{metric}_lower"] = row[f"{metric}_upper"] = 0.0
            row[f"represented_{metric}_pct"] = 0.0
        for category, group in cancer.groupby("burden_category", sort=False):
            selected = group[group.proteoform_key.isin(keys)]
            # Missing entire rows and partial cohort measurements are unknown.
            known = selected.set_index("proteoform_key").prevalence_p95.to_dict()
            if "complete_measurement" in selected:
                for item in selected.itertuples(index=False):
                    if not item.complete_measurement:
                        known[item.proteoform_key] = None
            lower, upper, missing = union_bounds(known.get(key) for key in keys)
            detail = {
                "axis": axis,
                "step": step,
                "burden_category": category,
                "lower": lower,
                "upper": upper,
                "missing_proteins": missing,
            }
            for metric in ("mortality", "incidence"):
                weight = float(group.iloc[0][f"world_{metric}_pct"])
                if not math.isfinite(weight) or not 0 <= weight <= 100:
                    raise ValueError("Invalid cancer burden share")
                row[f"{metric}_lower"] += weight / 100 * lower
                row[f"{metric}_upper"] += weight / 100 * upper
                row[f"represented_{metric}_pct"] += weight
                detail[f"world_{metric}_pct"] = weight
            cancer_rows.append(detail)
        summaries.append(row)

    keys, aa = set(), 0
    append("protein", 0, "No proteins", keys, retained.iloc[:0], aa)
    for i, p in enumerate(proteins.sort_values("rank").itertuples(index=False), 1):
        keys.add(p.proteoform_key)
        aa += p.assembled_aa
        append("protein", i, p.name, keys, retained[retained.proteoform_key.isin(keys)], aa)
    keys, hits = set(), retained.iloc[:0]
    append("segment", 0, "No segments", keys, hits, 0)
    native = [layer for layer in layers if layer["kind"] == "cta_segment"]
    for i, layer in enumerate(native, 1):
        keys.add(layer["proteoform_key"])
        hits = pd.concat([hits, retained[retained.segment_id.eq(layer["segment_id"])]])
        append("segment", i, layer["name"], keys, hits, layer["end_aa"])
    return pd.DataFrame(summaries), pd.DataFrame(cancer_rows)


def ciwd_frequencies(audit):
    """One coherent CIWD table; alternate allotype sources remain in the audit."""
    return {
        row.allele: float(row.published_global_frequency)
        for row in audit.itertuples(index=False)
        if str(row.global_source_label).startswith("CIWD")
        and pd.notna(row.published_global_frequency)
    }
