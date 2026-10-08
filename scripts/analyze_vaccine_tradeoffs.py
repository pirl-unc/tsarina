"""Allocation-only sensitivity analysis, using saved positive-MS assignments.

Run with a saved vaccine website's data.json and downloads directory. This does
not assemble or validate a new vaccine; junction constraints are not optimized.
Requires scipy>=1.9 in addition to Tsarina's vaccine reporting dependencies.
"""

import argparse
import json
import sys
from hashlib import sha256
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

from tsarina.vaccine_coverage import carrier_reach, ciwd_frequencies
from tsarina.vaccine_sequences import peptide_occurrences
from tsarina.version import __version__

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input-root", type=Path, default=ROOT / "docs/vaccine-results/downloads")
parser.add_argument("--report-json", type=Path, default=ROOT / "docs/vaccine-results/data.json")
parser.add_argument("--output-dir", type=Path, required=True)
parser.add_argument(
    "--tolerances",
    type=float,
    nargs="+",
    default=[0.0, 0.05, 0.1, 0.25],
    help="Maximum loss in percentage points for each modeled cancer-expression bound",
)
parser.add_argument(
    "--global-union-only",
    action="store_true",
    help="Relax target-specific and restriction-tier constraints for comparison",
)
args = parser.parse_args()
if any(t < 0 for t in args.tolerances):
    parser.error("Coverage-loss tolerances must be nonnegative")
DATA = json.loads(args.report_json.read_text())
OUTPUT = args.output_dir
OUTPUT.mkdir(parents=True, exist_ok=True)


def analyze_mode(mode):
    source = args.input_root / mode

    def table_path(name):
        plain = source / name
        return plain if plain.exists() else source / (name + ".gz")

    rank = pd.read_csv(table_path("ranking.csv")).set_index("proteoform_key")
    intervals = pd.read_csv(table_path("specific_intervals.csv"))
    assignments = pd.read_csv(table_path("ms_assignments.csv")).drop_duplicates(
        ["peptide", "allele", "evidence_tier"]
    )
    cancer = pd.read_csv(table_path("cancer_summary.csv"))
    cancer = cancer[cancer.complete_measurement.eq(True)]
    ligands = pd.read_csv(table_path("ligands.csv"))
    retained = ligands[ligands.assembled.eq(True)]
    frequencies = ciwd_frequencies(pd.DataFrame(DATA["designs"][mode]["hla"]))
    required_alleles = set(retained.allele)
    published = DATA["designs"][mode]
    overhead_nt = published["total_nt"] - 3 * published["protein_aa"]
    aa_limit = min(published["max_aa"], (published["max_nt"] - overhead_nt) // 3)
    segments = []
    for interval in intervals.itertuples(index=False):
        hits = []
        for peptide, group in assignments.groupby("peptide", sort=False):
            for a, b in peptide_occurrences(interval.sequence, peptide):
                for h in group.itertuples(index=False):
                    hits.append(
                        (interval.start + a, interval.start + b, peptide, h.allele, h.evidence_tier)
                    )
        if not hits:
            continue
        a, b = min(h[0] for h in hits), max(h[1] for h in hits)
        segments.append(
            {
                "segment_id": f"{interval.proteoform_key}:{interval.start}-{interval.end}",
                "key": interval.proteoform_key,
                "name": interval.name,
                "start": a,
                "end": b,
                "length_aa": b - a,
                "peptides": {h[2] for h in hits},
                "alleles": {h[3] for h in hits},
                "hits": hits,
            }
        )
    keys = sorted({s["key"] for s in segments})
    peptides = sorted(set().union(*(s["peptides"] for s in segments)))
    ns, nk, np_ = len(segments), len(keys), len(peptides)
    print(mode, "candidate segments", ns, "proteoforms", nk, "peptides", np_, flush=True)
    key_index = {k: ns + i for i, k in enumerate(keys)}
    pep_index = {p: ns + nk + i for i, p in enumerate(peptides)}
    rows, lows, highs = [], [], []

    def add(row, low=-np.inf, high=np.inf):
        rows.append(row)
        lows.append(low)
        highs.append(high)

    # Reserve one initiating methionine; no terminal padding or linkers are credited.
    add({i: s["length_aa"] for i, s in enumerate(segments)}, high=aa_limit - 1)
    for key in keys:
        members = [i for i, s in enumerate(segments) if s["key"] == key]
        for i in members:
            add({i: 1, key_index[key]: -1}, high=0)
        add({key_index[key]: 1, **dict.fromkeys(members, -1)}, high=0)
    for peptide in peptides:
        add(
            {
                pep_index[peptide]: 1,
                **{i: -1 for i, s in enumerate(segments) if peptide in s["peptides"]},
            },
            high=0,
        )
    for allele in sorted(required_alleles):
        add({i: 1 for i, s in enumerate(segments) if allele in s["alleles"]}, low=1)

    # Preserve the strongest current restriction tier for every global allele,
    # and every common-target allele. Extra observed peptides cannot conceal
    # weaker presenters for PRAME, MAGEA4 or XAGE1A/B.
    tier_order = {"monoallelic_ms": 0, "sample_allele_ms": 1, "unrestricted_ms": 2}
    restrictions = []
    for target, frame in [
        ("global", retained),
        *[
            (name, retained[retained.name.eq(name)])
            for name in ("PRAME", "MAGEA4", "XAGE1A/XAGE1B")
        ],
    ]:
        if args.global_union_only:
            continue
        for allele, group in frame.groupby("allele"):
            required_tier = min(tier_order[t] for t in group.evidence_tier)
            options = {
                i: 1
                for i, s in enumerate(segments)
                if (target == "global" or s["name"] == target)
                and any(h[3] == allele and tier_order[h[4]] <= required_tier for h in s["hits"])
            }
            add(options, low=1)
            restrictions.append({"target": target, "allele": allele, "tier": required_tier})

    coverage_rows = {"mortality": {}, "incidence": {}}
    nvars = ns + nk + np_
    # Indicator levels encode the maximum observed prevalence per cancer exactly.
    for _category, group in cancer.groupby("burden_category", sort=False):
        values = group.set_index("proteoform_key").prevalence_p95.to_dict()
        levels = sorted({float(values.get(k, 0)) for k in keys if values.get(k, 0) > 0})
        previous = 0.0
        for level in levels:
            j = nvars
            nvars += 1
            add({j: 1, **{key_index[k]: -1 for k in keys if values.get(k, 0) >= level}}, high=0)
            for metric in coverage_rows:
                coverage_rows[metric][j] = (level - previous) * float(
                    group.iloc[0][f"world_{metric}_pct"]
                )
            previous = level

    def coverage(chosen, metric):
        return sum(
            group.loc[group.proteoform_key.isin(chosen), "prevalence_p95"].max()
            * group.iloc[0][f"world_{metric}_pct"]
            for _, group in cancer.groupby("burden_category")
            if group.proteoform_key.isin(chosen).any()
        )

    baseline_keys = set(retained.proteoform_key)
    baseline = {m: coverage(baseline_keys, m) for m in coverage_rows}

    def evidence_stats(frame):
        strongest = frame.assign(tier=frame.evidence_tier.map(tier_order))
        strongest = strongest.sort_values("tier").drop_duplicates("peptide")
        return {
            "ms_peptides": frame.peptide.nunique(),
            "pmhc": len(frame.drop_duplicates(["peptide", "allele"])),
            "alleles": frame.allele.nunique(),
            "typed_alleles": frame.loc[
                frame.evidence_tier.ne("unrestricted_ms"), "allele"
            ].nunique(),
            "measured_alleles": frame.loc[
                frame.evidence_tier.eq("monoallelic_ms"), "allele"
            ].nunique(),
            "measured_peptides": int(strongest.tier.eq(0).sum()),
            "typed_inferred_peptides": int(strongest.tier.eq(1).sum()),
            "untyped_inferred_peptides": int(strongest.tier.eq(2).sum()),
            "carrier_proxy": carrier_reach(set(frame.allele), frequencies),
        }

    current = evidence_stats(retained)
    current.update(
        {
            "proteoforms": len(baseline_keys),
            "pieces": published["native_pieces"],
            "protein_aa": published["protein_aa"],
            "total_nt": published["total_nt"],
        }
    )
    current["proteins"] = []
    for protein in published["proteins"]:
        key = protein["proteoform_key"]
        own = retained[retained.proteoform_key.eq(key)]
        other = retained[retained.proteoform_key.ne(key)]
        current["proteins"].append(
            {
                "name": protein["name"],
                "aa": protein["assembled_aa"],
                "pieces": protein["assembled_pieces"],
                "ms_peptides": own.peptide.nunique(),
                "unique_ms_peptides": len(set(own.peptide) - set(other.peptide)),
                "mortality_loss_pp": baseline["mortality"]
                - coverage(baseline_keys - {key}, "mortality"),
                "incidence_loss_pp": baseline["incidence"]
                - coverage(baseline_keys - {key}, "incidence"),
                "global_alleles_lost": sorted(set(own.allele) - set(other.allele)),
                "carrier_proxy_loss_pp": 100
                * (
                    current["carrier_proxy"]["reach"]
                    - carrier_reach(set(other.allele), frequencies)["reach"]
                ),
            }
        )
    common_keys = set(
        retained.loc[retained.name.isin(["PRAME", "MAGEA4", "XAGE1A/XAGE1B"]), "proteoform_key"]
    )
    top_three = {f"{m}_lower_pp": coverage(common_keys, m) for m in baseline}
    prame_all = ligands[ligands.name.eq("PRAME")]
    prame_selected = retained[retained.name.eq("PRAME")]
    prame = {
        "retained": evidence_stats(prame_selected),
        "all_supported": evidence_stats(prame_all),
        "additional_peptides": sorted(set(prame_all.peptide) - set(prame_selected.peptide)),
        "additional_target_alleles": sorted(set(prame_all.allele) - set(prame_selected.allele)),
        "additional_global_alleles": sorted(set(prame_all.allele) - required_alleles),
    }
    # Exact priority: more distinct MS peptides; fewer proteins; fewer pieces;
    # then less native sequence. Integer weight bounds prevent priority reversal.
    protein_weight = ns + 1
    peptide_weight = (nk + 1) * protein_weight + ns + 1
    objective = np.zeros(nvars)
    objective[:ns] = [1 + s["length_aa"] / (aa_limit * (ns + 1)) for s in segments]
    for j in key_index.values():
        objective[j] = protein_weight
    for j in pep_index.values():
        objective[j] = -peptide_weight
    provenance = {
        table_path(name).name: sha256(table_path(name).read_bytes()).hexdigest()
        for name in (
            "ranking.csv",
            "specific_intervals.csv",
            "ms_assignments.csv",
            "cancer_summary.csv",
            "ligands.csv",
        )
    }
    provenance["report_data_json"] = sha256(args.report_json.read_bytes()).hexdigest()
    comparisons = []
    for tolerance in args.tolerances:
        test_rows = [*rows, *coverage_rows.values()]
        lower = [*lows, *(baseline[m] - tolerance for m in coverage_rows)]
        upper = [*highs, np.inf, np.inf]
        rr, cc, vv = [], [], []
        for i, row in enumerate(test_rows):
            for j, value in row.items():
                rr.append(i)
                cc.append(j)
                vv.append(value)
        matrix = coo_matrix((vv, (rr, cc)), shape=(len(test_rows), nvars)).tocsc()
        print(mode, "solving tolerance pp", tolerance, flush=True)
        solution = milp(
            objective,
            integrality=np.ones(nvars),
            bounds=Bounds(0, 1),
            constraints=LinearConstraint(matrix, lower, upper),
            options={"time_limit": 45, "mip_rel_gap": 0.0},
        )
        if solution.x is None:
            print("NO SOLUTION", solution.message, flush=True)
            comparisons.append(
                {"tolerance_pp": tolerance, "status": solution.status, "message": solution.message}
            )
            continue
        chosen = [s for i, s in enumerate(segments) if solution.x[i] > 0.5]
        chosen_keys = {s["key"] for s in chosen}
        chosen_peptides = set().union(*(s["peptides"] for s in chosen))
        chosen_alleles = set().union(*(s["alleles"] for s in chosen))
        aa = 1 + sum(s["length_aa"] for s in chosen)
        actual_coverage = {m: coverage(chosen_keys, m) for m in baseline}
        assert aa <= aa_limit
        assert required_alleles <= chosen_alleles
        for requirement in restrictions:
            assert any(
                (requirement["target"] == "global" or s["name"] == requirement["target"])
                and any(
                    h[3] == requirement["allele"] and tier_order[h[4]] <= requirement["tier"]
                    for h in s["hits"]
                )
                for s in chosen
            )
        assert all(actual_coverage[m] + 1e-7 >= baseline[m] - tolerance for m in baseline)
        assert len(chosen_peptides) == round(sum(solution.x[j] for j in pep_index.values()))
        assert len(chosen_keys) == round(sum(solution.x[j] for j in key_index.values()))
        # Recheck every proposed peptide against its unchanged CTA-specific interval.
        for s in chosen:
            protein = rank.loc[s["key"], "sequence"]
            for a, b, peptide, _allele, _tier in s["hits"]:
                assert s["start"] <= a < b <= s["end"]
                assert protein[a:b] == peptide
        details = []
        for key in sorted(chosen_keys):
            pieces = [s for s in chosen if s["key"] == key]
            hits = [h for s in pieces for h in s["hits"]]
            strongest = {}
            order = {"monoallelic_ms": 0, "sample_allele_ms": 1, "unrestricted_ms": 2}
            for _, _, peptide, _allele, tier in hits:
                strongest[peptide] = min(order[tier], strongest.get(peptide, 3))
            details.append(
                {
                    "name": pieces[0]["name"],
                    "aa": sum(s["length_aa"] for s in pieces),
                    "pieces": len(pieces),
                    "ms_peptides": len(strongest),
                    "alleles": len({h[3] for h in hits}),
                    "measured_peptides": sum(t == 0 for t in strongest.values()),
                    "typed_inferred_peptides": sum(t == 1 for t in strongest.values()),
                    "untyped_inferred_peptides": sum(t == 2 for t in strongest.values()),
                }
            )
        evidence = evidence_stats(
            pd.DataFrame(
                [h for s in chosen for h in s["hits"]],
                columns=["start", "end", "peptide", "allele", "evidence_tier"],
            )
        )
        assert evidence["ms_peptides"] == len(chosen_peptides)
        result = {
            "tolerance_pp": tolerance,
            "status": solution.status,
            "message": solution.message,
            "mip_gap": solution.mip_gap,
            "aa_upper_bound_with_start_m": aa,
            "complete_rna_nt_upper_bound": 3 * aa + overhead_nt,
            "proteoforms": len(chosen_keys),
            "pieces": len(chosen),
            "ms_peptides": len(chosen_peptides),
            "alleles": len(chosen_alleles),
            "evidence": evidence,
            "carrier_proxy": carrier_reach(chosen_alleles, frequencies),
            "mortality_lower_pp": actual_coverage["mortality"],
            "incidence_lower_pp": actual_coverage["incidence"],
            "proteins": details,
            "segments": [
                {k: v for k, v in s.items() if k not in {"hits", "peptides", "alleles"}}
                for s in chosen
            ],
        }
        print(
            json.dumps(
                {
                    k: v
                    for k, v in result.items()
                    if k not in {"segments", "proteins", "carrier_proxy"}
                }
            ),
            flush=True,
        )
        print("PRAME", next((p for p in details if p["name"] == "PRAME"), None), flush=True)
        comparisons.append(result)
        (OUTPUT / f"{mode}.json").write_text(
            json.dumps(
                {
                    "mode": mode,
                    "analyzer_source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
                    "calculation_versions": {
                        "python": sys.version.split()[0],
                        "tsarina_source": __version__,
                        **{p: version(p) for p in ("scipy", "numpy", "pandas")},
                    },
                    "input_report_versions": published["versions"],
                    "baseline": baseline,
                    "current": current,
                    "top_three": top_three,
                    "prame": prame,
                    "required_alleles": sorted(required_alleles),
                    "input_sha256": provenance,
                    "allocation_only": True,
                    "global_union_only": args.global_union_only,
                    "max_aa": aa_limit,
                    "max_nt": published["max_nt"],
                    "rna_overhead_nt": overhead_nt,
                    "protected_restrictions": restrictions,
                    "comparisons": comparisons,
                },
                indent=2,
            )
            + "\n"
        )


for mode in ("strict", "loose"):
    analyze_mode(mode)
