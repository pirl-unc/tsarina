"""Portable friendly results, also renderable from verified saved manifests."""

from __future__ import annotations

import gzip
import json
import shutil
from collections import Counter
from datetime import date
from hashlib import sha256
from pathlib import Path
from xml.etree import ElementTree

import pandas as pd

from .regions import allele_frequency_audit
from .vaccine_coverage import carrier_reach, ciwd_frequencies, coverage_tables

FRAME_KEYS = (
    "ranking",
    "cancer_cohorts",
    "cancer_summary",
    "specific_intervals",
    "ligands",
    "funnel",
)


def records(frame):
    return json.loads(frame.to_json(orient="records"))


def load_saved_report(path):
    """Verify every manifest artifact, then restore table types for rendering."""
    path = Path(path)
    path = path / "manifest.json" if path.is_dir() else path
    result = json.loads(path.read_text())
    for name, digest in result["artifact_sha256"].items():
        artifact = path.parent / name
        if not artifact.resolve().is_relative_to(path.parent.resolve()):
            raise ValueError(f"Artifact escapes report directory: {name}")
        if not artifact.is_file() or sha256(artifact.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Saved report artifact failed verification: {name}")
    for key in FRAME_KEYS:
        result[key] = pd.DataFrame(result[key])
    result["source_tables"] = {
        key: pd.DataFrame(value) for key, value in result["source_tables"].items()
    }
    return result


def render_saved_reports(reports, output_dir, *, analysis_date=None):
    """Render any named collection of completed reports without model execution."""
    results = {name: load_saved_report(path) for name, path in reports.items()}
    return write_vaccine_website(results, output_dir, reports=reports, analysis_date=analysis_date)


def _payload(result, destination):
    design = result["design"]
    if not design:
        raise ValueError(
            "A results website requires a construct; consult the scientific failure audit"
        )
    ranking = result["ranking"]
    contributing = result["funnel"][result["funnel"].assembled_aa.gt(0)]
    selected = ranking[ranking.proteoform_key.isin(contributing.proteoform_key)]
    keys = set(selected.proteoform_key)
    layers = pd.DataFrame(design["layers"])
    native = layers[layers.kind.eq("cta_segment")]
    ligands = result["ligands"]
    retained = ligands[ligands.assembled.eq(True)]
    cancer = result["cancer_summary"][result["cancer_summary"].proteoform_key.isin(keys)]
    cohorts = result["cancer_cohorts"][result["cancer_cohorts"].proteoform_key.isin(keys)]
    audit = allele_frequency_audit(result["alleles"])
    frequencies = ciwd_frequencies(audit)
    cumulative, cancer_bounds = coverage_tables(
        contributing, ligands, cancer, design["layers"], frequencies
    )
    cumulative.to_csv(destination / "cumulative-coverage.csv", index=False)
    cancer_bounds.to_csv(destination / "cumulative-cancer.csv", index=False)
    audit.to_csv(destination / "hla_panel.csv", index=False)
    source_tables = result["source_tables"]
    tissue = source_tables.get("protein_ms_map", pd.DataFrame())
    proteins = []
    for p in records(contributing.sort_values("rank")):
        key = p["proteoform_key"]
        hits = retained[retained.proteoform_key.eq(key)]
        p["alleles"] = sorted(set(hits.allele))
        p["layers"] = records(native[native.proteoform_key.eq(key)])
        p["peptides"] = []
        for peptide, rows in hits.groupby("peptide"):
            pmids = {
                v.strip().removesuffix(".0")
                for text in rows.ms_pmids.fillna("")
                for v in str(text).split(";")
                if v.strip().removesuffix(".0").isdigit()
            }
            p["peptides"].append(
                {"peptide": peptide, "alleles": sorted(set(rows.allele)), "pmids": sorted(pmids)}
            )
        p["tissue_hits"] = (
            records(tissue[tissue.proteoform_key.eq(key)]) if not tissue.empty else []
        )
        p["map_file"] = f"protein-map-{sha256(key.encode()).hexdigest()[:12]}.svg"
        proteins.append(p)
    cancers = []
    for category, group in cancer.groupby("burden_category"):
        row = group.iloc[0]
        cancers.append(
            {
                "category": category,
                "mortality_pct": float(row.world_mortality_pct),
                "incidence_pct": float(row.world_incidence_pct),
                "n_samples": int(row.n_samples),
                "proteins": {
                    r["proteoform_key"]: {
                        "prevalence": r["prevalence_p95"],
                        "sample_count": r["n_measured_samples"],
                    }
                    for r in records(group)
                },
            }
        )
    cancers.sort(key=lambda c: -c["mortality_pct"])
    screen = source_tables.get("selection_screen", pd.DataFrame())
    candidate_funnel = source_tables.get("candidate_funnel", pd.DataFrame())
    reasons = Counter(screen.selection_reason) if not screen.empty else Counter()
    pairs = retained.drop_duplicates(["peptide", "allele"])
    config = result["config"]
    data = {
        "config": config,
        "status": result["status"],
        "proteins": proteins,
        "layers": records(layers),
        "cancers": cancers,
        "hla": records(source_tables["hla_support_counts"].merge(audit, on="allele")),
        "tiers": pairs.evidence_tier.value_counts().to_dict(),
        "ms_peptides": retained.peptide.nunique(),
        "pmhc": len(pairs),
        "supported_alleles": retained.allele.nunique(),
        "panel_size": len(result["alleles"]),
        "protein_aa": design["length_aa"],
        "total_nt": design["length_nt"],
        "native_pieces": len(native),
        "linker_aa": int(
            (
                layers[layers.kind.eq("linker")].end_aa - layers[layers.kind.eq("linker")].start_aa
            ).sum()
        ),
        "initial_binders": design["initial_objective"][0],
        "final_binders": design["final_objective"][0],
        "max_aa": config["max_length_aa"],
        "max_nt": config["max_length_nt"],
        "poly_a": config["poly_a_length"],
        "reference_mortality_pct": sum(c["mortality_pct"] for c in cancers),
        "sample_count": int(cohorts.drop_duplicates("cancer_code").n_samples.sum()),
        "cohort_count": int(cohorts.cancer_code.nunique()),
        "cancer_cohort_codes": sorted(set(cohorts.cancer_code)),
        "burden_sources": sorted(set(cancer["source"].dropna())) if "source" in cancer else [],
        "sequences": {
            "protein": design["protein"],
            "cds": design["cds_nt"],
            "full": design["full_nt"],
        },
        "coverage": records(cumulative),
        "cancer_bounds": records(cancer_bounds),
        "panel_reach": carrier_reach(result["alleles"], frequencies),
        "versions": result["versions"],
        "expression_sources": records(source_tables.get("expression_availability", pd.DataFrame())),
        "provenance": result["provenance"],
        "candidate_funnel": records(candidate_funnel),
        "screen": records(screen),
        "ranking": records(ranking.drop(columns=["sequence"], errors="ignore")),
        "filter_counts": {
            "screened": sum(
                count
                for reason, count in reasons.items()
                if reason not in {"not_screened", "non_positive_score", "excluded_gene_pattern"}
            ),
            "no_specific": reasons["no_cta_specific_sequence"],
            "no_normal_filtered": reasons["no_sequence_after_normal_ms"],
            "normal_ms_affected": int(
                candidate_funnel.specific_aa.gt(candidate_funnel.normal_ms_filtered_aa).sum()
            )
            if "normal_ms_filtered_aa" in candidate_funnel
            else 0,
            "no_ms": reasons["no_qualified_panel_ms_ligands"],
            "lower_supported": reasons["lower_rank_supported"],
            "too_long": reasons["minimum_segment_exceeds_length"],
        },
    }
    scientific = {
        "settings": config,
        "model_and_data_versions": result["versions"],
        "source_provenance": result["provenance"],
        "coordinate_system": result["coordinate_system"],
        "objective": result["objective"],
        "prediction_callbacks": result["prediction_callbacks"],
        "population_method": "CIWD Table A2 known-allele proxy, HWE within locus and linkage equilibrium across loci; no panel normalization. Other-source frequencies are excluded.",
        "cancer_method": "Marginal p95 union bounds, no patient-overlap data or clinical coverage calculation.",
        "tissue_method": f"Normal-MS policy: {config.get('normal_ms_policy', 'audit')}. Exclusion uses primary donor-resolved Atlas nonmalignant heart/brain/lung HLA-I records, with {config.get('normal_ms_min_donors', 1)} distinct donor(s). Tissue flags alone do not trigger this gate. Remaining observations are audited, not proof of safe regions.",
    }
    (destination / "method-details.json").write_text(
        json.dumps(scientific, indent=2, allow_nan=False) + "\n"
    )
    from .vaccine_website_figures import coverage_figures

    coverage_figures(data, result, destination)
    data["figure_dimensions"] = {}
    for file in destination.glob("*.svg"):
        box = ElementTree.parse(file).getroot().get("viewBox")
        if box:
            _, _, width, height = map(float, box.split())
            data["figure_dimensions"][file.name] = [width, height]
    # Preserve complete provenance rows while keeping static bundles compact.
    for file in destination.glob("*.csv"):
        if file.stat().st_size >= 1024 * 1024:
            file.with_suffix(".csv.gz").write_bytes(gzip.compress(file.read_bytes(), mtime=0))
            file.unlink()
    data["downloads"] = sorted(p.name for p in destination.iterdir() if p.is_file())
    return data


def write_vaccine_website(results, output_dir, *, reports=None, analysis_date=None):
    """Write a portable static website for named result dictionaries.

    Output must be separate from report roots. Assets use relative URLs and
    need no framework, CDN, server-side code or scientific packages to view.
    Serve the folder over HTTP, or publish it with any static host.
    """
    out = Path(output_dir)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Website output directory must be empty to prevent stale artifacts")
    if analysis_date is not None:
        date.fromisoformat(analysis_date)
    reports = reports or {}
    for path in reports.values():
        root = Path(path)
        root = root.parent if root.is_file() else root
        if out.resolve() == root.resolve() or root.resolve().is_relative_to(out.resolve()):
            raise ValueError("Website output must not overwrite a scientific report directory")
    out.mkdir(parents=True, exist_ok=True)
    assets = Path(__file__).parent / "data" / "vaccine-website"
    shutil.copytree(assets, out, dirs_exist_ok=True)
    payload = {"analysis_date": analysis_date or date.today().isoformat(), "designs": {}}
    for name, result in results.items():
        if not name or any(
            c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_"
            for c in name
        ):
            raise ValueError(
                "Report names must contain only letters, numbers, hyphens or underscores"
            )
        dest = out / "downloads" / name
        dest.mkdir(parents=True, exist_ok=True)
        if name in reports:
            root = Path(reports[name])
            root = root.parent if root.is_file() else root
            for file in root.iterdir():
                if file.is_file() and file.suffix in {".csv", ".fasta", ".svg", ".png"}:
                    shutil.copy2(file, dest / file.name)
        # Single source of truth also works before a report manifest is written.
        for key in FRAME_KEYS:
            result[key].to_csv(dest / f"{key}.csv", index=False)
        for key, table in result["source_tables"].items():
            table.to_csv(dest / f"{key}.csv", index=False)
        if result["design"]:
            design = result["design"]
            for key, value in (
                ("protein", design["protein"]),
                ("cds", design["cds_nt"]),
                ("full", design["full_nt"]),
            ):
                (dest / f"{key}.fasta").write_text(f">cta_vaccine_{key}\n{value}\n")
            for key in ("layers", "cleavage", "junctions"):
                pd.DataFrame(design[key]).to_csv(dest / f"{key}.csv", index=False)
            pd.DataFrame([r for r in design["junctions"] if r["below_threshold"]]).to_csv(
                dest / "junction-binders.csv", index=False
            )
        payload["designs"][name] = _payload(result, dest)
    (out / "data.json").write_text(
        json.dumps(payload, allow_nan=False, separators=(",", ":")) + "\n"
    )
    return out / "index.html"
