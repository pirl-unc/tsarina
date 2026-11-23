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

"""Read-only environment status and catalog discovery CLIs.

Exposes:
- ``tsarina status [--json]`` — environment, package versions, OncoRef policy,
  Hitlist index/EvidenceStore state, and predictor availability.
- ``tsarina catalog {ctas,panels,mutations,viruses} [--json]`` — browse the
  built-in CTA, HLA panel, hotspot mutation, and oncogenic virus catalogs
  without dropping into Python.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import sys
from pathlib import Path


def _pkg_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def _plain_table(headers: list[str], rows: list[list[str]]) -> str:
    widths = [
        max(len(h), *(len(r[i]) for r in rows)) if rows else len(h) for i, h in enumerate(headers)
    ]
    lines = [
        "  ".join(h.ljust(w) for h, w in zip(headers, widths)),
        "  ".join("-" * w for w in widths),
    ]
    lines.extend("  ".join(cell.ljust(w) for cell, w in zip(row, widths)) for row in rows)
    return "\n".join(lines)


def _check_mhcflurry_weights() -> dict[str, object]:
    version = _pkg_version("mhcflurry")
    if version is None:
        return {"installed": False, "version": None, "models_ready": False}
    models_ready = False
    models_path = None
    try:
        from mhcflurry.downloads import get_path

        path = get_path("models_class1_presentation", test_exists=False)
        if path and Path(path).exists():
            models_ready = True
            models_path = str(path)
    except Exception:
        models_ready = False
    return {
        "installed": True,
        "version": version,
        "models_ready": models_ready,
        "models_path": models_path,
    }


def _check_pyensembl_releases(
    releases: tuple[tuple[str, int], ...] = (("human", 112), ("mouse", 115)),
) -> dict[str, object]:
    version = _pkg_version("pyensembl")
    if version is None:
        return {"installed": False, "version": None, "cached_releases": {}}
    cached: dict[str, bool] = {}
    try:
        from pyensembl import EnsemblRelease

        for species, release in releases:
            key = f"{species}:{release}"
            try:
                genome = EnsemblRelease(release=release, species=species)
                _ = genome.db
                cached[key] = True
            except Exception:
                cached[key] = False
    except Exception:
        pass
    return {"installed": True, "version": version, "cached_releases": cached}


def collect_status() -> dict[str, object]:
    """Gather read-only status across tsarina and its upstream dependencies."""
    from .cta_selection import CTA_TIERS, CTA_selection, selection_provenance
    from .downloads import data_dir, info, list_datasets
    from .evidence_store import resolve_evidence_package
    from .loader import ms_evidence_provenance
    from .version import __version__

    # OncoRef policy & CTA counts
    sel_prov = selection_provenance()
    cta_df = CTA_selection(category="CTA")
    tier_counts = (
        {
            tier: int((cta_df["cta_tier"] == tier).sum())
            for tier in CTA_TIERS
            if "cta_tier" in cta_df.columns
        }
        if not cta_df.empty
        else {}
    )
    overlay_meta = ms_evidence_provenance()
    overlay_matches = overlay_meta.get("cta_selection_policy_ids") == sel_prov.get(
        "cta_selection_policy_ids"
    ) and overlay_meta.get("cta_selection_source_snapshot_ids") == sel_prov.get(
        "cta_selection_source_snapshot_ids"
    )

    # Hitlist index & registrations (strictly read-only)
    ddir = data_dir()
    obs_path = ddir / "observations.parquet"
    meta_path = ddir / "observations_meta.json"
    obs_meta: dict[str, object] = {}
    if meta_path.exists():
        try:
            obs_meta = json.loads(meta_path.read_text())
        except Exception:
            obs_meta = {}

    registered: dict[str, str] = {}
    registered_symlinks: dict[str, str] = {}
    for name in sorted(list_datasets()):
        try:
            detail = info(name)
            registered[name] = str(detail.get("status", "unknown"))
            sym = detail.get("symlink_status")
            if sym and sym != "regular_file":
                registered_symlinks[name] = str(sym)
        except Exception:
            registered[name] = "unknown"

    ev_pkg = resolve_evidence_package(None)

    return {
        "tsarina_version": __version__,
        "packages": {
            name: _pkg_version(name)
            for name in (
                "oncoref",
                "hitlist",
                "mhcgnomes",
                "pyensembl",
                "mhcflurry",
                "mhctools",
                "pepsickle",
                "orfdata",
                "dnachisel",
            )
        },
        "oncoref": {
            "version": _pkg_version("oncoref"),
            "policy_ids": sel_prov.get("cta_selection_policy_ids", []),
            "reference_ids": sel_prov.get("cta_selection_reference_ids", []),
            "human_cta_count": len(cta_df),
            "human_cta_tiers": tier_counts,
            "ms_overlay_matches_selection": bool(overlay_matches),
            "ms_overlay_policy_ids": overlay_meta.get("cta_selection_policy_ids", []),
        },
        "hitlist": {
            "version": _pkg_version("hitlist"),
            "data_dir": str(ddir),
            "observations_built": obs_path.exists(),
            "observations_rows": obs_meta.get("n_rows"),
            "observations_peptides": obs_meta.get("n_peptides"),
            "observations_artifact_version": obs_meta.get("artifact_version"),
            "evidence_package": str(ev_pkg) if ev_pkg is not None else None,
            "registered_datasets": registered,
            "registered_symlink_status": registered_symlinks,
        },
        "pyensembl": _check_pyensembl_releases(),
        "mhcflurry": _check_mhcflurry_weights(),
    }


def format_status(status: dict[str, object]) -> str:
    """Format :func:`collect_status` for terminal display."""
    oncoref_info = status["oncoref"]
    hitlist_info = status["hitlist"]
    mhcflurry_info = status["mhcflurry"]
    pyensembl_info = status["pyensembl"]
    pkgs = status["packages"]

    policies = ", ".join(oncoref_info["policy_ids"]) or "unknown"
    tiers = ", ".join(f"{k}={v}" for k, v in oncoref_info["human_cta_tiers"].items())
    overlay_state = (
        "current"
        if oncoref_info["ms_overlay_matches_selection"]
        else f"predates selection ({', '.join(oncoref_info['ms_overlay_policy_ids']) or 'unknown'})"
    )

    if hitlist_info["observations_built"]:
        rows = hitlist_info["observations_rows"]
        rows_str = f"{int(rows):,} rows" if rows is not None else "built"
        obs_str = f"ready ({rows_str}, artifact v{hitlist_info['observations_artifact_version']})"
    else:
        obs_str = "not built (run `tsarina build observations`)"

    reg = hitlist_info["registered_datasets"]
    reg_sym = hitlist_info.get("registered_symlink_status", {})
    reg_str = (
        ", ".join(f"{k}={v} ({reg_sym[k]})" if k in reg_sym else f"{k}={v}" for k, v in reg.items())
        if reg
        else "none"
    )

    flurry_str = (
        f"{mhcflurry_info['version']} (weights ready)"
        if mhcflurry_info["installed"] and mhcflurry_info["models_ready"]
        else f"{mhcflurry_info['version']} (weights not downloaded)"
        if mhcflurry_info["installed"]
        else "not installed"
    )

    cached_rel = pyensembl_info.get("cached_releases", {})
    ensembl_str = (
        f"{pyensembl_info['version']} ("
        + ", ".join(f"{k}={'ready' if v else 'missing'}" for k, v in cached_rel.items())
        + ")"
        if pyensembl_info["installed"]
        else "not installed"
    )

    optional_str = ", ".join(
        f"{name}={pkgs[name] or 'not installed'}"
        for name in ("mhcgnomes", "mhctools", "pepsickle", "orfdata", "dnachisel")
    )

    lines = [
        f"tsarina {status['tsarina_version']}",
        f"  OncoRef:        {oncoref_info['version']} | policy {policies}",
        f"  Human CTAs:     {oncoref_info['human_cta_count']} ({tiers})",
        f"  MS overlay:     {overlay_state}",
        f"  Hitlist:        {hitlist_info['version']} | data_dir={hitlist_info['data_dir']}",
        f"  Observations:   {obs_str}",
        f"  EvidenceStore:  {hitlist_info['evidence_package'] or 'not configured'}",
        f"  Registered:     {reg_str}",
        f"  PyEnsembl:      {ensembl_str}",
        f"  MHCflurry:      {flurry_str}",
        f"  Other deps:     {optional_str}",
    ]
    return "\n".join(lines)


def build_status_parser(sub: argparse._SubParsersAction) -> argparse.ArgumentParser:
    p = sub.add_parser(
        "status",
        help="Show installed versions, OncoRef policy, Hitlist index, and predictor readiness",
        description=(
            "Inspect tsarina's local environment (read-only; downloads and builds nothing): "
            "active OncoRef selection policy, Hitlist observations/EvidenceStore status, "
            "cached Ensembl releases, and MHCflurry model readiness."
        ),
    )
    p.add_argument("--json", action="store_true", help="Machine-readable JSON output")
    return p


def handle_status(args: argparse.Namespace) -> None:
    status = collect_status()
    if args.json:
        print(json.dumps(status, indent=2, default=str))
    else:
        print(format_status(status))


# ── tsarina catalog ─────────────────────────────────────────────────────────


def _catalog_overview(args: argparse.Namespace) -> None:
    from .alleles import panel_names
    from .cta_selection import CTA_TIERS, CTA_selection
    from .fusions import cancer_fusions
    from .mutations import HOTSPOT_MUTATIONS
    from .viral import ONCOGENIC_VIRUSES, orfdata_viral_subsets

    ctas = CTA_selection(category="CTA")
    fusions_df = cancer_fusions()
    defining_fusions = (
        int(fusions_df["is_defining"].astype(bool).sum())
        if "is_defining" in fusions_df.columns
        else 0
    )
    try:
        orfdata_subsets = [
            {"name": s["name"], "virus_count": s["virus_count"]} for s in orfdata_viral_subsets()
        ]
    except ImportError:
        orfdata_subsets = []

    summary = {
        "ctas": {
            "human_cta_count": len(ctas),
            "tiers": {
                t: int((ctas["cta_tier"] == t).sum())
                for t in CTA_TIERS
                if "cta_tier" in ctas.columns
            },
        },
        "panels": list(panel_names()),
        "mutations_count": len(HOTSPOT_MUTATIONS),
        "fusions_count": len(fusions_df),
        "defining_fusions_count": defining_fusions,
        "viruses": sorted(ONCOGENIC_VIRUSES),
        "orfdata_viral_subsets": orfdata_subsets,
    }
    if getattr(args, "json", False):
        print(json.dumps(summary, indent=2))
        return
    tiers_str = ", ".join(f"{k}={v}" for k, v in summary["ctas"]["tiers"].items())
    print("Built-in tsarina catalogs:")
    print(f"  ctas       {summary['ctas']['human_cta_count']} human CTAs ({tiers_str})")
    print(f"  panels     {len(summary['panels'])} HLA panels ({', '.join(summary['panels'])})")
    print(f"  mutations  {summary['mutations_count']} curated hotspot mutations")
    print(
        f"  fusions    {summary['fusions_count']} curated cancer fusions "
        f"({summary['defining_fusions_count']} defining)"
    )
    vir_extra = ""
    if orfdata_subsets:
        subset_desc = ", ".join(f"{s['name']}={s['virus_count']}" for s in orfdata_subsets)
        vir_extra = f"; + orfdata GDC subsets: {subset_desc}"
    print(
        f"  viruses    {len(summary['viruses'])} built-in oncogenic viruses "
        f"({', '.join(summary['viruses'])}){vir_extra}"
    )
    print("\nInspect any catalog with `tsarina catalog <ctas|panels|mutations|fusions|viruses>`.")


def _catalog_ctas(args: argparse.Namespace) -> None:
    from .cta_selection import CTA_selection, CTA_selection_coverage, CTA_sex_restricted_antigens

    species = args.species
    is_human = species.lower() in ("human", "homo sapiens", "homo_sapiens")
    sex = getattr(args, "sex", None)
    sex_elig = getattr(args, "sex_eligibility", "all_adult")

    category = args.category
    coverage_info: dict[str, object] | None = None
    if not is_human:
        coverage_info = CTA_selection_coverage(species=species)
    if sex is not None:
        try:
            df = CTA_sex_restricted_antigens(sex=sex, eligibility=sex_elig, species=species)
        except ValueError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            sys.exit(1)
        if category is not None:
            df = df.loc[df["category"].eq(category)].copy()
    else:
        if category is None and (
            is_human or (coverage_info and coverage_info.get("selection_status") == "evaluated")
        ):
            category = "CTA"
        df = CTA_selection(category=category, species=species)
    if args.tier:
        if "cta_tier" in df.columns:
            df = df[df["cta_tier"] == args.tier].copy()
        else:
            df = df.iloc[:0].copy()

    # Attach bundled MS safety restriction for human CTAs when available.
    ms_lookup: dict[str, str] = {}
    if is_human and not df.empty:
        try:
            from .evidence import CTA_evidence

            ev = CTA_evidence()
            if {"Symbol", "ms_restriction"} <= set(ev.columns):
                ms_lookup = dict(zip(ev["Symbol"].astype(str), ev["ms_restriction"].astype(str)))
        except Exception:
            ms_lookup = {}

    if "label" in df.columns:
        sort_cols = (
            ["sex_restriction_class", "label"]
            if sex is not None and "sex_restriction_class" in df.columns
            else ["cta_tier", "label"]
            if "cta_tier" in df.columns
            else ["label"]
        )
        df = df.sort_values(sort_cols)

    total_rows = len(df)
    if args.top is not None:
        if args.top <= 0:
            print("Error: --top must be a positive integer.", file=sys.stderr)
            sys.exit(1)
        df = df.head(args.top)

    records: list[dict[str, object]] = []
    for row in df.itertuples(index=False):
        label = str(getattr(row, "label", "") or getattr(row, "symbol", ""))
        gene_ids = getattr(row, "gene_ids", None)
        if isinstance(gene_ids, (list, tuple)):
            gene_ids_str = ";".join(str(g) for g in gene_ids)
        else:
            gene_ids_str = "" if gene_ids is None else str(gene_ids)
        rec = {
            "label": label,
            "category": str(getattr(row, "category", "") or ""),
            "tier": str(getattr(row, "cta_tier", "") or getattr(row, "subtype", "") or ""),
            "gene_ids": gene_ids_str,
        }
        if sex is not None:
            rec["sex_restriction_class"] = str(getattr(row, "sex_restriction_class", "") or "")
        if ms_lookup:
            rec["ms_restriction"] = ms_lookup.get(label, "")
        records.append(rec)

    if args.json:
        payload: dict[str, object] = {
            "species": species,
            "total": total_rows,
            "shown": len(records),
            "ctas": records,
        }
        if sex is not None:
            payload["sex"] = sex
            payload["sex_eligibility"] = sex_elig
        if not is_human:
            payload["coverage"] = (
                coverage_info
                if coverage_info is not None
                else CTA_selection_coverage(species=species)
            )
        print(json.dumps(payload, indent=2, default=str))
        return

    if not records:
        print(f"No CTAs matched (species={species}, category={category}, tier={args.tier}).")
        return

    headers = ["label", "category", "tier"]
    if sex is not None:
        headers.append("sex_restriction_class")
    headers.append("gene_ids")
    if ms_lookup:
        headers.append("ms_restriction")
    rows = [[str(r.get(h, "")) for h in headers] for r in records]
    print(_plain_table(headers, rows))
    if total_rows > len(records):
        print(f"\nShowing {len(records)} of {total_rows} rows (pass --top to adjust).")


def _catalog_panels(args: argparse.Namespace) -> None:
    from .alleles import PANEL_DEFINITIONS, get_panel, panel_names
    from .regions import allele_frequency_audit

    if args.panel:
        if args.panel not in PANEL_DEFINITIONS:
            print(
                f"Error: unknown panel '{args.panel}'. Available: {', '.join(panel_names())}.",
                file=sys.stderr,
            )
            sys.exit(1)
        alleles = get_panel(args.panel)
        freq_df = allele_frequency_audit(alleles).set_index("allele")
        records = [
            {
                "allele": a,
                "weighted_allele_frequency": float(freq_df.loc[a, "coverage_frequency"])
                if a in freq_df.index
                else 0.0,
                "frequency_source": str(freq_df.loc[a, "coverage_frequency_source"])
                if a in freq_df.index
                else "missing",
            }
            for a in alleles
        ]
        if args.json:
            print(json.dumps({"panel": args.panel, "alleles": records}, indent=2))
            return
        rows = [
            [
                r["allele"],
                f"{100.0 * float(r['weighted_allele_frequency']):.1f}%",
                str(r["frequency_source"]),
            ]
            for r in records
        ]
        print(_plain_table(["allele", "freq", "source"], rows))
        return

    records = []
    for name in panel_names():
        alleles = get_panel(name)
        n_a = sum(1 for a in alleles if a.startswith("HLA-A*"))
        n_b = sum(1 for a in alleles if a.startswith("HLA-B*"))
        n_c = sum(1 for a in alleles if a.startswith("HLA-C*"))
        records.append(
            {
                "panel": name,
                "allele_count": len(alleles),
                "loci": f"A={n_a}, B={n_b}, C={n_c}",
                "sample_alleles": ", ".join(alleles[:4]) + ("..." if len(alleles) > 4 else ""),
            }
        )
    if args.json:
        print(json.dumps(records, indent=2))
        return
    rows = [
        [str(r["panel"]), str(r["allele_count"]), r["loci"], r["sample_alleles"]] for r in records
    ]
    print(_plain_table(["panel", "alleles", "loci", "examples"], rows))


def _catalog_mutations(args: argparse.Namespace) -> None:
    from .mutations import HOTSPOT_MUTATIONS

    if args.json:
        print(json.dumps(HOTSPOT_MUTATIONS, indent=2))
        return
    rows = [
        [
            str(m["label"]),
            str(m["gene"]),
            f"{m['ref_aa']}{m['protein_position']}{m['alt_aa']}",
            ", ".join(m.get("cancer_types", [])),
            str(m.get("frequency_note", "")),
        ]
        for m in HOTSPOT_MUTATIONS
    ]
    print(_plain_table(["label", "gene", "change", "cancers", "notes"], rows))


def _catalog_fusions(args: argparse.Namespace) -> None:
    from .fusions import cancer_fusions

    df = cancer_fusions(
        cancer_code=args.cancer,
        gene=args.gene,
        fusion_family=args.family,
        defining_only=args.defining_only,
        rnaseq_detectable_only=args.rnaseq_only,
    )
    total_rows = len(df)
    if args.top is not None:
        if args.top <= 0:
            print("Error: --top must be a positive integer.", file=sys.stderr)
            sys.exit(1)
        df = df.head(args.top)

    records: list[dict[str, object]] = []
    for row in df.itertuples(index=False):
        g5 = str(getattr(row, "gene_5prime", "") or "")
        g3 = str(getattr(row, "gene_3prime", "") or "")
        records.append(
            {
                "fusion": f"{g5}-{g3}",
                "cancer_code": str(getattr(row, "cancer_code", "") or ""),
                "fusion_family": str(getattr(row, "fusion_family", "") or ""),
                "gene_5prime": g5,
                "gene_3prime": g3,
                "frequency": str(getattr(row, "frequency", "") or ""),
                "is_defining": bool(getattr(row, "is_defining", False)),
                "rnaseq_detectable": str(getattr(row, "rnaseq_detectable", "") or ""),
                "mechanism": str(getattr(row, "mechanism", "") or ""),
                "confidence": str(getattr(row, "confidence", "") or ""),
                "pmid": str(getattr(row, "pmid", "") or ""),
            }
        )

    if args.json:
        print(
            json.dumps(
                {"total": total_rows, "shown": len(records), "fusions": records},
                indent=2,
            )
        )
        return

    if not records:
        print("No cancer fusions matched the requested filters.")
        return

    headers = ["fusion", "cancer", "family", "freq", "defining", "rnaseq", "mechanism"]
    rows = [
        [
            str(r["fusion"]),
            str(r["cancer_code"]),
            str(r["fusion_family"]),
            str(r["frequency"]),
            "yes" if r["is_defining"] else "no",
            str(r["rnaseq_detectable"]),
            str(r["mechanism"]),
        ]
        for r in records
    ]
    print(_plain_table(headers, rows))
    if total_rows > len(records):
        print(f"\nShowing {len(records)} of {total_rows} rows (pass --top to adjust).")


def _catalog_viruses(args: argparse.Namespace) -> None:
    from .viral import ONCOGENIC_VIRUSES, orfdata_viral_subsets

    source = getattr(args, "source", "builtin")
    subset_filter = getattr(args, "subset", None)
    if subset_filter is not None and source == "builtin":
        source = "orfdata"

    builtin_records = [{"key": key, **meta} for key, meta in sorted(ONCOGENIC_VIRUSES.items())]

    if source == "builtin":
        if args.json:
            print(json.dumps(builtin_records, indent=2))
            return
        rows = [
            [
                str(r["key"]),
                str(r["name"]),
                str(r["taxonomy_id"]),
                ", ".join(r.get("key_oncoproteins", [])),
                ", ".join(r.get("cancers", [])),
            ]
            for r in builtin_records
        ]
        print(_plain_table(["key", "name", "taxon", "oncoproteins", "cancers"], rows))
        try:
            subsets = orfdata_viral_subsets()
        except ImportError:
            subsets = ()
        if subsets:
            desc = ", ".join(f"{s['name']} ({s['virus_count']})" for s in subsets)
            print(f"\nTip: pass --source orfdata for IARC GDC viral subsets: {desc}.")
        return

    try:
        orf_subsets = list(orfdata_viral_subsets())
    except ImportError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

    if subset_filter is not None:
        orf_subsets = [s for s in orf_subsets if s["name"] == subset_filter]
        if not orf_subsets:
            valid = ", ".join(str(s["name"]) for s in orfdata_viral_subsets())
            print(
                f"Error: unknown orfdata viral subset '{subset_filter}'. Available: {valid}.",
                file=sys.stderr,
            )
            sys.exit(1)

    if args.json:
        if source == "orfdata":
            print(json.dumps(orf_subsets, indent=2))
        else:
            print(
                json.dumps(
                    {"builtin": builtin_records, "orfdata_subsets": orf_subsets},
                    indent=2,
                )
            )
        return

    if source == "all":
        rows = [
            [
                str(r["key"]),
                str(r["name"]),
                str(r["taxonomy_id"]),
                ", ".join(r.get("key_oncoproteins", [])),
                ", ".join(r.get("cancers", [])),
            ]
            for r in builtin_records
        ]
        print("Built-in UniProt oncogenic viruses:")
        print(_plain_table(["key", "name", "taxon", "oncoproteins", "cancers"], rows))
        print("\norfdata IARC GDC viral subsets:")

    orf_rows = [
        [
            str(s["name"]),
            str(s["virus_count"]),
            str(s["title"]),
            ", ".join(s["abbreviations"][:6]) + ("..." if len(s["abbreviations"]) > 6 else ""),
        ]
        for s in orf_subsets
    ]
    print(_plain_table(["subset", "viruses", "title", "examples"], orf_rows))
    if subset_filter is not None and orf_subsets:
        print("\nMembers: " + ", ".join(orf_subsets[0]["abbreviations"]))


def build_catalog_parser(sub: argparse._SubParsersAction) -> argparse.ArgumentParser:
    p = sub.add_parser(
        "catalog",
        help="Browse built-in CTA, HLA panel, mutation, fusion, and virus catalogs",
        description=(
            "Inspect the built-in target and allele catalogs that feed `tsarina panel`, "
            "`tsarina personalize`, and `tsarina vaccine`."
        ),
    )
    p.add_argument("--json", action="store_true", help="Machine-readable JSON output")
    p.set_defaults(func=_catalog_overview)
    cat_sub = p.add_subparsers(dest="catalog_command")

    p_ctas = cat_sub.add_parser("ctas", help="List current CTAs and tiers (human, mouse, dog)")
    p_ctas.add_argument(
        "--species", default="human", help="Species (default: human; e.g. mouse, dog)"
    )
    p_ctas.add_argument(
        "--category",
        default=None,
        help="Filter by OncoRef selection category (default: 'CTA' for human, all for non-human)",
    )
    p_ctas.add_argument(
        "--tier",
        choices=("core", "core_with_extended_exposure", "extended_only"),
        default=None,
        help="Filter human CTAs by cta_tier",
    )
    p_ctas.add_argument(
        "--sex",
        choices=("female", "male"),
        default=None,
        help="Filter evaluated-species antigens to those eligible for this biological sex",
    )
    p_ctas.add_argument(
        "--sex-eligibility",
        choices=("all_adult", "sex_organ_only"),
        default="all_adult",
        help="Eligibility policy when --sex is given (default: all_adult)",
    )
    p_ctas.add_argument("--top", type=int, default=50, help="Maximum rows to show (default 50)")
    p_ctas.add_argument("--json", action="store_true", help="Machine-readable JSON output")
    p_ctas.set_defaults(func=_catalog_ctas)

    p_panels = cat_sub.add_parser("panels", help="List built-in population HLA panels and alleles")
    p_panels.add_argument(
        "--panel", default=None, help="Show all alleles in a specific named panel"
    )
    p_panels.add_argument("--json", action="store_true", help="Machine-readable JSON output")
    p_panels.set_defaults(func=_catalog_panels)

    p_muts = cat_sub.add_parser("mutations", help="List curated hotspot cancer mutations")
    p_muts.add_argument("--json", action="store_true", help="Machine-readable JSON output")
    p_muts.set_defaults(func=_catalog_mutations)

    p_fus = cat_sub.add_parser("fusions", help="List curated OncoRef cancer fusions")
    p_fus.add_argument("--cancer", default=None, help="Filter by cancer code, e.g. SARC_EWS, LUAD")
    p_fus.add_argument("--gene", default=None, help="Filter by 5-prime or 3-prime gene symbol")
    p_fus.add_argument("--family", default=None, help="Filter by fusion family, e.g. EWSR1-ETS")
    p_fus.add_argument(
        "--defining-only",
        action="store_true",
        help="Show only cancer-defining fusions",
    )
    p_fus.add_argument(
        "--rnaseq-only",
        action="store_true",
        help="Show only directly RNA-seq detectable fusions",
    )
    p_fus.add_argument("--top", type=int, default=50, help="Maximum rows to show, default 50")
    p_fus.add_argument("--json", action="store_true", help="Machine-readable JSON output")
    p_fus.set_defaults(func=_catalog_fusions)

    p_vir = cat_sub.add_parser("viruses", help="List supported oncogenic viruses and GDC subsets")
    p_vir.add_argument(
        "--source",
        choices=("builtin", "orfdata", "all"),
        default="builtin",
        help="Catalog source: builtin UniProt panel, orfdata IARC GDC subsets, or all",
    )
    p_vir.add_argument(
        "--subset",
        default=None,
        help="Inspect a specific orfdata viral subset, e.g. cancer-causing",
    )
    p_vir.add_argument("--json", action="store_true", help="Machine-readable JSON output")
    p_vir.set_defaults(func=_catalog_viruses)

    return p
