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
    for name in sorted(list_datasets()):
        try:
            registered[name] = str(info(name).get("status", "unknown"))
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
    reg_str = ", ".join(f"{k}={v}" for k, v in reg.items()) if reg else "none"

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
    from .mutations import HOTSPOT_MUTATIONS
    from .viral import ONCOGENIC_VIRUSES

    ctas = CTA_selection(category="CTA")
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
        "viruses": sorted(ONCOGENIC_VIRUSES),
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
        f"  viruses    {len(summary['viruses'])} oncogenic viruses ({', '.join(summary['viruses'])})"
    )
    print("\nInspect any catalog with `tsarina catalog <ctas|panels|mutations|viruses>`.")


def _catalog_ctas(args: argparse.Namespace) -> None:
    from .cta_selection import CTA_selection, CTA_selection_coverage

    species = args.species
    category = args.category
    if category is None and species.lower() in ("human", "homo sapiens", "homo_sapiens"):
        category = "CTA"

    df = CTA_selection(category=category, species=species)
    if args.tier:
        if "cta_tier" in df.columns:
            df = df[df["cta_tier"] == args.tier].copy()
        else:
            df = df.iloc[:0].copy()

    # Attach bundled MS safety restriction for human CTAs when available.
    ms_lookup: dict[str, str] = {}
    if species.lower() in ("human", "homo sapiens", "homo_sapiens") and not df.empty:
        try:
            from .evidence import CTA_evidence

            ev = CTA_evidence()
            if {"Symbol", "ms_restriction"} <= set(ev.columns):
                ms_lookup = dict(zip(ev["Symbol"].astype(str), ev["ms_restriction"].astype(str)))
        except Exception:
            ms_lookup = {}

    if "label" in df.columns:
        df = df.sort_values(["cta_tier", "label"] if "cta_tier" in df.columns else ["label"])

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
        if species.lower() not in ("human", "homo sapiens", "homo_sapiens"):
            payload["coverage"] = CTA_selection_coverage(species=species)
        print(json.dumps(payload, indent=2, default=str))
        return

    if not records:
        print(f"No CTAs matched (species={species}, category={category}, tier={args.tier}).")
        return

    headers = ["label", "category", "tier", "gene_ids"]
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


def _catalog_viruses(args: argparse.Namespace) -> None:
    from .viral import ONCOGENIC_VIRUSES

    records = [{"key": key, **meta} for key, meta in sorted(ONCOGENIC_VIRUSES.items())]
    if args.json:
        print(json.dumps(records, indent=2))
        return
    rows = [
        [
            str(r["key"]),
            str(r["name"]),
            str(r["taxonomy_id"]),
            ", ".join(r.get("key_oncoproteins", [])),
            ", ".join(r.get("cancers", [])),
        ]
        for r in records
    ]
    print(_plain_table(["key", "name", "taxon", "oncoproteins", "cancers"], rows))


def build_catalog_parser(sub: argparse._SubParsersAction) -> argparse.ArgumentParser:
    p = sub.add_parser(
        "catalog",
        help="Browse built-in CTA, HLA panel, hotspot mutation, and oncogenic virus catalogs",
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

    p_vir = cat_sub.add_parser("viruses", help="List supported oncogenic viruses")
    p_vir.add_argument("--json", action="store_true", help="Machine-readable JSON output")
    p_vir.set_defaults(func=_catalog_viruses)

    return p
