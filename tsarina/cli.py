# PYTHON_ARGCOMPLETE_OK
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

"""CLI for tsarina.

Usage::

    # MS evidence datasets (IEDB / CEDAR / viral proteomes + the index)
    tsarina data list                         # installed datasets
    tsarina data available                    # full catalog (installed or not)
    tsarina data fetch hpv16                  # auto-download a viral proteome
    tsarina data preflight hpv16              # read-only fetch plan (JSON)
    tsarina data register iedb /path/to/file  # register a manual download
    tsarina data path iedb                    # print path to a dataset
    tsarina data remove iedb                  # unregister (keeps the file)
    tsarina build observations [--force]      # build the observations index

    # Environment & discovery
    tsarina status                             # installed versions, policy, index & models
    tsarina catalog ctas                       # browse built-in CTAs, panels, mutations, viruses

    # Analysis
    tsarina personalize --hla HLA-A*02:01,... --cta PRAME=87.3 -o targets.csv
    tsarina hits --gene PRAME --allele HLA-A*24:02
    tsarina panel -o panel.csv
"""

from __future__ import annotations

import argparse
import contextlib
import difflib
import functools
import json
import os
import sys
from pathlib import Path

from datacache.cli_help import HelpArgumentParser, enable_argcomplete

from .downloads import (
    DATASET_CATEGORIES,
    FETCHABLE_DATASETS,
    MANUAL_DATASETS,
    OfflineError,
    available_dataset_catalog,
    available_datasets,
    data_dir,
    fetch,
    fetch_all_data_assets,
    get_path,
    info,
    is_offline,
    list_datasets,
    preflight,
    prune_cache_bundles,
    refresh,
    register,
    remove,
    resolve_data_dir,
    set_offline,
    verify_cache_files,
)

# ── data subcommands ────────────────────────────────────────────────────────


def _fail(message: str) -> None:
    """Print ``Error: <message>`` to stderr and exit 1 (no traceback)."""
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(1)


def _format_size(size: int) -> str:
    """Format byte count using decimal SI units without boundary rollover quirks (#226)."""
    if size < 0:
        return f"-{_format_size(-size)}"
    if size < 1000:
        return f"{size} B"
    value = float(size)
    for unit in ("KB", "MB", "GB", "TB"):
        value /= 1000.0
        if round(value, 1) < 1000.0 or unit == "TB":
            return f"{value:.1f} {unit}"
    return f"{value:.1f} TB"


def _unknown_dataset(name: str) -> None:
    """Exit 1 with the closest catalog matches for an unknown dataset name."""
    close = difflib.get_close_matches(name, available_datasets(), n=5, cutoff=0.5)
    suggestion = f" Did you mean: {', '.join(close)}?" if close else ""
    _fail(
        f"unknown dataset '{name}'.{suggestion}\n"
        "  See `tsarina data available`; fetch with `tsarina data fetch <name>` "
        "or register a local file with `tsarina data register <name> <path>`."
    )


def _require_known(name: str) -> None:
    if name not in available_datasets():
        _unknown_dataset(name)


def _guarded(handler):
    """Run a ``data`` handler so every failure is one ``Error:`` line, exit 1."""

    @functools.wraps(handler)
    def run(args: argparse.Namespace) -> None:
        try:
            with _cache_root_scope(getattr(args, "cache_root", None)):
                handler(args)
        except KeyError as e:
            _unknown_dataset(str(e.args[0]) if e.args else "")
        except (ValueError, OSError, RuntimeError, OfflineError) as e:
            _fail(str(e))

    return run


@contextlib.contextmanager
def _cache_root_scope(cache_root: str | Path | None):
    if not cache_root:
        yield
        return
    import hitlist.downloads as hitlist_downloads

    prev_openvax = os.environ.get("OPENVAX_DATA_CACHE")
    prev_hitlist = os.environ.get("HITLIST_DATA_DIR")
    prev_override = getattr(hitlist_downloads, "_override_data_dir", None)
    os.environ["OPENVAX_DATA_CACHE"] = str(Path(cache_root).expanduser())
    os.environ.pop("HITLIST_DATA_DIR", None)
    if prev_override is not None:
        hitlist_downloads._override_data_dir = None
    hitlist_downloads.reset_data_dir_cache()
    try:
        yield
    finally:
        if prev_openvax is None:
            os.environ.pop("OPENVAX_DATA_CACHE", None)
        else:
            os.environ["OPENVAX_DATA_CACHE"] = prev_openvax
        if prev_hitlist is None:
            os.environ.pop("HITLIST_DATA_DIR", None)
        else:
            os.environ["HITLIST_DATA_DIR"] = prev_hitlist
        if prev_override is not None:
            hitlist_downloads._override_data_dir = prev_override
        hitlist_downloads.reset_data_dir_cache()


@contextlib.contextmanager
def _offline_scope(offline: bool):
    if not offline:
        yield
        return
    prev_env = os.environ.get("DATACACHE_OFFLINE")
    os.environ["DATACACHE_OFFLINE"] = "1"
    set_offline(True)
    try:
        yield
    finally:
        set_offline(None)
        if prev_env is None:
            os.environ.pop("DATACACHE_OFFLINE", None)
        else:
            os.environ["DATACACHE_OFFLINE"] = prev_env


def _dir_size_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    if path.is_file():
        try:
            return path.stat().st_size
        except OSError:
            return 0
    total = 0
    for child in path.rglob("*"):
        if child.is_file():
            try:
                total += child.stat().st_size
            except OSError:
                continue
    return total


def _registered_entry(name: str, meta: dict) -> dict:
    """One registered dataset as a plain dict matching the uniform JSON schema (#223/#224).

    ``status`` is datacache's file state (via hitlist's ``info``): ``installed``,
    ``missing``, ``corrupt`` or ``inaccessible``. When hitlist>=2.12.0 reports
    ``symlink_status`` / ``resolved_path``, they are included verbatim.
    """
    path = meta.get("path")
    exists = bool(path and Path(str(path)).exists())
    detail = info(name)
    status = detail.get("status", "installed" if exists else "missing")
    entry: dict[str, object] = {
        "name": name,
        "status": status,
        "exists": exists,
        "size_bytes": meta.get("size_bytes", 0),
        "path": str(path) if path is not None else None,
        "source": "registered" if meta.get("source") == "registered" else "fetched",
        "description": meta.get("description", ""),
    }
    if "symlink_status" in detail:
        entry["symlink_status"] = detail["symlink_status"]
    if "resolved_path" in detail:
        entry["resolved_path"] = detail["resolved_path"]
    return entry


def _data_list(args: argparse.Namespace) -> None:
    if getattr(args, "all", False):
        # `--all` is deprecated in favor of the `available` subcommand, matching
        # hitlist's `data list` / `data available` split.
        print(
            "Note: `tsarina data list --all` is deprecated and will be removed in tsarina 2.0; "
            "use `tsarina data available`.",
            file=sys.stderr,
        )
        _data_available(args)
        return
    entries = [_registered_entry(n, m) for n, m in sorted(list_datasets().items())]
    if getattr(args, "json", False):
        print(json.dumps(entries, indent=2))
        return
    if not entries:
        print("No datasets registered.")
        print(f"Data directory: {data_dir()}")
        print("Run 'tsarina data available' to see known datasets.")
        return
    name_width = max(12, max(len(e["name"]) for e in entries))
    print(f"{'Name':<{name_width}} {'Status':<10} {'Size':>10}  {'Source':<11} Description")
    print("-" * max(72, name_width + 44))
    for e in entries:
        print(
            f"{e['name']:<{name_width}} {e['status']:<10} {_format_size(e['size_bytes']):>10}  "
            f"{e['source']:<11} {e['description']}"
        )
    if any(e["status"] != "installed" for e in entries):
        print(
            "\nA status other than 'installed' means the registered file is unusable. "
            "Re-register it with `tsarina data register`, or run "
            "`tsarina data fetch <name>` to download it again."
        )
    print(f"\nData directory: {data_dir()}")
    print("`tsarina data available` shows the full catalog.")


def _data_available(args: argparse.Namespace) -> None:
    category = getattr(args, "category", None)
    grouped = getattr(args, "grouped", False)
    catalog = available_dataset_catalog(category=category)
    registered = {n: _registered_entry(n, m) for n, m in list_datasets().items()}
    entries = []
    for name, cat_meta in sorted(catalog.items()):
        desc = str(cat_meta.get("description") or "")
        if name in registered:
            entry = dict(registered[name])
            if not entry.get("description"):
                entry["description"] = desc
        else:
            source = (
                "fetchable"
                if name in FETCHABLE_DATASETS
                else "manual"
                if name in MANUAL_DATASETS
                else "catalog"
            )
            entry = {
                "name": name,
                "status": "not installed",
                "exists": False,
                "size_bytes": 0,
                "path": None,
                "source": source,
                "description": desc,
            }
        for key in ("category", "group", "species", "kind"):
            if key in cat_meta:
                entry[key] = cat_meta[key]
        entries.append(entry)
    if getattr(args, "json", False):
        print(json.dumps(entries, indent=2))
        return
    name_width = max(12, max((len(e["name"]) for e in entries), default=12))
    if grouped:
        by_cat: dict[str, list[dict]] = {}
        for e in entries:
            by_cat.setdefault(str(e.get("category") or "other"), []).append(e)
        ordered_cats = [c for c in DATASET_CATEGORIES if c in by_cat] + [
            c for c in sorted(by_cat) if c not in DATASET_CATEGORIES
        ]
        for idx, cat_name in enumerate(ordered_cats):
            cat_entries = by_cat[cat_name]
            if idx:
                print()
            print(f"[{cat_name}] ({len(cat_entries)})")
            for e in cat_entries:
                status = "" if e["status"] == "not installed" else f" [{e['status']}]"
                print(f"  {e['name']:<{name_width}} {e['description']}{status}")
        return
    print(f"{'Name':<{name_width}} {'Status':<14} Description")
    print("-" * max(72, name_width + 28))
    for e in entries:
        status = "" if e["status"] == "not installed" else e["status"]
        print(f"{e['name']:<{name_width}} {status:<14} {e['description']}")


@_guarded
def _data_register(args: argparse.Namespace) -> None:
    p = register(args.name, args.path, description=args.description)
    print(f"Registered '{args.name}' -> {p}")


def _fetch_preview(name: str, *, force: bool = False) -> None:
    """Describe what ``fetch`` would do via ``hitlist.preflight``, without downloading."""
    if name in MANUAL_DATASETS:
        print(f"'{name}' is a manual download: {MANUAL_DATASETS[name]['download_url']}")
        print(f"Register it with: tsarina data register {name} <path>")
        return

    plan = preflight(name, force=force)
    desc = available_datasets().get(name, "")
    header = f"{name} ({desc})" if desc else name
    print(f"Dataset:       {header}")
    if name in FETCHABLE_DATASETS:
        print(f"Source:        {FETCHABLE_DATASETS[name]['url']}")
    else:
        kind = plan.get("kind") or "catalog asset"
        print(f"Source:        managed {kind} (hitlist catalog)")
    print(f"Destination:   {plan['path']}")
    files = plan.get("files")
    if isinstance(files, list) and files:
        n_download = sum(1 for f in files if f.get("would_download"))
        print(f"Bundle files:  {len(files)} ({n_download} would download)")
    expected_size = plan.get("expected_size_bytes")
    if expected_size is not None:
        print(f"Expected size: {_format_size(int(expected_size))}")
    free = plan.get("free_bytes")
    if free is not None:
        space_note = ""
        if plan.get("enough_space") is False:
            space_note = " (insufficient space)"
        print(
            f"Free space:    {_format_size(int(free))} on the data directory's volume{space_note}"
        )

    registered = list_datasets().get(name)
    if not plan.get("would_download"):
        status = "already downloaded; `fetch` would reuse it (--force re-downloads)"
    elif plan.get("exists") and force:
        status = "already downloaded; --force would re-download"
    elif registered and (state := _registered_entry(name, registered)["status"]) != "installed":
        status = f"registered file {registered['path']} is {state}; would download"
    elif plan.get("stale"):
        status = "cached file is stale; would download"
    else:
        status = "not downloaded; would download"
    print(f"Status:        {status}")

    terms = plan.get("terms")
    if isinstance(terms, str) and terms.strip():
        first_line = next(line.strip() for line in terms.splitlines() if line.strip())
        print(f"Terms:         {first_line}")


@_guarded
def _data_dirs(args: argparse.Namespace) -> None:
    import datacache

    from .version import __version__

    resolved_dir, origin = resolve_data_dir()
    exists = resolved_dir.exists()
    size_bytes = _dir_size_bytes(resolved_dir)
    offline_enabled = bool(is_offline())
    offline_src = datacache.offline_source()
    payload: dict[str, object] = {
        "package": "tsarina",
        "version": __version__,
        "data_dir": str(resolved_dir),
        "origin": origin,
        "exists": exists,
        "size_bytes": size_bytes,
        "offline": offline_enabled,
        "offline_source": offline_src,
    }
    if getattr(args, "json", False):
        print(json.dumps(payload, indent=2))
        return
    print(f"tsarina data directory (v{__version__}; delegated to hitlist)")
    state = _format_size(size_bytes) if exists else "not created"
    print(f"  data_dir : {resolved_dir} [{origin}; {state}]")
    if offline_enabled:
        detail = f" ({offline_src})" if offline_src else ""
        print(f"  offline  : enabled{detail}")
    else:
        print("  offline  : disabled")


@_guarded
def _data_verify(args: argparse.Namespace) -> None:
    report = verify_cache_files(
        include_unregistered=bool(getattr(args, "include_unregistered", False)),
    )
    if getattr(args, "json", False):
        print(json.dumps(report, indent=2, default=str))
    else:
        print(
            f"Verified {report['n_verified']}/{report['n_files']} file(s); "
            f"{report['n_failed']} failed."
        )
        for row in report["files"]:
            if row["status"] in ("corrupt", "invalid", "inaccessible") or (
                row["status"] == "missing" and "registered dataset" in row.get("kinds", [])
            ):
                print(f"  FAIL  {row['status']:<12} {row['path']}")
    if not report["ok"]:
        sys.exit(1)


@_guarded
def _data_prune(args: argparse.Namespace) -> None:
    apply = bool(getattr(args, "yes", False))
    keep = int(getattr(args, "keep", 1))
    report = prune_cache_bundles(keep=keep, apply=apply)
    if getattr(args, "json", False):
        print(json.dumps(report, indent=2, default=str))
        return
    print(f"tsarina data cache: {report['root']}")
    if report["n_prunable_bundles"] == 0:
        print("Nothing to prune.")
        return
    size_text = _format_size(int(report["reclaimable_bytes"]))
    if report["dry_run"]:
        print(
            f"Dry run — would remove {report['n_prunable_bundles']} superseded bundle(s) "
            f"across {report['n_stores']} store(s), freeing {size_text}. "
            "Re-run with --yes to delete."
        )
    else:
        print(
            f"Deleted {report['n_deleted_bundles']} superseded bundle(s) "
            f"across {report['n_stores']} store(s), freed {size_text}."
        )


@_guarded
def _data_fetch(args: argparse.Namespace) -> None:
    _require_known(args.name)
    if args.dry_run:
        _fetch_preview(args.name, force=args.force)
        return
    registered = list_datasets().get(args.name)
    if registered and (state := _registered_entry(args.name, registered)["status"]) != "installed":
        print(
            f"Note: registered file {registered['path']} is {state}; downloading a replacement.",
            file=sys.stderr,
        )
    with _offline_scope(bool(getattr(args, "offline", False))):
        print(f"Ready: {fetch(args.name, force=args.force)}")


@_guarded
def _data_preflight(args: argparse.Namespace) -> None:
    _require_known(args.name)
    plan = preflight(args.name, force=args.force)
    print(json.dumps(plan, indent=2, default=str))


@_guarded
def _data_path(args: argparse.Namespace) -> None:
    _require_known(args.name)
    resolved = str(get_path(args.name))
    if getattr(args, "json", False):
        print(json.dumps({"name": args.name, "path": resolved}, indent=2))
    else:
        print(resolved)


@_guarded
def _data_remove(args: argparse.Namespace) -> None:
    delete = getattr(args, "delete", False)
    meta = list_datasets().get(args.name)
    if meta and delete:
        path = meta.get("path")
        size = _format_size(meta.get("size_bytes", 0))
        print(f"Deleting {path} ({size}); other projects may still use this file.", file=sys.stderr)
    if not remove(args.name, delete_file=delete):
        _fail(f"no dataset '{args.name}' was registered (nothing to do).")
    if delete:
        print(f"Unregistered and deleted '{args.name}'.")
    else:
        print(f"Unregistered '{args.name}' (file kept on disk).")


@_guarded
def _data_refresh(args: argparse.Namespace) -> None:
    _require_known(args.name)
    with _offline_scope(bool(getattr(args, "offline", False))):
        print(f"Refreshed: {refresh(args.name)}")


@_guarded
def _data_info(args: argparse.Namespace) -> None:
    _require_known(args.name)
    meta = dict(info(args.name))
    path = meta.get("path")
    exists = bool(path and Path(str(path)).exists())
    meta.setdefault("name", args.name)
    meta["exists"] = exists
    meta.setdefault("size_bytes", 0)
    meta.setdefault("path", str(path) if path is not None else None)
    meta.setdefault(
        "source",
        "fetchable"
        if args.name in FETCHABLE_DATASETS
        else "manual"
        if args.name in MANUAL_DATASETS
        else "catalog",
    )
    meta.setdefault("description", available_datasets().get(args.name, ""))
    print(json.dumps(meta, indent=2, default=str))


@_guarded
def _data_fetch_all(args: argparse.Namespace) -> None:
    with _offline_scope(bool(getattr(args, "offline", False))):
        paths = fetch_all_data_assets(force=args.force)
    print(f"Fetched {len(paths)} data-asset file(s).")


def _build_observations(args: argparse.Namespace) -> None:
    from .indexing import ensure_index_built

    path = ensure_index_built(force=args.force, verbose=True, report_current=True)
    print(f"Observations index: {path}")


def _data_build(args: argparse.Namespace) -> None:
    # Deprecated: mirror hitlist's `data build` -> `build observations` move.
    print(
        "Note: `tsarina data build` is deprecated and will be removed in tsarina 2.0; "
        "use `tsarina build observations`.",
        file=sys.stderr,
    )
    _build_observations(args)


def _add_cache_root_arg(parser: argparse.ArgumentParser, *, suppress_default: bool = True) -> None:
    parser.add_argument(
        "--cache-root",
        default=argparse.SUPPRESS if suppress_default else None,
        help="Override the shared OpenVax cache root for this command; sets OPENVAX_DATA_CACHE.",
    )


def _build_data_parser(sub: argparse._SubParsersAction) -> None:
    data_parser = sub.add_parser(
        "data",
        help="Datasets: list/fetch/register IEDB, CEDAR, viral proteomes, HPA/DepMap and mirrored assets",
    )
    _add_cache_root_arg(data_parser, suppress_default=False)
    # Bare `tsarina data` defaults to listing.
    data_parser.set_defaults(func=_data_list)
    data_sub = data_parser.add_subparsers(dest="data_command")

    p_dirs = data_sub.add_parser(
        "dirs", help="Show cache directory location, size, and offline status (read-only)"
    )
    _add_cache_root_arg(p_dirs)
    p_dirs.add_argument("--json", action="store_true", help="Machine-readable output")
    p_dirs.set_defaults(func=_data_dirs)

    p_list = data_sub.add_parser("list", help="Show installed datasets")
    _add_cache_root_arg(p_list)
    p_list.add_argument(
        "--all", action="store_true", help="Deprecated: use `tsarina data available`."
    )
    p_list.add_argument("--json", action="store_true", help="Machine-readable output")
    p_list.set_defaults(func=_data_list)

    p_avail = data_sub.add_parser("available", help="Show all known datasets (installed or not)")
    _add_cache_root_arg(p_avail)
    p_avail.add_argument(
        "--category",
        choices=DATASET_CATEGORIES,
        default=None,
        help="Filter datasets to one structured category",
    )
    p_avail.add_argument(
        "--grouped",
        action="store_true",
        help="Group datasets by structured category",
    )
    p_avail.add_argument("--json", action="store_true", help="Machine-readable output")
    p_avail.set_defaults(func=_data_available)

    p_reg = data_sub.add_parser("register", help="Register a local file")
    _add_cache_root_arg(p_reg)
    p_reg.add_argument("name", help="Dataset name (e.g. iedb, cedar; see `data available`)")
    p_reg.add_argument(
        "path",
        help="Path to the downloaded file (IEDB/CEDAR: the MHC ligand full CSV export)",
    )
    p_reg.add_argument("--description", "-d", help="Optional description")
    p_reg.set_defaults(func=_data_register)

    p_fetch = data_sub.add_parser("fetch", help="Download a fetchable dataset")
    _add_cache_root_arg(p_fetch)
    p_fetch.add_argument("name", help="Dataset name (e.g. iedb, hpv16, ebv)")
    p_fetch.add_argument("--force", "-f", action="store_true", help="Re-download")
    p_fetch.add_argument(
        "--dry-run",
        action="store_true",
        help="Show source, destination and free disk space; download nothing",
    )
    p_fetch.add_argument(
        "--offline",
        action="store_true",
        help="Disallow network downloads (also honors OPENVAX_OFFLINE / DATACACHE_OFFLINE)",
    )
    p_fetch.set_defaults(func=_data_fetch)

    p_preflight = data_sub.add_parser(
        "preflight",
        help="Read-only fetch plan (destination, expected size, free space, terms) as JSON",
    )
    _add_cache_root_arg(p_preflight)
    p_preflight.add_argument("name", help="Dataset name (e.g. iedb, hpv16, uniprot/human, depmap)")
    p_preflight.add_argument(
        "--force", "-f", action="store_true", help="Plan as if --force were passed to fetch"
    )
    p_preflight.add_argument(
        "--json", action="store_true", help="Machine-readable output (default)"
    )
    p_preflight.set_defaults(func=_data_preflight)

    p_refresh = data_sub.add_parser("refresh", help="Re-download a fetchable dataset")
    _add_cache_root_arg(p_refresh)
    p_refresh.add_argument("name", help="Dataset name")
    p_refresh.add_argument(
        "--offline",
        action="store_true",
        help="Disallow network downloads (also honors OPENVAX_OFFLINE / DATACACHE_OFFLINE)",
    )
    p_refresh.set_defaults(func=_data_refresh)

    p_fetch_all = data_sub.add_parser("fetch-all", help="Fetch all mirrored data-asset files")
    _add_cache_root_arg(p_fetch_all)
    p_fetch_all.add_argument("--force", "-f", action="store_true", help="Re-download")
    p_fetch_all.add_argument(
        "--offline",
        action="store_true",
        help="Disallow network downloads (also honors OPENVAX_OFFLINE / DATACACHE_OFFLINE)",
    )
    p_fetch_all.set_defaults(func=_data_fetch_all)

    p_info = data_sub.add_parser("info", help="Detailed metadata for a dataset (JSON)")
    _add_cache_root_arg(p_info)
    p_info.add_argument("name", help="Dataset name")
    p_info.add_argument("--json", action="store_true", help="Machine-readable output (default)")
    p_info.set_defaults(func=_data_info)

    p_path = data_sub.add_parser("path", help="Print path to a dataset")
    _add_cache_root_arg(p_path)
    p_path.add_argument("name", help="Dataset name")
    p_path.add_argument("--json", action="store_true", help="Machine-readable output")
    p_path.set_defaults(func=_data_path)

    p_verify = data_sub.add_parser(
        "verify",
        help="Verify cached/registered files and trusted mirrored-asset hashes (read-only)",
    )
    _add_cache_root_arg(p_verify)
    p_verify.add_argument(
        "--include-unregistered",
        action="store_true",
        help="Include unregistered files found in the cache directory",
    )
    p_verify.add_argument("--json", action="store_true", help="Machine-readable output")
    p_verify.set_defaults(func=_data_verify)

    p_prune = data_sub.add_parser(
        "prune", help="Prune superseded datacache bundle generations; dry-run by default"
    )
    _add_cache_root_arg(p_prune)
    p_prune.add_argument(
        "--keep",
        type=int,
        default=1,
        help="Number of newest generations to keep per store (default: 1)",
    )
    p_prune.add_argument(
        "--yes", action="store_true", help="Actually delete superseded bundles (default is dry-run)"
    )
    p_prune.add_argument("--json", action="store_true", help="Machine-readable output")
    p_prune.set_defaults(func=_data_prune)

    p_rm = data_sub.add_parser("remove", help="Unregister a dataset (keeps file by default)")
    _add_cache_root_arg(p_rm)
    p_rm.add_argument("name", help="Dataset name")
    p_rm.add_argument("--delete", action="store_true", help="Also delete the cached file")
    p_rm.set_defaults(func=_data_remove)

    p_build = data_sub.add_parser(
        "build",
        help="DEPRECATED — use `tsarina build observations`",
    )
    _add_cache_root_arg(p_build)
    p_build.add_argument("--force", "-f", action="store_true", help="Rebuild from scratch")
    p_build.set_defaults(func=_data_build)


def _build_usage_error(args: argparse.Namespace) -> None:
    print("Usage: tsarina build observations [--force]", file=sys.stderr)
    sys.exit(2)


def _build_build_parser(sub: argparse._SubParsersAction) -> None:
    build_parser = sub.add_parser(
        "build", help="Build cached artifacts (observations index from IEDB/CEDAR)"
    )
    build_parser.set_defaults(func=_build_usage_error)
    build_sub = build_parser.add_subparsers(dest="build_command")
    p_obs = build_sub.add_parser(
        "observations",
        help="Build the hitlist observations index (IEDB + CEDAR; cached after first run)",
    )
    p_obs.add_argument("--force", "-f", action="store_true", help="Rebuild from scratch")
    p_obs.set_defaults(func=_build_observations)


# ── Main entry point ────────────────────────────────────────────────────────


def build_parser() -> HelpArgumentParser:
    from . import cli_catalog, cli_hits, cli_personalize, cli_spanning, cli_vaccine
    from .version import __version__

    parser = HelpArgumentParser(
        prog="tsarina",
        description="tsarina: cancer-testis antigens, viral targets, and shared cancer immunotherapy peptides",
    )
    parser.add_argument("--version", action="version", version=f"tsarina {__version__}")
    sub = parser.add_subparsers(dest="command")

    cli_catalog.build_status_parser(sub).set_defaults(func=cli_catalog.handle_status)
    cli_catalog.build_catalog_parser(sub)
    _build_data_parser(sub)
    _build_build_parser(sub)
    cli_personalize.build_parser(sub).set_defaults(func=cli_personalize.handle)
    cli_hits.build_parser(sub).set_defaults(func=cli_hits.handle)
    cli_spanning.build_parser(sub).set_defaults(func=cli_spanning.handle)
    cli_vaccine.build_parser(sub).set_defaults(func=cli_vaccine.handle)
    cli_vaccine.build_report_parser(sub).set_defaults(func=cli_vaccine.handle_report)
    return parser


def main() -> None:
    parser = build_parser()
    enable_argcomplete(parser)

    argv = sys.argv[1:]
    deprecated_spanning = bool(argv and argv[0] == "spanning")
    if deprecated_spanning:
        argv = ["panel", *argv[1:]]

    args = parser.parse_args(argv)
    if deprecated_spanning:
        args.deprecated_spanning = True

    func = getattr(args, "func", None)
    if func is None:
        parser.print_help()
        sys.exit(1)
    with _cache_root_scope(getattr(args, "cache_root", None)):
        func(args)


if __name__ == "__main__":
    main()
