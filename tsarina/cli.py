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
    tsarina data register iedb /path/to/file  # register a manual download
    tsarina data path iedb                    # print path to a dataset
    tsarina data remove iedb                  # unregister (keeps the file)
    tsarina build observations [--force]      # build the observations index

    # Curation reference data (HPA RNA/IHC + NCBI gene_info), version-pinned
    tsarina reference list                    # versions + cache status
    tsarina reference fetch                   # fetch all (pinned HPA version)
    tsarina reference fetch hpa_normal_tissue --hpa-version v23
    tsarina reference path hpa_rna_consensus  # print cached path

    # Analysis
    tsarina personalize --hla HLA-A*02:01,... --cta PRAME=87.3 -o targets.csv
    tsarina hits --gene PRAME --allele HLA-A*24:02
    tsarina panel -o panel.csv
"""

from __future__ import annotations

import argparse
import json
import sys

from .downloads import (
    available_datasets,
    data_dir,
    fetch,
    fetch_all_data_assets,
    get_path,
    info,
    list_datasets,
    refresh,
    register,
    remove,
)

# ── data subcommands ────────────────────────────────────────────────────────


def _data_list(args: argparse.Namespace) -> None:
    if getattr(args, "all", False):
        # `--all` is deprecated in favor of the `available` subcommand, matching
        # hitlist's `data list` / `data available` split.
        print(
            "Note: `tsarina data list --all` is deprecated; use `tsarina data available`.",
            file=sys.stderr,
        )
        _data_available(args)
        return
    datasets = list_datasets()
    if not datasets:
        print("No datasets registered.")
        print(f"Data directory: {data_dir()}")
        print("Run 'tsarina data available' to see known datasets.")
        print("HPA/NCBI curation reference data is separate — see `tsarina reference list`.")
        return
    print(f"{'Name':<12} {'Size':>12}  {'Source':<14} Description")
    print("-" * 72)
    for name, meta in sorted(datasets.items()):
        size = meta.get("size_bytes", 0)
        if size > 1_000_000_000:
            size_str = f"{size / 1e9:.1f} GB"
        elif size > 1_000_000:
            size_str = f"{size / 1e6:.1f} MB"
        elif size > 1_000:
            size_str = f"{size / 1e3:.1f} KB"
        else:
            size_str = f"{size} B"
        source = meta.get("source", "")
        source_label = "registered" if source == "registered" else "fetched"
        desc = meta.get("description", "")
        print(f"{name:<12} {size_str:>12}  {source_label:<14} {desc}")
    print(f"\nData directory: {data_dir()}")
    print("`tsarina data available` shows the full catalog. HPA/NCBI curation")
    print("reference data is separate — see `tsarina reference list`.")


def _data_available(args: argparse.Namespace) -> None:
    datasets = available_datasets()
    registered = set(list_datasets().keys())
    print(f"{'Name':<12} {'Status':<12} Description")
    print("-" * 72)
    for name, desc in sorted(datasets.items()):
        status = "installed" if name in registered else ""
        print(f"{name:<12} {status:<12} {desc}")


def _data_register(args: argparse.Namespace) -> None:
    try:
        p = register(args.name, args.path, description=args.description)
        print(f"Registered '{args.name}' -> {p}")
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def _data_fetch(args: argparse.Namespace) -> None:
    try:
        p = fetch(args.name, force=args.force)
        print(f"Ready: {p}")
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def _unknown_dataset(name: str) -> None:
    """Print a friendly 'unknown dataset' error (lists the catalog) and exit 1.

    Matches the `tsarina reference` error style instead of leaking a KeyError repr.
    """
    known = ", ".join(sorted(available_datasets())) or "(none)"
    print(
        f"Error: unknown dataset '{name}'. Known: {known}.\n"
        f"  See `tsarina data available`; fetch with `tsarina data fetch {name}` "
        f"or register a local file with `tsarina data register {name} <path>`.",
        file=sys.stderr,
    )
    sys.exit(1)


def _data_path(args: argparse.Namespace) -> None:
    try:
        p = get_path(args.name)
        print(str(p))
    except KeyError:
        _unknown_dataset(args.name)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def _data_remove(args: argparse.Namespace) -> None:
    if not remove(args.name, delete_file=getattr(args, "delete", False)):
        print(f"No dataset '{args.name}' was registered (nothing to do).", file=sys.stderr)
        sys.exit(1)
    if getattr(args, "delete", False):
        print(f"Unregistered and deleted '{args.name}'.")
    else:
        print(f"Unregistered '{args.name}' (file kept on disk).")


def _data_refresh(args: argparse.Namespace) -> None:
    try:
        p = refresh(args.name)
        print(f"Refreshed: {p}")
    except (ValueError, KeyError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def _data_info(args: argparse.Namespace) -> None:
    try:
        print(json.dumps(info(args.name), indent=2, default=str))
    except KeyError:
        _unknown_dataset(args.name)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def _data_fetch_all(args: argparse.Namespace) -> None:
    paths = fetch_all_data_assets(force=args.force)
    print(f"Fetched {len(paths)} data-asset file(s) into {data_dir()}.")


def _build_observations(args: argparse.Namespace) -> None:
    from .indexing import ensure_index_built

    path = ensure_index_built(force=args.force, verbose=True)
    print(f"Observations index: {path}")


def _data_build(args: argparse.Namespace) -> None:
    # Deprecated: mirror hitlist's `data build` -> `build observations` move.
    print(
        "Note: `tsarina data build` is deprecated; use `tsarina build observations`.",
        file=sys.stderr,
    )
    _build_observations(args)


def _fmt_bytes(size: int | None) -> str:
    if not size:
        return "-"
    for unit, scale in (("GB", 1e9), ("MB", 1e6), ("KB", 1e3)):
        if size >= scale:
            return f"{size / scale:.1f} {unit}"
    return f"{size} B"


# ── reference subcommands (HPA / NCBI curation data) ────────────────────────


def _reference_list(args: argparse.Namespace) -> None:
    from . import reference_data

    rows = reference_data.status()
    print(f"{'Dataset':<20} {'Cached':<8} {'Version':<10} {'Size':>9}  Description")
    print("-" * 92)
    for r in rows:
        cached = "yes" if r["cached"] else "no"
        ver = r["cached_version"] or f"({r['default_version']})"
        print(
            f"{r['name']:<20} {cached:<8} {ver:<10} {_fmt_bytes(r['bytes']):>9}  {r['description']}"
        )
    print(f"\nCache directory: {reference_data.cache_dir()}")
    print(
        f"Default HPA version: {reference_data.DEFAULT_HPA_VERSION} "
        "(only release serving both RNA consensus and normal_tissue)"
    )


def _reference_fetch(args: argparse.Namespace) -> None:
    from . import reference_data

    names = args.names or list(reference_data.REFERENCE_DATASETS)
    # Attempt every dataset; a single flaky URL shouldn't abort a "fetch all".
    failures = []
    for name in names:
        # --hpa-version only applies to HPA datasets; the rolling NCBI dataset
        # has no such version and would otherwise fail a "fetch all" run.
        version = args.hpa_version if reference_data.is_hpa_dataset(name) else None
        try:
            path = reference_data.download(name, version=version, force=args.force)
        except reference_data.ReferenceDataError as e:
            print(f"Error: {name}: {e}", file=sys.stderr)
            failures.append(name)
            continue
        print(f"Ready: {name} -> {path}")
    if failures:
        print(f"Failed to fetch: {', '.join(failures)}", file=sys.stderr)
        sys.exit(1)


def _reference_path(args: argparse.Namespace) -> None:
    from . import reference_data

    version = args.hpa_version if reference_data.is_hpa_dataset(args.name) else None
    try:
        print(reference_data.local_path(args.name, version=version))
    except reference_data.ReferenceDataError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def _build_reference_parser(sub: argparse._SubParsersAction) -> None:
    ref = sub.add_parser(
        "reference",
        help="Versioned HPA/NCBI curation reference data (RNA, IHC, gene_info)",
    )
    ref_sub = ref.add_subparsers(dest="reference_command")
    ref_sub.add_parser("list", help="Show reference datasets, versions, and cache status")
    p_fetch = ref_sub.add_parser("fetch", help="Fetch reference dataset(s); default: all")
    p_fetch.add_argument("names", nargs="*", help="Dataset name(s); default: all")
    p_fetch.add_argument("--hpa-version", default=None, help="HPA release (e.g. v23)")
    p_fetch.add_argument("--force", "-f", action="store_true", help="Re-download")
    p_path = ref_sub.add_parser("path", help="Print the cache path for a dataset")
    p_path.add_argument("name", help="Dataset name")
    p_path.add_argument("--hpa-version", default=None, help="HPA release (e.g. v23)")


def _handle_reference(args: argparse.Namespace) -> None:
    handlers = {
        "list": _reference_list,
        "fetch": _reference_fetch,
        "path": _reference_path,
    }
    handlers[getattr(args, "reference_command", None) or "list"](args)


def _build_data_parser(sub: argparse._SubParsersAction) -> None:
    data_parser = sub.add_parser(
        "data", help="MS evidence datasets (IEDB/CEDAR/viral) + observations index"
    )
    data_sub = data_parser.add_subparsers(dest="data_command")

    p_list = data_sub.add_parser("list", help="Show installed datasets")
    p_list.add_argument(
        "--all", action="store_true", help="Deprecated: use `tsarina data available`."
    )

    data_sub.add_parser("available", help="Show all known datasets (installed or not)")

    p_reg = data_sub.add_parser("register", help="Register a local file")
    p_reg.add_argument("name", help="Dataset name (e.g. iedb, cedar)")
    p_reg.add_argument("path", help="Path to the file")
    p_reg.add_argument("--description", "-d", help="Optional description")

    p_fetch = data_sub.add_parser("fetch", help="Download a fetchable dataset")
    p_fetch.add_argument("name", help="Dataset name (e.g. iedb, hpv16, ebv)")
    p_fetch.add_argument("--force", "-f", action="store_true", help="Re-download")

    p_refresh = data_sub.add_parser("refresh", help="Re-download a fetchable dataset")
    p_refresh.add_argument("name", help="Dataset name")

    p_fetch_all = data_sub.add_parser("fetch-all", help="Fetch all mirrored data-asset files")
    p_fetch_all.add_argument("--force", "-f", action="store_true", help="Re-download")

    p_info = data_sub.add_parser("info", help="Detailed metadata for a dataset (JSON)")
    p_info.add_argument("name", help="Dataset name")

    p_path = data_sub.add_parser("path", help="Print path to a dataset")
    p_path.add_argument("name", help="Dataset name")

    p_rm = data_sub.add_parser("remove", help="Unregister a dataset (keeps file by default)")
    p_rm.add_argument("name", help="Dataset name")
    p_rm.add_argument("--delete", action="store_true", help="Also delete the cached file")

    p_build = data_sub.add_parser(
        "build",
        help="DEPRECATED — use `tsarina build observations`",
    )
    p_build.add_argument("--force", "-f", action="store_true", help="Rebuild from scratch")


def _handle_data(args: argparse.Namespace) -> None:
    handlers = {
        "list": _data_list,
        "available": _data_available,
        "register": _data_register,
        "fetch": _data_fetch,
        "refresh": _data_refresh,
        "fetch-all": _data_fetch_all,
        "info": _data_info,
        "path": _data_path,
        "remove": _data_remove,
        "build": _data_build,
    }
    # Bare `tsarina data` defaults to listing, mirroring `tsarina reference`.
    handlers[getattr(args, "data_command", None) or "list"](args)


def _build_build_parser(sub: argparse._SubParsersAction) -> None:
    build_parser = sub.add_parser(
        "build", help="Build cached artifacts (observations index from IEDB/CEDAR)"
    )
    build_sub = build_parser.add_subparsers(dest="build_command")
    p_obs = build_sub.add_parser(
        "observations",
        help="Build the hitlist observations index (IEDB + CEDAR; cached after first run)",
    )
    p_obs.add_argument("--force", "-f", action="store_true", help="Rebuild from scratch")


def _handle_build(args: argparse.Namespace) -> None:
    handlers = {"observations": _build_observations}
    cmd = getattr(args, "build_command", None)
    if cmd is None:
        print("Usage: tsarina build observations [--force]", file=sys.stderr)
        sys.exit(2)
    handlers[cmd](args)


# ── Main entry point ────────────────────────────────────────────────────────


def main() -> None:
    from . import cli_hits, cli_personalize, cli_spanning
    from .version import __version__

    parser = argparse.ArgumentParser(
        prog="tsarina",
        description="tsarina: cancer-testis antigens, viral targets, and shared cancer immunotherapy peptides",
    )
    parser.add_argument("--version", action="version", version=f"tsarina {__version__}")
    sub = parser.add_subparsers(dest="command")

    _build_data_parser(sub)
    _build_build_parser(sub)
    _build_reference_parser(sub)
    cli_personalize.build_parser(sub)
    cli_hits.build_parser(sub)
    cli_spanning.build_parser(sub)

    argv = sys.argv[1:]
    deprecated_spanning = bool(argv and argv[0] == "spanning")
    if deprecated_spanning:
        argv = ["panel", *argv[1:]]

    args = parser.parse_args(argv)
    if deprecated_spanning:
        args.deprecated_spanning = True

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "data":
        _handle_data(args)
    elif args.command == "build":
        _handle_build(args)
    elif args.command == "reference":
        _handle_reference(args)
    elif args.command == "personalize":
        cli_personalize.handle(args)
    elif args.command == "hits":
        cli_hits.handle(args)
    elif args.command == "panel":
        cli_spanning.handle(args)


if __name__ == "__main__":
    main()
