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
    tsarina data list [--all]                 # installed datasets (--all = full catalog)
    tsarina data fetch hpv16                  # auto-download a viral proteome
    tsarina data register iedb /path/to/file  # register a manual download
    tsarina data path iedb                    # print path to a dataset
    tsarina data remove iedb                  # unregister (keeps the file)
    tsarina data build [--force]              # build the observations index

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
import sys

from .downloads import (
    available_datasets,
    data_dir,
    fetch,
    get_path,
    list_datasets,
    register,
    remove,
)

# ── data subcommands ────────────────────────────────────────────────────────


def _data_list(args: argparse.Namespace) -> None:
    if getattr(args, "all", False):
        _data_list_all(args)
        return
    datasets = list_datasets()
    if not datasets:
        print("No datasets registered.")
        print(f"Data directory: {data_dir()}")
        print("Run 'tsarina data list --all' to see known datasets.")
        print("(HPA/NCBI curation reference data is separate: `tsarina reference list`.)")
        return
    print(f"{'Name':<12} {'Size':>12}  {'Source':<14} Description")
    print("-" * 72)
    for name, info in sorted(datasets.items()):
        size = info.get("size_bytes", 0)
        if size > 1_000_000_000:
            size_str = f"{size / 1e9:.1f} GB"
        elif size > 1_000_000:
            size_str = f"{size / 1e6:.1f} MB"
        elif size > 1_000:
            size_str = f"{size / 1e3:.1f} KB"
        else:
            size_str = f"{size} B"
        source = info.get("source", "")
        source_label = "registered" if source == "registered" else "fetched"
        desc = info.get("description", "")
        print(f"{name:<12} {size_str:>12}  {source_label:<14} {desc}")
    print(f"\nData directory: {data_dir()}")
    print("(`tsarina data list --all` shows the full catalog; HPA/NCBI curation")
    print(" reference data is separate: `tsarina reference list`.)")


def _data_list_all(args: argparse.Namespace) -> None:
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


def _data_path(args: argparse.Namespace) -> None:
    try:
        p = get_path(args.name)
        print(str(p))
    except (KeyError, FileNotFoundError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def _data_remove(args: argparse.Namespace) -> None:
    remove(args.name)
    print(f"Unregistered '{args.name}' (file not deleted).")


def _data_build(args: argparse.Namespace) -> None:
    from .indexing import ensure_index_built

    path = ensure_index_built(force=args.force, verbose=True)
    print(f"Observations index: {path}")


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

    p_list = data_sub.add_parser("list", help="Show installed datasets (--all for full catalog)")
    p_list.add_argument(
        "--all", action="store_true", help="Show the full catalog, not just installed datasets."
    )

    p_reg = data_sub.add_parser("register", help="Register a local file")
    p_reg.add_argument("name", help="Dataset name (e.g. iedb, cedar)")
    p_reg.add_argument("path", help="Path to the file")
    p_reg.add_argument("--description", "-d", help="Optional description")

    p_fetch = data_sub.add_parser("fetch", help="Download a fetchable dataset")
    p_fetch.add_argument("name", help="Dataset name (e.g. hpv16, ebv)")
    p_fetch.add_argument("--force", "-f", action="store_true", help="Re-download")

    p_path = data_sub.add_parser("path", help="Print path to a dataset")
    p_path.add_argument("name", help="Dataset name")

    p_rm = data_sub.add_parser("remove", help="Unregister a dataset (keeps file)")
    p_rm.add_argument("name", help="Dataset name")

    p_build = data_sub.add_parser(
        "build",
        help="Build the hitlist observations index (IEDB + CEDAR; cached after first run)",
    )
    p_build.add_argument("--force", "-f", action="store_true", help="Rebuild from scratch")


def _handle_data(args: argparse.Namespace) -> None:
    handlers = {
        "list": _data_list,
        "register": _data_register,
        "fetch": _data_fetch,
        "path": _data_path,
        "remove": _data_remove,
        "build": _data_build,
    }
    # Bare `tsarina data` defaults to listing, mirroring `tsarina reference`.
    handlers[getattr(args, "data_command", None) or "list"](args)


# ── Main entry point ────────────────────────────────────────────────────────


def main() -> None:
    from . import cli_hits, cli_personalize, cli_spanning

    parser = argparse.ArgumentParser(
        prog="tsarina",
        description="tsarina: cancer-testis antigens, viral targets, and shared cancer immunotherapy peptides",
    )
    sub = parser.add_subparsers(dest="command")

    _build_data_parser(sub)
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
