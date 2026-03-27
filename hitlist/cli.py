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

"""CLI for hitlist.

Usage::

    hitlist data list                           # show registered datasets
    hitlist data available                      # show all known datasets
    hitlist data register iedb /path/to/file    # register a manual download
    hitlist data fetch hpv16                    # auto-download a viral proteome
    hitlist data fetch hpv16 --force            # re-download
    hitlist data path iedb                      # print path to registered file
    hitlist data remove iedb                    # unregister (keeps the file)
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
    datasets = list_datasets()
    if not datasets:
        print("No datasets registered.")
        print(f"Data directory: {data_dir()}")
        print("Run 'hitlist data available' to see known datasets.")
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


def _build_data_parser(sub: argparse._SubParsersAction) -> None:
    data_parser = sub.add_parser("data", help="Manage external datasets")
    data_sub = data_parser.add_subparsers(dest="data_command")

    data_sub.add_parser("list", help="Show registered datasets")
    data_sub.add_parser("available", help="Show all known datasets")

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


def _handle_data(args: argparse.Namespace) -> None:
    handlers = {
        "list": _data_list,
        "available": _data_available,
        "register": _data_register,
        "fetch": _data_fetch,
        "path": _data_path,
        "remove": _data_remove,
    }
    if args.data_command is None:
        print("Usage: hitlist data {list,available,register,fetch,path,remove}", file=sys.stderr)
        sys.exit(1)
    handlers[args.data_command](args)


# ── Main entry point ────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="hitlist",
        description="hitlist: cancer-testis antigens, viral targets, and shared cancer immunotherapy peptides",
    )
    sub = parser.add_subparsers(dest="command")

    _build_data_parser(sub)

    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "data":
        _handle_data(args)


if __name__ == "__main__":
    main()
