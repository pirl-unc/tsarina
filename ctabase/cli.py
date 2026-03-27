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

"""CLI for ctabase data management.

Usage::

    ctabase-data list                           # show registered datasets
    ctabase-data available                      # show all known datasets
    ctabase-data register iedb /path/to/file    # register a manual download
    ctabase-data fetch hpv16                    # auto-download a viral proteome
    ctabase-data fetch hpv16 --force            # re-download
    ctabase-data path iedb                      # print path to registered file
    ctabase-data remove iedb                    # unregister (keeps the file)
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


def _cmd_list(args: argparse.Namespace) -> None:
    datasets = list_datasets()
    if not datasets:
        print("No datasets registered.")
        print(f"Data directory: {data_dir()}")
        print("Run 'ctabase-data available' to see known datasets.")
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


def _cmd_available(args: argparse.Namespace) -> None:
    datasets = available_datasets()
    registered = set(list_datasets().keys())
    print(f"{'Name':<12} {'Status':<12} Description")
    print("-" * 72)
    for name, desc in sorted(datasets.items()):
        status = "installed" if name in registered else ""
        print(f"{name:<12} {status:<12} {desc}")


def _cmd_register(args: argparse.Namespace) -> None:
    try:
        p = register(args.name, args.path, description=args.description)
        print(f"Registered '{args.name}' -> {p}")
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def _cmd_fetch(args: argparse.Namespace) -> None:
    try:
        p = fetch(args.name, force=args.force)
        print(f"Ready: {p}")
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def _cmd_path(args: argparse.Namespace) -> None:
    try:
        p = get_path(args.name)
        print(str(p))
    except (KeyError, FileNotFoundError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def _cmd_remove(args: argparse.Namespace) -> None:
    remove(args.name)
    print(f"Unregistered '{args.name}' (file not deleted).")


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="ctabase-data",
        description="ctabase data management",
    )
    sub = parser.add_subparsers(dest="command")

    # list
    sub.add_parser("list", help="Show registered datasets")

    # available
    sub.add_parser("available", help="Show all known datasets")

    # register
    p_reg = sub.add_parser("register", help="Register a local file")
    p_reg.add_argument("name", help="Dataset name (e.g. iedb, cedar)")
    p_reg.add_argument("path", help="Path to the file")
    p_reg.add_argument("--description", "-d", help="Optional description")

    # fetch
    p_fetch = sub.add_parser("fetch", help="Download a fetchable dataset")
    p_fetch.add_argument("name", help="Dataset name (e.g. hpv16, ebv)")
    p_fetch.add_argument("--force", "-f", action="store_true", help="Re-download")

    # path
    p_path = sub.add_parser("path", help="Print path to a dataset")
    p_path.add_argument("name", help="Dataset name")

    # remove
    p_rm = sub.add_parser("remove", help="Unregister a dataset (keeps file)")
    p_rm.add_argument("name", help="Dataset name")

    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
        sys.exit(1)

    handlers = {
        "list": _cmd_list,
        "available": _cmd_available,
        "register": _cmd_register,
        "fetch": _cmd_fetch,
        "path": _cmd_path,
        "remove": _cmd_remove,
    }
    handlers[args.command](args)


if __name__ == "__main__":
    main()
