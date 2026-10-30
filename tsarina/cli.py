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

    # Analysis
    tsarina personalize --hla HLA-A*02:01,... --cta PRAME=87.3 -o targets.csv
    tsarina hits --gene PRAME --allele HLA-A*24:02
    tsarina panel -o panel.csv
"""

from __future__ import annotations

import argparse
import difflib
import functools
import json
import shutil
import sys

from hitlist.cli_help import ColorArgumentParser

from .downloads import (
    FETCHABLE_DATASETS,
    MANUAL_DATASETS,
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


def _fail(message: str) -> None:
    """Print ``Error: <message>`` to stderr and exit 1 (no traceback)."""
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(1)


def _format_size(size: int) -> str:
    for unit, scale in (("GB", 1e9), ("MB", 1e6), ("KB", 1e3)):
        if size >= scale:
            return f"{size / scale:.1f} {unit}"
    return f"{size} B"


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
            handler(args)
        except KeyError as e:
            _unknown_dataset(str(e.args[0]) if e.args else "")
        except (ValueError, OSError, RuntimeError) as e:
            _fail(str(e))

    return run


def _registered_entry(name: str, meta: dict) -> dict:
    """One registered dataset as a plain dict.

    ``status`` is datacache's file state (via hitlist's ``info``): ``installed``,
    ``missing``, ``corrupt`` or ``inaccessible``.
    """
    path = meta.get("path")
    return {
        "name": name,
        "status": info(name)["status"],
        "size_bytes": meta.get("size_bytes", 0),
        "source": "registered" if meta.get("source") == "registered" else "fetched",
        "path": path,
        "description": meta.get("description", ""),
    }


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
    entries = [_registered_entry(n, m) for n, m in sorted(list_datasets().items())]
    if getattr(args, "json", False):
        print(json.dumps(entries, indent=2))
        return
    if not entries:
        print("No datasets registered.")
        print(f"Data directory: {data_dir()}")
        print("Run 'tsarina data available' to see known datasets.")
        return
    print(f"{'Name':<12} {'Status':<10} {'Size':>10}  {'Source':<11} Description")
    print("-" * 72)
    for e in entries:
        print(
            f"{e['name']:<12} {e['status']:<10} {_format_size(e['size_bytes']):>10}  "
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
    registered = {n: _registered_entry(n, m) for n, m in list_datasets().items()}
    entries = [
        {
            "name": name,
            "status": registered[name]["status"] if name in registered else "not installed",
            "description": desc,
        }
        for name, desc in sorted(available_datasets().items())
    ]
    if getattr(args, "json", False):
        print(json.dumps(entries, indent=2))
        return
    print(f"{'Name':<12} {'Status':<14} Description")
    print("-" * 72)
    for e in entries:
        status = "" if e["status"] == "not installed" else e["status"]
        print(f"{e['name']:<12} {status:<14} {e['description']}")


@_guarded
def _data_register(args: argparse.Namespace) -> None:
    p = register(args.name, args.path, description=args.description)
    print(f"Registered '{args.name}' -> {p}")


def _fetch_preview(name: str) -> None:
    """Describe what ``fetch`` would do, without downloading."""
    if name in FETCHABLE_DATASETS:
        ds = FETCHABLE_DATASETS[name]
        dest = data_dir() / ds["filename"]
        # Stopgap until datacache provides a read-only fetch plan (openvax/datacache#135).
        free = shutil.disk_usage(data_dir()).free
        registered = list_datasets().get(name)
        print(f"Dataset:     {name} ({ds['description']})")
        print(f"Source:      {ds['url']}")
        print(f"Destination: {dest}")
        print(f"Free space:  {_format_size(free)} on the data directory's volume")
        if dest.exists():
            print("Status:      already downloaded; `fetch` would reuse it (--force re-downloads)")
        elif registered and (state := _registered_entry(name, registered)["status"]) != "installed":
            print(f"Status:      registered file {registered['path']} is {state}; would download")
        else:
            print("Status:      not downloaded; would download")
    elif name in MANUAL_DATASETS:
        print(f"'{name}' is a manual download: {MANUAL_DATASETS[name]['download_url']}")
        print(f"Register it with: tsarina data register {name} <path>")
    else:
        print(f"'{name}' is a mirrored asset or built dataset; no preview is available.")


@_guarded
def _data_fetch(args: argparse.Namespace) -> None:
    _require_known(args.name)
    if args.dry_run:
        _fetch_preview(args.name)
        return
    registered = list_datasets().get(args.name)
    if registered and (state := _registered_entry(args.name, registered)["status"]) != "installed":
        print(
            f"Note: registered file {registered['path']} is {state}; downloading a replacement.",
            file=sys.stderr,
        )
    print(f"Ready: {fetch(args.name, force=args.force)}")


@_guarded
def _data_path(args: argparse.Namespace) -> None:
    _require_known(args.name)
    print(str(get_path(args.name)))


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
    print(f"Refreshed: {refresh(args.name)}")


@_guarded
def _data_info(args: argparse.Namespace) -> None:
    _require_known(args.name)
    print(json.dumps(info(args.name), indent=2, default=str))


@_guarded
def _data_fetch_all(args: argparse.Namespace) -> None:
    paths = fetch_all_data_assets(force=args.force)
    print(f"Fetched {len(paths)} data-asset file(s).")


def _build_observations(args: argparse.Namespace) -> None:
    from .indexing import ensure_index_built

    path = ensure_index_built(force=args.force, verbose=True, report_current=True)
    print(f"Observations index: {path}")


def _data_build(args: argparse.Namespace) -> None:
    # Deprecated: mirror hitlist's `data build` -> `build observations` move.
    print(
        "Note: `tsarina data build` is deprecated; use `tsarina build observations`.",
        file=sys.stderr,
    )
    _build_observations(args)


def _build_data_parser(sub: argparse._SubParsersAction) -> None:
    data_parser = sub.add_parser(
        "data", help="MS evidence datasets (IEDB/CEDAR/viral) + observations index"
    )
    data_sub = data_parser.add_subparsers(dest="data_command")

    p_list = data_sub.add_parser("list", help="Show installed datasets")
    p_list.add_argument(
        "--all", action="store_true", help="Deprecated: use `tsarina data available`."
    )
    p_list.add_argument("--json", action="store_true", help="Machine-readable output")

    p_avail = data_sub.add_parser("available", help="Show all known datasets (installed or not)")
    p_avail.add_argument("--json", action="store_true", help="Machine-readable output")

    p_reg = data_sub.add_parser("register", help="Register a local file")
    p_reg.add_argument("name", help="Dataset name (e.g. iedb, cedar)")
    p_reg.add_argument("path", help="Path to the file")
    p_reg.add_argument("--description", "-d", help="Optional description")

    p_fetch = data_sub.add_parser("fetch", help="Download a fetchable dataset")
    p_fetch.add_argument("name", help="Dataset name (e.g. iedb, hpv16, ebv)")
    p_fetch.add_argument("--force", "-f", action="store_true", help="Re-download")
    p_fetch.add_argument(
        "--dry-run",
        action="store_true",
        help="Show source, destination and free disk space; download nothing",
    )

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
    # Bare `tsarina data` defaults to listing.
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
    from . import cli_hits, cli_personalize, cli_spanning, cli_vaccine
    from .version import __version__

    parser = ColorArgumentParser(
        prog="tsarina",
        description="tsarina: cancer-testis antigens, viral targets, and shared cancer immunotherapy peptides",
    )
    parser.add_argument("--version", action="version", version=f"tsarina {__version__}")
    sub = parser.add_subparsers(dest="command")

    _build_data_parser(sub)
    _build_build_parser(sub)
    cli_personalize.build_parser(sub)
    cli_hits.build_parser(sub)
    cli_spanning.build_parser(sub)
    cli_vaccine.build_parser(sub)
    cli_vaccine.build_report_parser(sub)

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
    elif args.command == "personalize":
        cli_personalize.handle(args)
    elif args.command == "hits":
        cli_hits.handle(args)
    elif args.command == "panel":
        cli_spanning.handle(args)
    elif args.command == "vaccine":
        cli_vaccine.handle(args)
    elif args.command == "vaccine-report":
        cli_vaccine.handle_report(args)


if __name__ == "__main__":
    main()
