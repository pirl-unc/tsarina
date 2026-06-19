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

"""Shared CLI argument helpers, so the ``hits``/``personalize``/``panel``
subcommands spell common options identically (one source of truth instead of
three drifting copies)."""

from __future__ import annotations

import argparse
import sys

#: mhctools predictors accepted by every scoring-capable subcommand.
SUPPORTED_PREDICTORS = ("mhcflurry", "netmhcpan", "netmhcpan_el")


def split_csv(value: str) -> list[str]:
    """Parse a comma-separated argument into a trimmed, non-empty list."""
    return [s.strip() for s in value.split(",") if s.strip()]


def parse_lengths(value: str) -> tuple[int, ...]:
    """argparse type for ``--lengths`` (comma-separated integers)."""
    try:
        return tuple(int(x) for x in split_csv(value))
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"--lengths must be integers: {value!r}") from e


def add_iedb_cedar_args(parser: argparse.ArgumentParser) -> None:
    """Add the IEDB/CEDAR path overrides, spelled the same on every subcommand."""
    parser.add_argument(
        "--iedb",
        dest="iedb_path",
        default=None,
        help="Explicit IEDB ligand CSV path (default: hitlist data registry).",
    )
    parser.add_argument(
        "--cedar",
        dest="cedar_path",
        default=None,
        help="Explicit CEDAR ligand CSV path (default: hitlist data registry).",
    )


def add_predictor_arg(parser: argparse.ArgumentParser, *, context: str) -> None:
    """Add the ``--predictor`` choice with a context-specific help string."""
    parser.add_argument(
        "--predictor",
        choices=SUPPORTED_PREDICTORS,
        default="mhcflurry",
        help=f"mhctools predictor for {context}. Defaults to mhcflurry.",
    )


def warn_ignored_flags(mode_flag: str, ignored: list[str]) -> None:
    """Print a stderr warning that *mode_flag* makes *ignored* flags no-ops."""
    if ignored:
        joined = ", ".join(sorted(ignored))
        print(
            f"Warning: {mode_flag} is a mode switch; these flags are ignored: {joined}",
            file=sys.stderr,
        )
