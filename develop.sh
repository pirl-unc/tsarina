#!/usr/bin/env bash

set -e

# Install into the virtualenv that is already active, if there is one.
#
# This script used to create and activate ./.venv unconditionally.  Run from a
# shell that already had a virtualenv active, it therefore installed into a
# *different* environment than the developer was using -- reporting success
# while `tsarina` on PATH went on resolving to whatever stale copy was in the
# active environment.  That is how a 1.23.1 install survived two releases
# behind this checkout, answering `--serotype` queries with retired code.
# Pattern borrowed from hitlist's own develop.sh (pirl-unc/hitlist#400).
if [ -n "$VIRTUAL_ENV" ]; then
    echo "Installing into the active virtualenv: $VIRTUAL_ENV"
else
    VENV_DIR=".venv"
    if [ ! -d "$VENV_DIR" ]; then
        echo "Creating virtual environment at $VENV_DIR..."
        python -m venv "$VENV_DIR"
    fi
    # shellcheck disable=SC1091
    source "$VENV_DIR/bin/activate"
fi

if command -v uv &> /dev/null; then
    echo "Using uv to install package with development dependencies..."
    uv pip install -e ".[dev]"
    PIP_INSTALL=(uv pip install)
else
    echo "uv not found, falling back to regular pip..."
    pip install -e ".[dev]"
    PIP_INSTALL=(pip install)
fi

# Develop against sibling checkouts when they are present, so tsarina tracks
# the code its results actually depend on rather than the last PyPI release.
# Falls back silently to the pinned wheel for any repo that is absent.
#
# These are the packages that decide what tsarina reports: the MS evidence
# index and its curation (hitlist), the CTA definitions (oncoref), the MHC
# allele and serotype vocabulary (mhcgnomes), the gene models (pyensembl), and
# the download / parsing layer all of the above sit on.  A stale copy of any
# of them changes scientific output without changing tsarina, which is the
# failure this list exists to prevent: a mhcgnomes four releases old silently
# drops whole serological specificities from the evidence index.
#
# Prediction backends (mhcflurry, mhctools, topiary) are deliberately absent --
# they are optional extras, and pinning their dev trees here would make a
# scoring change look like a tsarina change.
SIBLINGS=(hitlist oncoref mhcgnomes pyensembl datacache gtfparse serializable sercol)
SIBLING_ROOT="${SIBLING_ROOT:-..}"
for sibling in "${SIBLINGS[@]}"; do
    sibling_dir="$SIBLING_ROOT/$sibling"
    if [ -d "$sibling_dir" ]; then
        echo "Installing sibling $sibling editable from $sibling_dir ..."
        # --no-deps is load-bearing.  Resolving one sibling's dependencies can
        # replace another sibling's editable install with a released wheel:
        # installing sercol pulled serializable==0.4.1 over the editable
        # serializable 1.1.0 checkout, silently undoing the line above it.
        # tsarina's own ".[dev]" install above resolves the closure; a sibling
        # that needs a dependency its last release did not is the one case
        # wanting a manual install, and the report below names it.
        "${PIP_INSTALL[@]}" -e "$sibling_dir" --no-deps
    fi
done

# Say where everything landed.  The failure being guarded against is not an
# install error -- it is an install that succeeds into the wrong place, so the
# only useful confirmation is the resolved path, not an exit code.
echo
python - <<'PY'
import importlib
import os
import shutil
import subprocess
from pathlib import Path

import tsarina

print(f"import tsarina -> {tsarina.__version__}  ({tsarina.__file__})")

cli = shutil.which("tsarina")
if cli is None:
    print("WARNING: no `tsarina` on PATH")
else:
    reported = subprocess.run([cli, "--version"], capture_output=True, text=True).stdout.strip()
    print(f"`tsarina` on PATH -> {cli}")
    print(f"                     {reported}")
    if tsarina.__version__ not in reported:
        print(
            f"\nWARNING: the console script reports {reported!r} but the importable "
            f"package is {tsarina.__version__}.\n"
            "Another install is shadowing this one; `pip uninstall tsarina` until "
            "none remain, then re-run this script."
        )

SIBLING_NAMES = "hitlist oncoref mhcgnomes pyensembl datacache gtfparse serializable sercol"

print()
print("Resolved dependencies whose code decides tsarina's output:")
root = Path(os.environ.get("SIBLING_ROOT", "..")).resolve()
shadowed = []
for name in SIBLING_NAMES.split():
    checkout = root / name
    try:
        module = importlib.import_module(name)
    except ImportError:
        print(f"  {name:12s} NOT INSTALLED")
        continue
    version = getattr(module, "__version__", "?")
    path = Path(module.__file__ or "")
    live = not checkout.is_dir() or checkout in path.parents
    print(f"  {name:12s} {version:12s} {path}")
    if not live:
        shadowed.append(name)

if shadowed:
    print(
        f"\nWARNING: {', '.join(shadowed)} resolved outside their checkouts despite "
        "one being present.\nA released wheel is shadowing the dev tree — most likely "
        "pulled in as another\npackage's dependency. Re-run this script, and if it "
        "persists, install that\none last with `--no-deps`."
    )
PY
