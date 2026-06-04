#!/usr/bin/env bash

set -e

VENV_DIR=".venv"

if [ ! -d "$VENV_DIR" ]; then
    echo "Creating virtual environment at $VENV_DIR..."
    python -m venv "$VENV_DIR"
fi

# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

# Check if UV is installed and available in the PATH
if command -v uv &> /dev/null; then
    echo "Using uv to install package with development dependencies..."
    uv pip install -e ".[dev]"
    PIP_INSTALL=(uv pip install)
else
    echo "uv not found, falling back to regular pip..."
    pip install -e ".[dev]"
    PIP_INSTALL=(pip install)
fi

# Develop against the sibling hitlist checkout when present, so tsarina always
# tracks the latest hitlist code (its MS-evidence schema changes frequently).
# Falls back silently to the PyPI-pinned hitlist if the sibling is absent.
HITLIST_DIR="${HITLIST_DIR:-../hitlist}"
if [ -d "$HITLIST_DIR" ]; then
    echo "Installing sibling hitlist editable from $HITLIST_DIR ..."
    "${PIP_INSTALL[@]}" -e "$HITLIST_DIR"
fi
