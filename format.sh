#!/usr/bin/env bash

set -e

SOURCES="tsarina tests scripts/develop.py"

echo "Running ruff format..."
ruff format $SOURCES

echo "Formatting complete!"
