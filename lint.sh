#!/usr/bin/env bash

set -e

SOURCES="tsarina tests scripts/develop.py"

echo "Running ruff check..."
ruff check $SOURCES

echo "Running ruff format check..."
ruff format --check $SOURCES

echo "All checks passed!"
