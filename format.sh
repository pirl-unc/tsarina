#!/usr/bin/env bash

set -e

SOURCES="tsarina tests"

echo "Running ruff format..."
ruff format $SOURCES

echo "Formatting complete!"
