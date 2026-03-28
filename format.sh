#!/usr/bin/env bash

set -e

SOURCES="perseo tests"

echo "Running ruff format..."
ruff format $SOURCES

echo "Formatting complete!"
