#!/usr/bin/env bash

set -e

SOURCES="ctabase tests"

echo "Running ruff format..."
ruff format $SOURCES

echo "Formatting complete!"
