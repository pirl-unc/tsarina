#!/usr/bin/env bash

set -e

SOURCES="tsarina tests scripts/develop.py scripts/analyze_vaccine_tradeoffs.py examples/cta_vaccine_demo.py"

echo "Running ruff format..."
ruff format $SOURCES

echo "Formatting complete!"
