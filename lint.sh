#!/usr/bin/env bash

set -e

SOURCES="tsarina tests scripts/develop.py scripts/analyze_vaccine_tradeoffs.py scripts/regenerate_gene_ms_safety_evidence.py examples/cta_vaccine_demo.py"

echo "Running ruff check..."
ruff check $SOURCES

echo "Running ruff format check..."
ruff format --check $SOURCES

echo "All checks passed!"
