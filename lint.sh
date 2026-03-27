#!/usr/bin/env bash

set -e

SOURCES="perseus tests"

echo "Running ruff check..."
ruff check $SOURCES

echo "Running ruff format check..."
ruff format --check $SOURCES

echo "All checks passed!"
