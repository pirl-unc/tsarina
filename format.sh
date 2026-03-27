#!/usr/bin/env bash

set -e

SOURCES="perseus tests"

echo "Running ruff format..."
ruff format $SOURCES

echo "Formatting complete!"
