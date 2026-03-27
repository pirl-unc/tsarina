#!/usr/bin/env bash

set -e

echo "Running lint checks..."
./lint.sh

echo "Running tests..."
./test.sh

echo "Building distribution..."
python -m build

echo "Uploading to PyPI..."
twine upload dist/*

echo "Deploy complete!"
