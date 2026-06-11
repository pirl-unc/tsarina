#!/usr/bin/env bash

set -e

VERSION=$(python -c "from tsarina.version import __version__; print(__version__)")
echo "Deploying tsarina v${VERSION}"
echo ""

echo "==> Running lint checks..."
./lint.sh

echo ""
echo "==> Running tests..."
./test.sh

echo ""
echo "==> Cleaning old builds..."
# Build hermetically: a stale egg-info/SOURCES.txt or build/ tree can re-leak
# untracked files (e.g. the gitignored data/*.regen.csv sidecar) into the sdist.
rm -rf dist/* build/ ./*.egg-info

echo ""
echo "==> Building distribution..."
python -m build

echo ""
echo "==> Built artifacts:"
ls -lh dist/

echo ""
echo "==> Uploading to PyPI..."
twine upload dist/*

echo ""
echo "Deploy complete! https://pypi.org/project/tsarina/${VERSION}/"
