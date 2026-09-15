#!/usr/bin/env bash

set -eo pipefail

./lint.sh
./test.sh
python -m pip install --upgrade build
python -m pip install --upgrade twine
rm -rf dist
python -m build
git --version
python -m twine upload dist/*
git tag "$(python pyensembl/version.py)"
git push --tags
