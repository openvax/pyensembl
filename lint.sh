#!/bin/bash
set -o errexit

ruff check pyensembl/ tests/ \
&& \
echo "Passes ruff check"
