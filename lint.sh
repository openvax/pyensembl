#!/bin/bash
set -o errexit

ruff check pyensembl/ \
&& \
echo "Passes ruff check"
