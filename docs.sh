#!/usr/bin/env bash
set -euo pipefail
python -m mkdocs build --strict
python scripts/check_docs_links.py
