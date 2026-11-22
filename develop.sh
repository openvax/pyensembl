#!/usr/bin/env bash
set -euo pipefail

"${PYTHON:-python}" -m pip install -e "." "$@"
