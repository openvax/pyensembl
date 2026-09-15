#!/usr/bin/env bash
set -eo pipefail

./lint.sh
exec ./test.sh "$@"
