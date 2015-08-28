#!/bin/bash
set -o errexit

find . -name '*.py' \
  | xargs pylint \
  --errors-only \
  --disable=print-statement
  --ignored-classes=nose.tools

echo 'Passes pylint check'
