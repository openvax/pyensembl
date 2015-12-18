#!/bin/bash
set -o errexit


# getting false positives due to this issue with pylint:
# https://bitbucket.org/logilab/pylint/issues/701/false-positives-with-not-an-iterable-and

find . -name '*.py' \
  | xargs pylint \
  --errors-only \
  --disable=print-statement,unsubscriptable-object,not-an-iterable

echo 'Passes pylint check'
