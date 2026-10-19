# Python 3.11 support floor (#463)

## Specification

Require Python 3.11 or newer without an upper bound. Keep CI coverage through
3.14 and drop 3.9/3.10 jobs. Align classifiers, docs, and repository guidance.
Remove the obsolete pre-3.10 filelock branch and keep the supported filelock
constraint. Preserve genome/cache behavior and other dependency bounds. Bump
the patch version for release.

## Plan

- [x] Read the issue, package metadata, CI, release script, and guidance.
- [x] Check plan before implementation: this is a support-policy change; tests
  must still exercise the reference-backed suite using a writable copy of
  existing Ensembl caches.
- [x] Update metadata, CI, docs, version, and obsolete dependency marker.
- [x] Build isolated Python 3.11 and 3.14 environments; verify dependency
  resolution against released foundation packages.
- [x] Run `./lint.sh`, `./test.sh`, and documented example checks with a
  writable reference cache; review the diff.
- [ ] Open PR, pass CI, merge, deploy from clean main, verify PyPI and tag.

## Review

`./lint.sh` passed. `./test.sh -q` passed 711 tests on Python 3.11.12 and
711 tests on Python 3.14.6 (94% coverage); the three 3.14 resource warnings
come from existing SQLite connections. The first 3.11 run under the restricted
sandbox denied loopback sockets for 15 local HTTP tests; rerunning with local
socket access passed all 711 tests. All 16 documented examples passed on both
versions using a writable copy-on-write cache of the Ensembl reference data.
Both isolated environments resolved typechecks 0.1.2, serializable 1.3.2,
gtfparse 3.0.5, and datacache 1.28.0; dependency checks passed. Installed
metadata reports version 2.28.1 and `Requires-Python: >=3.11`.
