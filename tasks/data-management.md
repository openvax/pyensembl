# Data management through released datacache APIs

## Problem and scope

PyEnsembl delegates remote transfers and SQLite publication to datacache, but
still performs its own existence-only source checks and non-atomic local copies.
Those paths disagree: empty sources and directories may count as installed,
copied files require an original source even after import, and new downloads
have no visible provenance. Users cannot inspect the exact files and indexes
behind a selected genome without reaching into implementation details.

Use released datacache 1.16.1 (APIs available at the existing 1.14 floor). Preserve
cache paths, domain exceptions for missing files, index ownership, and shared
DNA identity/reference/pruning. Do not depend on the unrelated in-progress
datacache checkout or pending bundle enhancements.

## Design

1. Give DownloadCache a side-effect-free local-path resolver and inspection
   adapter. Delegate readability, regular-file checks, structured errors and
   provenance to datacache.inspect_file. Add the PyEnsembl requirement that data
   files are nonempty. Keep an explicit empty-files-ok option for the existing
   required_local_files_exist compatibility contract.
2. Delegate cache-hit validation to inspect_file and forced replacement to fetch_file,
   with explicit raw/decompression settings and provenance recording. Missing
   offline data retains the existing install hints; corrupt/inaccessible data
   raises its structured cause and requires explicit overwrite to repair.
3. Import local files through the same atomic fetch_file path using a file URI.
   Reuse a cached imported file offline even after the original is removed.
   Apply requested decompression consistently to copied and remote inputs;
   preserve attached sources and permit copying a file already at its target.
4. Centralize annotation source/index inventories in Genome. Route readiness,
   required-file checks and release selection through readable, nonempty file
   inspection. Keep the complete SQLite/schema check. Distinguish missing,
   incomplete, invalid and inaccessible sources; incomplete indexes remain
   'not indexed'. Never parse FASTA dictionaries just to inspect readiness.
5. Add Genome.inspect_data(): annotation/DNA readiness plus a role-keyed mapping
   of datacache FileInspection objects for configured sources and indexes.
   Represent a database with wrong completion metadata as an invalid index.
   DNA index status still uses its source fingerprint and existing DNA contract.
6. Add pyensembl inspect using existing genome selection flags, with a readable
   file table and --json for scripts. JSON serializes structured filesystem
   errors to messages and includes sizes, timestamps and recorded provenance.
   Inspection is offline and never creates/copies/downloads/indexes/repairs.
   Existing list output stays compact and gains accurate source states.
7. Update usage docs and bump 2.19.0 to 2.20.0. File discovered defects as an
   upstream issue and link it from the PR.

## Verification and delivery

- Reproduce existing empty/directory-source readiness, non-atomic copy,
  compressed-local-copy and deleted-original reuse failures with small fixtures.
- Verify real datacache gzip/plain acquisition, atomic preservation after a
  failed overwrite, read-only/offline reuse, provenance and stale/malformed
  receipts; prohibit writes/network during inspection tests.
- Exercise installed/list/reference selection against invalid/inaccessible
  sources and incomplete indexes; cover custom, Ensembl and reference-DNA data.
- Exercise human/JSON CLI inspection, absence, legacy files and selection flags.
- Run ./lint.sh and the full ./test.sh. Use an isolated environment with the
  released datacache wheel and isolated test outputs; keep real Ensembl fixtures
  available without modifying the user's cache. Check wheel/sdist artifacts and
  review the diff before submitting.
- Open a PR, wait for required CI, merge, deploy with ./deploy.sh from a clean
  main checkout, verify PyPI, then review the next dependency/urgency group.

## Verification results

- Filed the acquisition/readiness defects as openvax/pyensembl#425.
- The four initial regressions fail on 2.19.0 and pass after the changes.
- Released datacache 1.16.1: full ./test.sh passes (492 passed, 1 skipped;
  92% coverage), and ./lint.sh passes. HTTP-resume tests require an unsandboxed
  localhost bind; fixture data lives in an isolated copy-on-write cache.
- Minimum supported datacache 1.14.0: all 29 focused data-management/readiness
  tests pass, so the dependency floor need not change.
- Equality fixtures now use a nonempty annotation header. The cache-location
  test restores the caller's environment instead of deleting it (related #422).
