# Local FASTA identity (#467)

## Contract

New indexes must belong to the configured local source and its actual bytes.
Different source paths, changed bytes (including equal-length changes with a
restored mtime), transcript/protein roles, and copied sources must not reuse
another file's sequence dictionary. Existing legacy sources and indexes remain
usable in place, including read-only caches; no migration, metadata rewrite,
renaming, reinstall, or automatic deletion. Unknown legacy provenance retains
the historical assumptions, as in the remote-cache compatibility correction.

## Implementation

- Centralize FASTA index selection for SequenceData and Genome inspection.
  For new local entries, use a namespace derived from the original absolute
  source path and a SHA-256 of the actual indexed file. Keep a readable filename.
  Memoize byte hashes only within a process, keyed by complete stat identity;
  detect changes during hashing/parsing before publishing any index.
- Reuse a legacy basename index unchanged if no source-specific generation has
  been created for that source. New remote indexes keep their existing download
  filename identity, and existing remote caches remain unchanged.
- Refresh loaded local sequence dictionaries when their selected source
  generation changes. Permission-only changes must reuse identical bytes/indexes.
- Give new copied local sources their own original-path/mode-specific download
  destination. Reuse existing legacy copies in place; if an explicit refresh or
  changed local input needs new bytes, write a new destination and retain the
  old legacy files. Copies remain usable after the original disappears.
- Index readiness, inspection and deletion must agree with sequence lookup.
  Explicit deletion can find this source's generations even when the original
  has disappeared, without deleting another source's indexes or any FASTA.
- Document the new identity and the legacy compatibility boundary; patch bump.

## Verification

Reproduce the published 2.28.1 collision before the fix. Cover transcript and
protein collisions, separate instances and offline reuse, changed content with
unchanged size/mtime, loaded readers, compressed/plain and copied/attached
sources, missing originals, deletion isolation, interrupted/changed inputs,
long names and read-only legacy/current caches. Snapshot legacy bytes, names,
permissions, inodes and timestamps, forbidding acquisition/rebuild on reuse.
Keep full legacy fixtures independent of the new path selection code.

Run focused tests serially, ./lint.sh, ./test.sh and pip check in an isolated
repository environment with temporary data. Open a PR closing #467; merge only
after review/checks, then ./deploy.sh from clean main, verify PyPI artifacts.
Never change the shared environment or actual user caches.
