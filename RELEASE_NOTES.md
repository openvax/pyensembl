# 2.28.2

Fix local transcript/protein FASTA index collisions (#467). Newly created
indexes identify the original source path and file contents; new cached copies
also identify their original paths. Changed local inputs select fresh indexes,
including in loaded readers, without deleting older generations. Unchanged
gzip copies with decompression enabled are reused without copying again. Existing
basename-only sources and indexes stay usable in place, with no reinstall or
automatic migration; they retain their historical source assumptions.

# 2.28.0

Add `EnsemblRelease(...).with_mirror(URL)` to reopen the same source-specific
Genome installed by `pyensembl install --custom-mirror URL` (#470). The CLI
uses this shared constructor. Mirror downloads and indexes retain their own
identity; the default EnsemblRelease constructor continues to select canonical
sources. See the cache guide for installation, Python/Varcode use and migration.

# 2.27.8

Restore in-place reuse of existing basename-only remote caches and their
SQLite/FASTA indexes after the 2.27.5 compatibility regression. Existing
installations need no reinstall, migration, renaming, or reindexing. New
downloads retain URL-specific identity; legacy caches retain their historical
source assumptions. Explicit refreshes leave legacy files intact.

# 2.27.2

Transcript.protein_id, Transcript.protein and Transcript.protein_sequence now
return None when a noncoding-only or custom annotation has no CDS protein_id
column. This fixes the missing-table error reported in #459 without changing
canonical protein identities in ordinary annotations.
