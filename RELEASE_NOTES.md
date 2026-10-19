# 2.27.6

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
