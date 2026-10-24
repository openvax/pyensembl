# 2.32.3

Remove the obsolete `typechecks` runtime dependency (#456). Keep built-in
string/integer validation, including boolean/integer subclasses, TypeError for
invalid types and ValueError for empty column names. Validate chromosome types
before memoization so a cached integer cannot admit a float or NumPy integer
(#497). Add initial type annotations, a `py.typed` marker and focused mypy checks
to lint/CI. Existing on-disk caches and indexes remain usable without changes.

# 2.32.2

Require datacache 1.29.3 so explicit `install --overwrite` can recover a corrupt
SQLite annotation index (#495). Recovery preserves the database inode, hard
links, ownership and permissions; a failed rebuild leaves the old database
intact. Normal installation and queries retain existing cache behavior. No
reinstall, migration or automatic repair is required.

# 2.32.1

Finish #387's CLI error handling and fix #494. Missing annotation or reference
DNA files now produce a concise runtime error (exit 1), while invalid arguments
remain argparse errors (exit 2). Explicitly empty GTF/FASTA paths, mirrors and
server URLs are rejected before dataset selection. Malformed data, checksum
failures and unexpected internal errors retain their tracebacks. Existing
cached data, source paths and indexing behavior are unchanged.

# 2.32.0

- Finish the unified numbered/dated release catalogue with scoped, reproducible
  full-GTF comparison evidence, explicit `latest` selection, and filtered
  `available` output. Serialization pins the concrete selection; existing
  defaults and cache paths are unchanged.
- Verify new canonical dated GTF/cDNA/peptide downloads against upstream MD5s.
  Add explicit read-only `verify_sources()` / `pyensembl verify`, preserve
  checksum pins through serialization, and warn for known subset annotations.
  Existing caches need no new metadata, migration, reinstall, or reindex.
- Let reference-name lookup and installed listings reuse known explicit-accession
  dated cache layouts, while keeping installed numbered releases preferred.

# 2.31.2

Include required test helpers and fixtures in source distributions (#489).
Exclude developer-local data and generated cache files. CI and deployment
now extract the source archive and run offline tests before publication.
Wheel contents and installed genome caches are unchanged.

# 2.31.1

Use the species' main (nuclear) genetic code for arbitrary scaffold names and
reference lookups without a contig (#490). This corrects the unnecessarily
unresolved result in 2.31.0. Preserve nonstandard nuclear codes, recognized
mitochondrial/plastid assignments, and explicit source/reference overrides.
Unknown species or explicitly unset tables remain unresolved. Existing caches
and serialized metadata stay usable in place without migration or reindexing.

# 2.31.0

Expose offline reference and transcript genetic-code metadata (#481).
`Genome.genetic_code()` returns a compartment-specific NCBI table assignment
with provenance or an explicit unknown reason. `Transcript.translation_table`
uses annotation declarations before reference/taxonomy defaults. Retain taxon
and species-ancestor IDs, including metadata returned by the new Ensembl API.
One `reference_metadata` argument supports custom contigs and assignments.

Pin metadata through JSON, pickle and mirrors. Existing source paths, cached
files and indexes remain unchanged: old indexes can read source declarations
from the installed GTF into memory without downloads, migration or reindexing.

# 2.30.1

Support tomato SL4.0 in Ensembl Genomes 63 (#453). Use its renamed source
directory and published cDNA/peptide files without requesting a nonexistent
ncRNA FASTA. Missing required files still fail installation. Existing SL3.0
numbered and dated selections keep their URLs, cache paths, and indexes.

# 2.30.0

Remove the unused `EnsemblAnnotation` class, module, and package export. Use
`EnsemblRelease(date, assembly_accession=accession, **options)` for the same
explicit dated selection. The removed class's imports and serialized objects
are no longer supported. Existing cache files, directory layouts, and indexes
remain reusable in place; no reinstall or migration is needed.

Exercise numbered, dated, inferred and explicit assembly choices by installing
tiny GTF/FASTA datasets, querying the resulting indexes, and reopening the
caches offline. Cover human, dog breeds, providers, coverage, optional DNA,
custom labels/directories and serialization.

Fix reopening previously installed private DNA attachments rewriting the
current attachment manifest (#485). Reuse existing indexed/materialized
sources without writes, and track new download-only local attachments
individually so switching between them also preserves the manifest.

# 2.28.3

Share new reference DNA downloads and indexes across dated Ensembl releases
when accession, provider, masking, upstream MD5 and compressed size agree
(#454). Verify downloaded bytes against the manifest before publication.
Missing metadata keeps URLs separate. Existing private dated DNA stays usable
in place, without migration, reindexing or metadata rewrites; numbered-release
shared caches retain their paths. Both EnsemblRelease and EnsemblAnnotation
benefit. Clarify that ordinary dated releases need only EnsemblRelease;
EnsemblAnnotation is for explicit assembly, provider and coverage selection.

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
