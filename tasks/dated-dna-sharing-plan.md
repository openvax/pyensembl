# Dated reference DNA sharing (#454)

## Findings and contract

Official dated DNA has a genome/md5sum.txt manifest (not the dataset root).
Verified on 2026-10-19: human GRCh38 2025_12 and 2026_04 publish the same
three compressed DNA hashes; unmasked MD5 30d2f26286e2053454a74c723e70e701,
847651000 bytes. Accession alone is insufficient: masking, provider and
upstream file revisions must remain distinct. Numbered and dated formats
have different metadata and are not assumed interchangeable.

Existing release-private dated files, indexes and manifests must stay usable
in place, including read-only/offline reuse and repeated download/index calls.
No automatic migration, copy, rename, reindex, deletion or metadata rewrite.
Existing numbered shared identities and paths remain unchanged.

## Implementation

- Recognize only official HTTPS ftp.ebi.ac.uk ensemblorganisms URLs with
  versioned GCA/GCF accession, provider, date and standard masked DNA filenames.
- For new dated installs, share by accession + provider + mask + MD5 + exact
  compressed size. Fetch small metadata only during explicit acquisition;
  verify downloaded bytes against MD5 before decompression/publication.
  Missing, malformed or conflicting metadata isolates by full URL instead.
- Reuse existing SharedGenomeFasta locks, atomic writes, references, indexes
  and conservative pruning. Validate metadata-derived directory components.
- Prefer an existing dated release-private materialized file if it has no
  shared reference, preserving its current directory and implementation.
  New installs through EnsemblRelease and EnsemblAnnotation share equally.
- Clarify the API: ordinary date selection needs only EnsemblRelease;
  EnsemblAnnotation exists for explicit older/unregistered accessions,
  providers and annotation coverage. Retain compatible constructors.
- Patch version bump, release notes and cache/platform documentation.

## Verification and delivery

Use tiny local DNA with simulated official metadata to cover equal dates,
changed hashes/size, masks, providers, accession versions/GCF, malformed or
missing manifests, download checksum mismatch, offline read-only reuse,
legacy-private snapshots, concurrency, pruning and both constructors. Run
focused tests, ./lint.sh, ./test.sh, ./docs.sh, pip check and saved legacy
snapshot verification. Open a separate PR closing #454, pass CI, merge and
run ./deploy.sh from clean main; verify published artifacts. Never touch
actual user caches or the shared environment.
