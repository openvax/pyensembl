# Install and manage data

PyEnsembl downloads annotation and sequence files once and indexes them in a
local cache. After that, queries need no network access.

## Install data

Install from the command line:

```sh
pyensembl install --release 93 --species human
```

or from Python, for example in a notebook or pipeline:

```python
from pyensembl import EnsemblRelease

data = EnsemblRelease(93, species="human")
data.download()
data.index()
```

Both skip files that are already downloaded or indexed, so rerunning finishes
an interrupted install. `install` prints one progress line per step on stderr,
with progress bars in a terminal; add `--verbose` (`-v`) to see every download
and database step. In Python, pass `show_progress=True` to `download()` and
`index()`. Use `--overwrite` or `overwrite=True` to replace existing files.

[Choose a reference](assembly-selection.md) explains release, species and
assembly options. Whole-genome DNA is a [separate download](reference-dna.md).

### Recover a corrupt annotation index

If SQLite reports a corrupt annotation database, explicitly rerun your original
install command with `--overwrite`, or call `data.index(overwrite=True)` to
rebuild from the installed GTF. The CLI also redownloads remote source files
because `--overwrite` applies to the entire installation.

PyEnsembl 2.32.2 requires datacache 1.29.3 for this recovery. It builds and
validates a replacement before publishing it into the existing database,
preserving its inode, hard links, ownership and permissions. A failed rebuild
leaves the original database intact. Recovery needs permission to write the
database, create a child process and use SQLite's reset API. If the runtime
cannot support recovery, the operation fails while retaining the old database.

Ordinary installs and queries do not repair or migrate caches. Existing healthy
installations continue to work in place without reinstalling or reindexing.

## Check what is installed <a id="list-installed-genomes"></a>

```sh
pyensembl list
```

```text
Species  Assembly  Release  Annotation   Reference DNA      Location
human    GRCh38    81       indexed      toplevel, indexed  ~/Library/Caches/pyensembl/GRCh38/ensembl81
human    GRCh38    82       not indexed  -                  ~/Library/Caches/pyensembl/GRCh38/ensembl82
custom   GRCm38    mine1    indexed      -                  ~/Library/Caches/pyensembl/GRCm38/mine1
```

The Annotation column shows whether a release is ready to query:

| Annotation | Meaning |
| --- | --- |
| `indexed` | Everything is downloaded and indexed; queries need no network access or setup |
| `not indexed` | Files are downloaded, but the first query would spend minutes indexing them |
| `incomplete` | Some downloads are missing |
| `invalid` | A source is empty or not a regular file |
| `inaccessible` | A source cannot be read |

Run `pyensembl install` for that release to finish. Add `--species` for
non-human genomes; custom genomes need their original install options.

In Python, `installed()` is `True` when everything a genome is configured with,
reference DNA included, is downloaded and indexed. It only reads the cache.
`genome_for_reference_name` picks the newest installed release of an assembly,
else the newest downloaded one, else the newest supported Ensembl release:

```python
from pyensembl import EnsemblRelease, genome_for_reference_name
EnsemblRelease(93).installed()
genome_for_reference_name("GRCh38")

# Ensembl releases with any files in the cache, ready or not:
from pyensembl.shell import collect_all_installed_ensembl_releases
collect_all_installed_ensembl_releases()
```

## Cache location

PyEnsembl keeps all of its data under one directory, with a subdirectory per
genome (`<reference>/<annotation><version>`, e.g. `GRCh38/ensembl81`). By
default this is the platform cache directory that datacache chooses:
`~/.cache/pyensembl` on Linux, `~/Library/Caches/pyensembl` on macOS, and
`%LOCALAPPDATA%\pyensembl\pyensembl\Cache` on Windows. Releases that
PyEnsembl 2.16 or earlier installed on Windows stay in their old per-genome
directories and keep working. To use another location, set
`PYENSEMBL_CACHE_DIR`; the data then goes in its `pyensembl` subdirectory:

```sh
export PYENSEMBL_CACHE_DIR=/custom/cache/dir
```

In Python, set `os.environ["PYENSEMBL_CACHE_DIR"]` before creating a genome.

New annotation and transcript/protein FASTA downloads have a SHA-256
prefix identifying the complete source URL and decompression mode, followed by
the remote filename (for example, `<key>-annotation.gtf.gz`). Their SQLite and
FASTA indexes retain that prefix. URLs with different hosts, paths or query
parameters use separate new files, even in the same genome directory. Changing a
signed URL also selects a new cache entry. Query parameters and credentials
are included in the hash but are not copied into the filename.

**Existing installations remain usable in place.** Version 2.27.8 restores
reuse of basename-only downloads and their original indexes, correcting the
2.27.5 regression. No reinstall, migration, renaming, provenance rewrite or
reindexing is required, including for read-only shared caches. Files already
downloaded with URL-specific names also remain usable and take precedence
when both layouts exist. The shared reference-DNA store keeps its existing
behavior.

Legacy entries retain their historical basename lookup and source assumptions,
including installations from mirrors. Their original full URL cannot always
be recovered from absent or redacted provenance. The new URL-identity guarantee
applies to the new layout. If you explicitly need to acquire a different URL
despite an existing legacy basename, use `download(overwrite=True)` or
`pyensembl install --overwrite`: the download uses a URL-specific destination
and leaves the legacy files intact. New `--custom-mirror` downloads belong to
that mirror's URLs; use the same mirror to access them through the CLI.

## Local transcript and protein FASTA indexes

Starting with 2.28.2, new local FASTA indexes identify both the original absolute
source path and a SHA-256 of the indexed bytes. Two local files named
`sequences.fa` can safely use the same cache directory. Changed bytes select
a fresh index even if the file size and modification time stay the same;
older index generations are retained until explicit deletion. New cached
copies also identify their original paths and remain usable when the original
file disappears. Call `download()` or `index()` to refresh a copied input.

Existing basename-only local copies and pickle indexes continue to work in
place, including in read-only caches. No migration or reindexing is required.
Their original source identity cannot be recovered reliably, so they retain
the historical basename assumptions; the new isolation guarantee applies to
new entries. Explicitly deleting an index and reindexing adopts the new layout.
Refreshing a changed legacy copy writes a new destination and leaves the
legacy files intact.

## Reuse a mirror installation

For a flat directory containing the original Ensembl filenames, install with
an explicit mirror:

```sh
pyensembl install --release 75 --species human --custom-mirror https://github.com/openvax/ensembl-data/releases/download/GRCh37.75/
```

In later Python processes, select that same source with `with_mirror`:

```python
from pyensembl import EnsemblRelease

mirror = "https://github.com/openvax/ensembl-data/releases/download/GRCh37.75/"
data = EnsemblRelease(75, species="human").with_mirror(mirror)
assert data.installed()
print(data.gene_by_id("ENSG00000141510").name)
```

`with_mirror` returns a `Genome` whose GTF, transcript and protein URLs point
at the mirror. Pass this object explicitly to downstream tools; for example,
`Variant("1", 53513530, "A", "C", genome=data)` in Varcode. To reuse the data
through the CLI, repeat `--custom-mirror` for the same release and species.
Requested Ensembl reference DNA is also read from the mirror; explicitly
attached local DNA is preserved.

New mirror installations do **not** populate the source-specific cache used by
`EnsemblRelease(75)` or Varcode's `genome=75`. If those default constructors are
needed, install canonical sources with `pyensembl install --release 75 --species
human`. Existing URL-specific mirror installs from 2.27.5 onward need no
redownload or reindexing when opened with the same mirror. The legacy
basename-cache compatibility behavior described above still applies; this API
does not copy indexes or change existing files' source identity.

## Free disk space <a id="managing-disk-space"></a>

```sh
pyensembl delete-index-files --release 93  # keep downloads, reindex on use
pyensembl delete-all-files --release 93    # all of release 93's files
pyensembl prune --dry-run                  # list unused shared DNA
pyensembl prune
```

Add `--species` for non-human releases.

Compatible releases share one copy of Ensembl DNA, so deleting a release keeps
DNA that other releases still use, and `prune` removes DNA that no release
references. `prune` skips DNA that is being downloaded or indexed, never
touches local FASTA files, and deletes nothing if any release's DNA metadata is
malformed (`pyensembl list` shows which one).

`delete-index-files` keeps shared DNA indexes because other releases may use
them; rebuild one with `index_genome_fasta(overwrite=True)`. In Python,
`prune_genome_fastas(dry_run=True)` returns `(path, bytes)` candidates, and
`pyensembl list --check-genome-fasta` verifies each release's DNA index.

## Inspect files without installing <a id="inspect-data-without-installing"></a>

Inspect a selected genome's source files and indexes, including paths,
availability, sizes and any recorded download provenance:

```sh
pyensembl inspect --release 93
pyensembl inspect --release 93 --with-genome-fasta --json
# Custom genomes use the same --gtf / --transcript-fasta / --protein-fasta
# and --reference-name / --annotation-name options as install.
```

In Python, `EnsemblRelease(93).inspect_data()` returns annotation and optional
reference-DNA readiness, `installed`, and a `files` dictionary keyed by role
(for example `gtf`, `gtf_index`, `transcript_fasta_1`). Values are datacache
`FileInspection` objects with `path`, `status`, `error`, `size`, `mtime`,
`source_url`, `fetched_at`, `recorded_sha256` and `verified` fields. CLI JSON
is an array of reports, with errors rendered as strings. Add
`--check-genome-fasta` (Python: `check_genome_fasta=True`) to validate existing
DNA indexes more thoroughly.

Inspection only reads: no network, copying, directory creation, indexing or
pickle deserialization. It checks that source files are readable, regular and
nonempty and that SQLite indexes are complete; it does not validate biological
contents. Recorded provenance is advisory, not a trusted checksum, and files
from older versions have none. Downloads are published atomically, so a failed
overwrite keeps the previous file; replace an invalid file with
`pyensembl install --overwrite`.

## Share a cache

Reads take no locks and write nothing, so a fully installed and indexed cache
can be read-only for other users. A download or index build locks only the
file it writes; registering, deleting and pruning releases briefly lock the
whole cache.

New files and lock files follow your umask, so set
`umask 002` (or default ACLs) before installing into a group-shared cache.
Files downloaded before PyEnsembl 2.13.1 were readable only by their owner;
share them with `chmod -R g+rX "$PYENSEMBL_CACHE_DIR/pyensembl"` (or the
platform cache directory). `dna_cache` may be a symlink, e.g. to a larger disk.

## How reference DNA is stored <a id="how-the-shared-dna-cache-works"></a>

Ensembl DNA is stored once per upstream file under `pyensembl/dna_cache/`:

```text
pyensembl/dna_cache/
  homo_sapiens/ftp.ensembl.org/GRCh38-GCA_000001405.18/
    toplevel/unmasked/fasta/<file key>/
      sequence.fa        uncompressed, even when downloaded as .fa.gz
      sequence.fa.fai
      object.json        full identity of the upstream file
      index.json
```

For numbered releases, before downloading PyEnsembl reads Ensembl's small README and CHECKSUMS files
to see whether another release already has the same file. The versioned
assembly accession distinguishes assembly patches, and the 16-character file
key (a SHA-256 prefix of the assembly, Ensembl's checksum and compressed size)
distinguishes upstream revisions. These are metadata checks: Ensembl's Unix
checksums are not cryptographic hashes. If the metadata is incomplete, the
assembly directory ends in `-unverified` and each release keeps its own copy.
Local FASTA files and custom mirrors are never shared.

For dated releases, new installs read the official `genome/md5sum.txt` and
compressed size. Matching versioned accession, provider, masking, MD5 and size
share one download and index across dates. Downloaded compressed bytes must
match the MD5 before the FASTA is published. The date is excluded from a
verified file identity, so a new annotation date does not duplicate unchanged
DNA. Changed checksums select a separate object; missing or ambiguous metadata
keeps each source URL separate. These objects live under
`dna_cache/ensemblorganisms/ftp.ebi.ac.uk/<accession>/<provider>/toplevel/<mask>/fasta/<file key>/`.

**Existing dated DNA installations remain usable in place.** Their private
files, indexes and manifests are neither migrated nor rewritten. Repeated
`download_genome_fasta()` and `index_genome_fasta()` calls reuse them, including
offline and read-only caches. New sharing therefore prevents future duplicates;
it does not automatically reclaim existing copies. Numbered and dated DNA
are not assumed interchangeable solely because their accessions match: their
file formats and upstream checksum schemes differ.

Downloads retry transient HTTP failures and stalls of five minutes, and are
checked against the upstream size. An interrupted DNA download resumes on the
next install on POSIX systems, using only bytes the server confirms come from
the same file; on Windows it starts over.
