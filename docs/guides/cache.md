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

**Annotation** is `indexed` when everything is downloaded and indexed, so
queries need no network access or setup. `not indexed` means the files are
downloaded but the first query would spend minutes indexing them, and
`incomplete` means some downloads are missing. `invalid` means a source is
empty or not a regular file; `inaccessible` means it cannot be read. Run
`pyensembl install` for that release to finish (add `--species` for non-human
genomes; custom genomes need their original install options).

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

## Free disk space <a id="managing-disk-space"></a>

```sh
pyensembl delete-index-files --release 93  # keep downloads, reindex on use
pyensembl delete-all-files --release 93    # all of release 93's files
pyensembl prune --dry-run                  # list unused shared DNA
pyensembl prune
```

Add `--species` for non-human releases. Compatible releases share one copy of
Ensembl DNA, so deleting a release keeps DNA that other releases still use, and
`prune` removes DNA that no release references. It skips DNA that is being
downloaded or indexed, never touches local FASTA files, and deletes nothing if
any release's DNA metadata is malformed (`pyensembl list` shows which one).
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

New files follow your umask, as do lock files on Python 3.10+, so set
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

Before downloading, PyEnsembl reads Ensembl's small README and CHECKSUMS files
to see whether another release already has the same file. The versioned
assembly accession distinguishes assembly patches, and the 16-character file
key (a SHA-256 prefix of the assembly, Ensembl's checksum and compressed size)
distinguishes upstream revisions. These are metadata checks: Ensembl's Unix
checksums are not cryptographic hashes. If the metadata is incomplete, the
assembly directory ends in `-unverified` and each release keeps its own copy.
Local FASTA files and custom mirrors are never shared.

Downloads retry transient HTTP failures and stalls of five minutes, and are
checked against the upstream size. An interrupted DNA download resumes on the
next install on POSIX systems, using only bytes the server confirms come from
the same file; on Windows it starts over.
