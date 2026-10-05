# Manage and inspect cached data

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

or

```python
import os

os.environ['PYENSEMBL_CACHE_DIR'] = '/custom/cache/dir'
# ... PyEnsembl API usage
```

To share a cache with a group, set a group-friendly umask such as `umask 002`
before installing: new files follow it. Files downloaded before PyEnsembl
2.13.1 (datacache 1.10.0) were readable only by the user who installed them;
share an existing cache with `chmod -R g+rX "$PYENSEMBL_CACHE_DIR/pyensembl"`
(or the platform cache directory).

## List installed genomes

To see which genomes are in the local cache, whether each is ready to use, and
any reference DNA:

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
empty or not a regular file; `inaccessible` means it cannot be read.
Run `pyensembl install` for
that release to finish (add `--species` for non-human genomes; custom
genomes need their original install options).

`install` prints progress on stderr, one line per step, with progress bars for
downloads, reading GTF and sequence files, and database builds when run in a
terminal.
Add `--verbose` (`-v`) to see every download and database step. In Python, pass
`show_progress=True` to `download()`, `index()`, `download_genome_fasta()`, or
`index_genome_fasta()`.

In Python, `installed()` is `True` when a genome is ready: everything it is
configured with, reference DNA included, is downloaded and indexed. It only
reads the cache, so checking never downloads or creates files.
`genome_for_reference_name` picks the newest installed release of an assembly,
else the newest downloaded one, else the newest Ensembl release:

```python
from pyensembl import EnsemblRelease, genome_for_reference_name
EnsemblRelease(93).installed()
genome_for_reference_name("GRCh38")

# Ensembl releases with any files in the cache, ready or not:
from pyensembl.shell import collect_all_installed_ensembl_releases
collect_all_installed_ensembl_releases()
```

## Inspect data without installing

Inspect a selected genome's configured source files and indexes, including
paths, availability, sizes and any recorded download provenance:

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
pickle deserialization. Source files must be readable, regular and nonempty;
SQLite indexes must have the current completed schema. This is a readiness
check, not a full validation of biological contents or pickle integrity.
An available file is not necessarily checksum-verified: provenance receipts
are advisory, not trusted expected checksums. Old files without receipts
remain usable; missing, malformed or stale receipts simply omit provenance.

New annotation downloads and copied local imports use datacache's atomic
publication and record provenance (URL credentials and query strings are
redacted). Failed overwrites preserve the previous destination. Existing
cache paths are unchanged. Invalid cached files are reported; use
`install --overwrite` to explicitly replace them. In Python,
`Genome(..., copy_local_files_to_cache=True)` makes an independent cached
import that remains usable after the original is removed, and
`decompress_on_download=True` applies to both downloaded and copied sources.
Local sources attached without copying are never modified.

For DNA disk usage, shared-cache identity, locking and pruning, see [Reference DNA](reference-dna.md#managing-disk-space).
