# Reference DNA (optional)

Reference DNA lets you read any genomic interval, including introns,
intergenic regions, and flanking sequence. It is **opt-in**: a normal
installation does not download a whole genome. Human DNA is about 1 GB to
download and takes several GB of disk once decompressed.

## Quick start

```sh
pyensembl install --release 93 --with-genome-fasta  # annotation and DNA
pyensembl install --release 93 --only-genome-fasta  # just the DNA
```

```python
from pyensembl import EnsemblRelease

release = EnsemblRelease(93, genome_fasta=True)
release.download_genome_fasta()  # does nothing if the DNA is already installed
with release:
    bases = release.sequence("7", 117_480_000, 117_480_100)
    tp53 = release.gene_by_id("ENSG00000141510")  # TP53, on the minus strand
    tp53_dna = release.sequence(tp53.contig, tp53.start, tp53.end, strand=tp53.strand)
```

`genome_fasta=True` only chooses the DNA. Nothing is downloaded until you call
`download_genome_fasta()` or `download()`, or run `pyensembl install`.
Python objects use reference DNA only when constructed with `genome_fasta`:
a plain `EnsemblRelease(93)` does not pick up DNA installed by the CLI, and its
error message names the call that does.

## Reading sequences

`sequence(contig, start, end, mask="upper", *, strand="+")`:

- Coordinates are **one-based and inclusive**, like the rest of PyEnsembl,
  with `1 <= start <= end <= contig length`.
- Bases are read from the plus strand; `strand="-"` returns the reverse
  complement, so a gene or transcript reads 5' to 3'.
- Contigs can be named as in the FASTA or as PyEnsembl reports them
  (`gene.contig`). If the names differ only by a `chr` prefix, the error
  suggests the right one.
- Results are uppercase; `mask="raw"` keeps soft-masked repeats in lowercase.
- Absent contigs and invalid ranges raise `ValueError`. Missing DNA raises
  `MissingGenomeFastaError`, a `ValueError` whose message explains how to
  install or enable it. Reads never download anything.

Related attributes and methods:

- `fasta` is a [pyfaidx](https://github.com/mdshw5/pyfaidx) reader with
  zero-based, half-open slices (`fasta[contig][start - 1:end].seq`), as used by
  Varcode. It needs the FASTA's own contig names, and is `None` when DNA is not
  configured or not installed.
- `genome_fasta_path` is the uncompressed FASTA on disk, or `None`.
- `download()` and `index()` include DNA when it is configured.
  `index_genome_fasta()` builds the DNA index before the first query needs it.
- `close()` closes the reader. Readers already handed out stay usable after
  `clear_cache()`.
- Attached DNA does not affect equality: genes and transcripts from the same
  release compare equal with or without it.

## Choosing DNA

| | Python (`EnsemblRelease`) | CLI (`pyensembl install`) |
|---|---|---|
| Ensembl's DNA | `genome_fasta=True` | `--with-genome-fasta` or `--only-genome-fasta` |
| A local FASTA | `genome_fasta="/data/ref.fa.gz"` | `--genome-fasta-path /data/ref.fa.gz` |
| Coverage | `genome_fasta_type="primary_assembly"` | `--genome-fasta-type primary_assembly` |
| Masking | `genome_fasta_mask="soft"` | `--masked soft` |

The default is unmasked **toplevel** DNA, which covers the patch and haplotype
contigs in Ensembl annotations. `primary_assembly` has the chromosomes and
unplaced/unlocalized sequences but no patches or haplotypes, and some older
releases and species don't provide it. Masking is `none`, `soft` (repeats in
lowercase), or `hard` (repeats replaced with `N`). See Ensembl's
[DNA file definitions](https://ftp.ensembl.org/pub/release-81/fasta/homo_sapiens/dna/README).

## Local FASTA files

```python
release = EnsemblRelease(93, genome_fasta="/data/my_reference.fa.gz")

# Custom annotations can attach DNA too, from a path or URL:
from pyensembl import Genome
custom = Genome("custom", "my_annotations", genome_fasta_path_or_url="/data/reference.fa")
```

Plain FASTA files are read in place; gzip and BGZF files are decompressed into
PyEnsembl's cache on first use. Indexes always live in the cache, so read-only
source directories work and your own files and indexes are never modified.
`index()` warns about annotation contigs that are missing from a local FASTA,
but matching contig names don't prove that the assembly matches.

## Managing disk space

```sh
pyensembl list --check-genome-fasta      # DNA for each release, verifying indexes
pyensembl delete-all-files --release 93  # release 93's files and DNA references
pyensembl prune --dry-run                # shared DNA that no installed release uses
pyensembl prune
```

Compatible releases share one copy of Ensembl DNA, so deleting a release keeps
DNA that other releases still use, and `prune` removes DNA that no release
references. It skips DNA that is being downloaded or indexed, never touches
local FASTA files, and deletes nothing if any release's DNA metadata is
malformed (`list` shows which one). `list` includes DNA-only installs and
shows each release's most recently installed DNA. `delete-index-files` keeps
shared DNA indexes because other releases may use them; rebuild one with
`index_genome_fasta(overwrite=True)`. In Python,
`prune_genome_fastas(dry_run=True)` returns `(path, bytes)` candidates.

## How the shared DNA cache works

Ensembl DNA is stored once per upstream file under `pyensembl/dna_cache/`, in
`<species>/<provider>/<reference>-<assembly accession>/<coverage>/<masking>/fasta/<file key>/`:

```text
pyensembl/dna_cache/
  homo_sapiens/ftp.ensembl.org/GRCh38-GCA_000001405.18/
    toplevel/unmasked/fasta/<file key>/
      sequence.fa        uncompressed, even when downloaded as .fa.gz
      sequence.fa.fai
      object.json        full identity of the upstream file
      index.json
```

Installing a release first reads Ensembl's small README and CHECKSUMS files to
see whether another release already downloaded the same file. The versioned
assembly accession distinguishes assembly patches. The 16-character file key
is a prefix of the SHA-256 of the file's identity (assembly, Ensembl's Unix
checksum, and compressed size), which separates upstream revisions of the same
file, and conflicting identities are never reused. These are metadata checks:
Ensembl's Unix checksums are not cryptographic hashes. If the metadata is
incomplete, the assembly directory ends in `-unverified` and each release keeps
its own copy. Local FASTA files and custom mirrors are never shared. Downloads
retry transient HTTP failures and are checked against the upstream size, and
installed releases work offline. An interrupted Ensembl DNA download resumes
where it stopped the next time you install; datacache appends only bytes that
Ensembl's server confirms come from the same file (its ETag). Resuming needs a
POSIX system; on Windows an interrupted download starts over. A download that
receives no data for five minutes is retried.

Reads take no locks and write nothing, so a fully installed and indexed cache
can be read-only for other users. A download or index build locks only the
file it writes; registering, deleting, and pruning releases briefly lock the
whole cache. Files follow your umask, as do lock files on Python 3.10+, so use
`umask 002` or default ACLs for a group-shared cache. `dna_cache` itself may be
a symlink, e.g. to a larger disk.

## Upgrading from 2.11.0

| 2.11.0 | 2.12.0 and later |
|---|---|
| `EnsemblRelease(81, download_genome_fasta=True)` | `EnsemblRelease(81, genome_fasta=True)` |
| `EnsemblRelease(81, genome_fasta_path="/data/ref.fa")` | `EnsemblRelease(81, genome_fasta="/data/ref.fa")` |

The old keywords still work but emit a `DeprecationWarning`. Objects pickled or
serialized by 2.11.0 still load. CLI flags are unchanged.
