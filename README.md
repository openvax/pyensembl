[![Tests](https://github.com/openvax/pyensembl/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/pyensembl/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/github/openvax/pyensembl/badge.svg?branch=main)](https://coveralls.io/github/openvax/pyensembl?branch=main)
<a href="https://pypi.python.org/pypi/pyensembl/">
<img src="https://img.shields.io/pypi/v/pyensembl.svg?maxAge=1000" alt="PyPI" />
</a>

# PyEnsembl

PyEnsembl is a Python interface to [Ensembl](http://www.ensembl.org) reference genome metadata such as exons and transcripts. PyEnsembl downloads [GTF](https://en.wikipedia.org/wiki/Gene_transfer_format) and [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files from the [Ensembl FTP server](https://ftp.ensembl.org/pub/) and loads them into a local database. PyEnsembl can also work with custom reference data specified using user-supplied GTF and FASTA files.

# Example Usage

```python
from pyensembl import EnsemblRelease

# release 77 uses human reference genome GRCh38
data = EnsemblRelease(77)

# will return ['HLA-A']
gene_names = data.gene_names_at_locus(contig=6, position=29945884)

# get all exons associated with HLA-A
exon_ids  = data.exon_ids_of_gene_name('HLA-A')
```

# Installation

PyEnsembl requires Python 3.9 or later. You can install PyEnsembl using [pip](https://pip.pypa.io/en/latest/quickstart.html):

```sh
pip install pyensembl
```

This should also install any required packages such as [datacache](https://github.com/openvax/datacache).

Before using PyEnsembl, run the following command to download and install
Ensembl data:

```
pyensembl install --release <list of Ensembl release numbers> --species <species-name>
```

For example, `pyensembl install --release 75 76 --species human` will download and install all
human reference data from Ensembl releases 75 and 76.

To install the newest supported Ensembl release for a reference assembly:

```sh
pyensembl install --reference-name GRCh37
```

Reference names are case-insensitive. This selects human release 75 for GRCh37;
the species is inferred from the reference. You can also specify `--release`
to select older releases for that assembly. Conflicting `--species` or
`--release` selections are rejected before downloading data. Deletion commands
still require an explicit `--release`.

Alternatively, you can create the `EnsemblRelease` object from inside a Python
process and call `ensembl_object.download()` followed by `ensembl_object.index()`.

Whole-genome reference DNA, for intronic and intergenic sequence, is optional
and not installed by default; see [Reference DNA](#reference-dna-optional).

## Annotation coverage

PyEnsembl uses Ensembl's complete `chr_patch_hapl_scaff` GTF for human GRCh38
from release 82, mouse GRCm38 releases 82–102, and zebrafish GRCz11 from release
92. These files include additional genes on assembly patches and haplotypes.
Other assemblies and earlier releases use the standard GTF filename.

Patch and haplotype contig names are preserved, for example
`CHR_HG2263_PATCH`. Gene-name searches can return additional genes on these
contigs; use stable gene IDs or a contig filter when selecting a particular locus.

After upgrading from versions before 2.10.17, rerun installation for each affected
release you use, for example `pyensembl install --release 97 --species human`.
The complete GTF creates a separate index, so an older index cannot hide the
additional genes. Existing source files and indexes are retained, and unchanged
FASTA files are reused. Custom mirrors must provide the complete GTF filename;
to use a deliberately restricted annotation, supply its GTF as custom data.

## Development Setup

For development, install PyEnsembl in editable mode with development dependencies:

```sh
git clone https://github.com/openvax/pyensembl.git
cd pyensembl
pip install -e .[dev]
```

This installs the package in development mode along with tools for testing, linting, and building:
- `pytest` for running tests
- `ruff` for code linting
- `pytest-cov` for coverage reporting
- `build` for package building

Run lint and tests with:
```sh
./lint.sh
./test.sh
```

Most tests need Ensembl data installed first; `.github/workflows/tests.yml`
lists the releases CI installs.

Species assembly ranges are checked against Ensembl's archive. After raising
`MAX_ENSEMBL_RELEASE`, recheck every assembly boundary on the live FTP servers:

```sh
PYENSEMBL_NETWORK_TESTS=1 ./test.sh tests/test_species_assemblies.py
```

## Cache Location

By default, PyEnsembl uses the platform-specific `Cache` folder
and caches the files into the `pyensembl` sub-directory.
You can override this default by setting the environment key `PYENSEMBL_CACHE_DIR`
as your preferred location for caching:

```sh
export PYENSEMBL_CACHE_DIR=/custom/cache/dir
```

or

```python
import os

os.environ['PYENSEMBL_CACHE_DIR'] = '/custom/cache/dir'
# ... PyEnsembl API usage
```

# Usage tips

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
`incomplete` means some downloads are missing. Run `pyensembl install` for
that release to finish (add `--species` for non-human genomes; custom
genomes need their original install options). Custom genomes are listed on
Linux and macOS, or wherever `PYENSEMBL_CACHE_DIR` is set.

`install` prints progress on stderr, one line per step; add `--verbose` (`-v`)
to see every download and database step.

To get the installed Ensembl releases in Python:

```python
from pyensembl.shell import collect_all_installed_ensembl_releases
collect_all_installed_ensembl_releases()
```

## List supported species

To see every species PyEnsembl knows about, with its assemblies and supported
Ensembl release ranges:

```sh
pyensembl available
```

## Load genome in Python

Here's an example Python snippet that loads fly genome data from Ensembl release v100:

```python
from pyensembl import EnsemblRelease
data = EnsemblRelease(release=100, species='drosophila_melanogaster')
```

## Data structures

### Gene

```python
gene = data.gene_by_id(gene_id='FBgn0011747')
```

### Transcript

```python
transcript = gene.transcripts[0]
```

### Protein information

```python
transcript.protein_id
transcript.protein_sequence
```

# Non-Ensembl Data

PyEnsembl also allows arbitrary genomes via the specification
of local file paths or remote URLs to both Ensembl and non-Ensembl GTF
and FASTA files. (Warning: GTF formats can vary, and handling of
non-Ensembl data is still very much in development.)

For example:

```python
from pyensembl import Genome
data = Genome(
    reference_name='GRCh38',
    annotation_name='my_genome_features',
    # annotation_version=None,
    gtf_path_or_url='/My/local/gtf/path_to_my_genome_features.gtf', # Path or URL of GTF file
    # transcript_fasta_paths_or_urls=None, # List of paths or URLs of FASTA files containing transcript sequences
    # protein_fasta_paths_or_urls=None, # List of paths or URLs of FASTA files containing protein sequences
    # cache_directory_path=None, # Where to place downloaded and cached files for this genome
)
# parse GTF and construct database of genomic features
data.index()
gene_names = data.gene_names_at_locus(contig=6, position=29945884)
```

# Reference DNA (optional)

Reference DNA lets you read any genomic interval, including introns,
intergenic regions, and flanking sequence. It is **opt-in**: a normal
installation does not download a whole genome. Human DNA is about 1 GB to
download and takes several GB of disk once decompressed.

## Quick start

```sh
pyensembl install --release 81 --with-genome-fasta  # annotation and DNA
pyensembl install --release 81 --only-genome-fasta  # just the DNA
```

```python
from pyensembl import EnsemblRelease

release = EnsemblRelease(81, genome_fasta=True)
release.download_genome_fasta()  # does nothing if the DNA is already installed
with release:
    bases = release.sequence("7", 117_480_000, 117_480_100)
    tp53 = release.genes_by_name("TP53")[0]  # annotated on the minus strand
    tp53_dna = release.sequence(tp53.contig, tp53.start, tp53.end, strand=tp53.strand)
```

`genome_fasta=True` only chooses the DNA. Nothing is downloaded until you call
`download_genome_fasta()` or `download()`, or run `pyensembl install`.
Python objects use reference DNA only when constructed with `genome_fasta`:
a plain `EnsemblRelease(81)` does not pick up DNA installed by the CLI, and its
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
release = EnsemblRelease(81, genome_fasta="/data/my_reference.fa.gz")

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
pyensembl delete-all-files --release 81  # release 81's files and DNA references
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
installed releases work offline.

Reads take no locks and write nothing, so a fully installed and indexed cache
can be read-only for other users. A download or index build locks only the
file it writes; registering, deleting, and pruning releases briefly lock the
whole cache. Files follow your umask, as do lock files on Python 3.10+, so use
`umask 002` or default ACLs for a group-shared cache. `dna_cache` itself may be
a symlink, e.g. to a larger disk.

Sharing needs release caches next to `dna_cache`. That is the case on Linux
and macOS, and whenever `PYENSEMBL_CACHE_DIR` is set. Windows' default cache
layout is different, so there each release keeps its own DNA unless
`PYENSEMBL_CACHE_DIR` is set.

## Upgrading from 2.11.0

| 2.11.0 | 2.12.0 and later |
|---|---|
| `EnsemblRelease(81, download_genome_fasta=True)` | `EnsemblRelease(81, genome_fasta=True)` |
| `EnsemblRelease(81, genome_fasta_path="/data/ref.fa")` | `EnsemblRelease(81, genome_fasta="/data/ref.fa")` |

The old keywords still work but emit a `DeprecationWarning`. Objects pickled or
serialized by 2.11.0 still load. CLI flags are unchanged.

# API

The `EnsemblRelease` object has methods to let you access all possible
combinations of the annotation features _gene_name_, _gene_id_,
_transcript_name_, _transcript_id_, _exon_id_ as well as the location of
these genomic elements (contig, start position, end position, strand).

## Genes

<dl>
<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20genes(">genes(contig=None, strand=None, biotype=None)</a></dt>
<dd>Returns a list of <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/gene.py#:~:text=class%20Gene(">Gene</a> objects, optionally restricted to a particular contig,
strand, or <code>gene_biotype</code>.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20genes_at_locus(">genes_at_locus(contig, position, end=None, strand=None)</a></dt>
<dd>Returns a list of <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/gene.py#:~:text=class%20Gene(">Gene</a> objects overlapping a particular position on a contig,
optionally extend into a range with the end parameter and restrict to
forward or backward strand by passing strand='+' or strand='-'.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_by_id(">gene_by_id(gene_id)</a></dt>
<dd>Return a <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/gene.py#:~:text=class%20Gene(">Gene</a> object for given Ensembl gene ID (e.g. "ENSG00000068793").</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_names(">gene_names(contig=None, strand=None)</a></dt>
<dd>Returns all gene names in the annotation database, optionally restricted
to a particular contig or strand.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20genes_by_name(">genes_by_name(gene_name)</a></dt>
<dd>Get all the unique genes with the given name (there might be multiple
due to copies in the genome), return a list containing a <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/gene.py#:~:text=class%20Gene(">Gene</a> object for each
distinct ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_by_protein_id(">gene_by_protein_id(protein_id)</a></dt>
<dd>Find <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/gene.py#:~:text=class%20Gene(">Gene</a> associated with the given Ensembl protein ID (e.g. "ENSP00000350283")</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_names_at_locus(">gene_names_at_locus(contig, position, end=None, strand=None)</a></dt>
<dd>Names of genes overlapping with the given locus, optionally restricted by strand.
(returns a list to account for overlapping genes)</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_name_of_gene_id(">gene_name_of_gene_id(gene_id)</a></dt>
<dd>Returns name of gene with given gene ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_name_of_transcript_id(">gene_name_of_transcript_id(transcript_id)</a></dt>
<dd>Returns name of gene associated with given transcript ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_name_of_transcript_name(">gene_name_of_transcript_name(transcript_name)</a></dt>
<dd>Returns name of gene associated with given transcript name.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_name_of_exon_id(">gene_name_of_exon_id(exon_id)</a></dt>
<dd>Returns name of gene associated with given exon ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_ids(">gene_ids(contig=None, strand=None, biotype=None)</a></dt>
<dd>Return all gene IDs in the annotation database, optionally restricted by
chromosome name, strand, or <code>gene_biotype</code>.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20gene_ids_of_gene_name(">gene_ids_of_gene_name(gene_name)</a></dt>
<dd>Returns all Ensembl gene IDs with the given name.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20nearest_gene(">nearest_gene(contig, position, end=None, strand=None)</a></dt>
<dd>Returns <code>(distance, <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/gene.py#:~:text=class%20Gene(">Gene</a>)</code> for the gene whose locus is nearest to the
position (or position..end interval) on the given contig — even when no
gene overlaps. Returns <code>(inf, None)</code> when no candidates exist.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20merged_gene_intervals(">merged_gene_intervals(contig, strand=None)</a></dt>
<dd>Returns the union of all gene loci on the contig as a sorted list of
non-overlapping <code>(start, end)</code> tuples. Adjacent intervals
(<code>end+1 == next start</code>) are merged into one.</dd>

</dl>

## Transcripts

<dl>
<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcripts(">transcripts(contig=None, strand=None, biotype=None)</a></dt>
<dd>Returns a list of <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/transcript.py#:~:text=class%20Transcript(">Transcript</a> objects for all transcript entries in the
Ensembl database, optionally restricted to a particular contig, strand, or
<code>transcript_biotype</code>.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcript_by_id(">transcript_by_id(transcript_id)</a></dt>
<dd>Construct a <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/transcript.py#:~:text=class%20Transcript(">Transcript</a> object for given Ensembl transcript ID (e.g. "ENST00000369985")</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcripts_by_name(">transcripts_by_name(transcript_name)</a></dt>
<dd>Returns a list of <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/transcript.py#:~:text=class%20Transcript(">Transcript</a> objects for every transcript matching the given name.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcript_names(">transcript_names(contig=None, strand=None)</a></dt>
<dd>Returns all transcript names in the annotation database.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcript_ids(">transcript_ids(contig=None, strand=None, biotype=None)</a></dt>
<dd>Returns all transcript IDs in the annotation database.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcript_ids_of_gene_id(">transcript_ids_of_gene_id(gene_id)</a></dt>
<dd>Return IDs of all transcripts associated with given gene ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcript_ids_of_gene_name(">transcript_ids_of_gene_name(gene_name)</a></dt>
<dd>Return IDs of all transcripts associated with given gene name.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcript_ids_of_transcript_name(">transcript_ids_of_transcript_name(transcript_name)</a></dt>
<dd>Find all Ensembl transcript IDs with the given name.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20transcript_ids_of_exon_id(">transcript_ids_of_exon_id(exon_id)</a></dt>
<dd>Return IDs of all transcripts associated with given exon ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20nearest_transcript(">nearest_transcript(contig, position, end=None, strand=None)</a></dt>
<dd>Returns <code>(distance, <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/transcript.py#:~:text=class%20Transcript(">Transcript</a>)</code> to the closest transcript on the contig.
Returns <code>(inf, None)</code> when no candidates exist.</dd>
</dl>

## Exons

<dl>
<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20exon_ids(">exon_ids(contig=None, strand=None)</a></dt>
<dd>Returns a list of exon IDs in the annotation database, optionally restricted
by the given chromosome and strand.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20exon_by_id(">exon_by_id(exon_id)</a></dt>
<dd>Construct an <a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/exon.py#:~:text=class%20Exon(">Exon</a> object for given Ensembl exon ID (e.g. "ENSE00001209410")</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20exon_ids_of_gene_id(">exon_ids_of_gene_id(gene_id)</a></dt>
<dd>Returns a list of exon IDs associated with a given gene ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20exon_ids_of_gene_name(">exon_ids_of_gene_name(gene_name)</a></dt>
<dd>Returns a list of exon IDs associated with a given gene name.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20exon_ids_of_transcript_id(">exon_ids_of_transcript_id(transcript_id)</a></dt>
<dd>Returns a list of exon IDs associated with a given transcript ID.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20exon_ids_of_transcript_name(">exon_ids_of_transcript_name(transcript_name)</a></dt>
<dd>Returns a list of exon IDs associated with a given transcript name.</dd>
</dl>

## Reference DNA

These need reference DNA; see [Reference DNA](#reference-dna-optional).

<dl>
<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20sequence(">sequence(contig, start, end, mask="upper", strand="+")</a></dt>
<dd>Returns the bases from <code>start</code> to <code>end</code> (one-based, inclusive) on the plus strand, or their reverse complement with <code>strand="-"</code>. <code>mask="raw"</code> keeps soft-masked lowercase.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20download_genome_fasta(">download_genome_fasta(overwrite=False)</a></dt>
<dd>Downloads the configured reference DNA without annotation data; does nothing if it is already installed.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20index_genome_fasta(">index_genome_fasta(overwrite=False)</a></dt>
<dd>Builds the DNA index now rather than on the first query.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20fasta(">fasta</a></dt>
<dd>A <a href="https://github.com/mdshw5/pyfaidx">pyfaidx</a> reader with zero-based, half-open slices, or <code>None</code> when DNA is not configured or not installed.</dd>

<dt><a href="https://github.com/openvax/pyensembl/blob/main/pyensembl/genome.py#:~:text=def%20genome_fasta_path(">genome_fasta_path</a></dt>
<dd>Path of the installed, uncompressed FASTA, or <code>None</code>.</dd>
</dl>
