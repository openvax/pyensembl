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

To see the genomes for which PyEnsembl has already downloaded and indexed metadata you can run:

```sh
pyensembl list
```

Or equivalently do this in Python:

```python
from pyensembl.shell import collect_all_installed_ensembl_releases
collect_all_installed_ensembl_releases()
```

## Reference DNA sequences (optional)

Reference DNA enables intronic, intergenic, and flanking sequence queries. It is
**opt-in**: a normal installation does not download a whole genome. Human DNA
needs roughly 1 GB compressed and several GB of disk space after decompression.

```python
from pyensembl import EnsemblRelease

release = EnsemblRelease(81, download_genome_fasta=True)
release.download_genome_fasta()  # DNA only; reuses an existing compatible cache
release.index_genome_fasta()
with release:
    bases = release.sequence("7", 117_480_000, 117_480_100)
```

`sequence(contig, start, end)` returns **plus-strand, one-based inclusive** bases,
including for loci annotated on the minus strand. Contig names must match the
FASTA exactly (`"1"` and `"chr1"` are distinct). Coordinates must be integers with
`1 <= start <= end <= contig length`. Missing contigs and invalid ranges raise
`ValueError`. Unconfigured or uninstalled DNA raises `MissingGenomeFastaError`,
a `ValueError` subclass. Reads never download missing remote data implicitly.
The default result is uppercase; `mask="raw"` preserves soft-masked lowercase.

The lazy `.fasta` reader also supports `fasta[contig][start-1:end].seq` for
consumers such as Varcode. `.fasta` is `None` when DNA is unconfigured.
`genome_fasta_path` reports the existing uncompressed file, or `None`.
`download()` and `index()` include DNA when configured, alongside annotation,
transcript, and protein files.

```sh
# Install all data, including DNA.
pyensembl install --release 81 --with-genome-fasta

# Install DNA alone; annotation/transcript/protein data are not needed.
pyensembl install --release 81 --only-genome-fasta

# An optional subset and soft masking.
pyensembl install --release 81 --only-genome-fasta \
    --genome-fasta-type primary_assembly --masked soft
```

The default downloaded DNA is unmasked **toplevel**, covering the patch and
haplotype contigs included in Ensembl annotations. `primary_assembly` includes
chromosomes and unplaced/unlocalized sequences, but excludes patches and
haplotypes. It is unavailable for some older releases and species. Mask choices
are `none`, `soft` (lowercase repeats), and `hard` (repeats replaced with N).
The corresponding Python options are `genome_fasta_type` and `genome_fasta_mask`.
See Ensembl's [DNA file definitions](https://ftp.ensembl.org/pub/release-81/fasta/homo_sapiens/dna/README).

### Attach a local FASTA

```python
release = EnsemblRelease(81, genome_fasta_path="/data/my_reference.fa.gz")
bases = release.sequence("7", 117_480_000, 117_480_100)

# Custom annotations can also attach DNA:
from pyensembl import Genome
custom = Genome("custom", "my_annotations", genome_fasta_path_or_url="/data/reference.fa")
```

```sh
pyensembl install --release 81 --genome-fasta-path /data/my_reference.fa
# Add --only-genome-fasta to skip annotation installation.
```

Local FASTAs take precedence over canonical downloads. Plain FASTAs are read
in place; gzip/BGZF files are streamed into an uncompressed cache copy.
Indexes always live in PyEnsembl's cache, so read-only source directories work
and user-owned files/indexes remain untouched. Full indexing warns about
annotation contigs missing from a local FASTA. Matching contig names alone do
not verify assembly identity. Supply the same local path when constructing
subsequent Python objects; CLI metadata does not silently change constructors.
Custom `Genome` sources also accept a FASTA URL.

### Shared DNA cache and cleanup

Canonical Ensembl downloads use `pyensembl/dna_cache/objects/` under the cache
location described above. Compatible releases reuse one uncompressed FASTA and
FAI index. Sharing uses the **versioned assembly accession**, species, provider,
file flavor/masking, and Ensembl file metadata (Unix checksum and compressed
size). A major assembly name alone is insufficient. Local FASTAs and custom
mirrors remain separate. If upstream identity metadata is incomplete, storage
is isolated by source URL. Ensembl's Unix checksums are not cryptographic hashes.

A new release installation fetches small metadata files before deciding whether
DNA can be reused. Registered releases work offline. Each release retains its
references, including different installed masking/flavor choices. Install,
index publication, release removal, and pruning use a shared lock.

```sh
pyensembl list --check-genome-fasta  # inspect source, presence, and index
pyensembl delete-all-files --release 81  # remove this release's references
pyensembl prune --orphan-genome-fastas --dry-run
pyensembl prune --orphan-genome-fastas
```

`list` includes DNA-only installations and shows the most recently installed
DNA choice for each release. Inspection does not download or rebuild files.
Deleting a release preserves shared DNA still referenced by any release.
Pruning removes only cache-owned objects with no release references; malformed
reference metadata aborts pruning. `delete-index-files` preserves shared DNA
indexes to keep other releases usable; call `index_genome_fasta(overwrite=True)`
to rebuild one. Local source files are never pruned. Python callers can use
`prune_genome_fastas(dry_run=True)` to inspect `(path, bytes)` candidates.

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
