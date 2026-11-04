# PyEnsembl

[![Tests](https://github.com/openvax/pyensembl/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/pyensembl/actions/workflows/tests.yml)
[![Coverage](https://coveralls.io/repos/github/openvax/pyensembl/badge.svg?branch=main)](https://coveralls.io/github/openvax/pyensembl?branch=main)
[![PyPI](https://img.shields.io/pypi/v/pyensembl.svg)](https://pypi.org/project/pyensembl/)

PyEnsembl looks up genes, transcripts, exons and their protein and transcript
sequences in a pinned [Ensembl](https://www.ensembl.org/) release or in your
own GTF and FASTA files. It downloads and indexes an annotation once; after
that, queries run locally. It is the gene annotation layer of the
[OpenVax](https://github.com/openvax) tools:
[Varcode](https://github.com/openvax/varcode),
[Isovar](https://github.com/openvax/isovar),
[Topiary](https://github.com/openvax/topiary) and
[Vaxrank](https://github.com/openvax/vaxrank) use it, and it builds on
[datacache](https://github.com/openvax/datacache) for downloads and
[gtfparse](https://github.com/openvax/gtfparse) for reading GTF files.

## Install

<a id="installation"></a>Use Python 3.11 or later. Install the package and one
annotation release:

```sh
python -m pip install pyensembl
pyensembl install --release 93 --species human
```

The second command selects human GRCh38 and downloads and indexes Ensembl
release 93's annotation, transcript and protein files, about 1.3 GB. It takes
several minutes. Whole-genome DNA is a separate, optional download.

## Quick start

<a id="example-usage"></a>After installing the data above:

```python
from pyensembl import EnsemblRelease

with EnsemblRelease(93, species="human") as data:
    gene = data.gene_by_id("ENSG00000141510")
    print(gene.name, gene.contig, gene.start, gene.end, gene.strand)
    transcript = data.transcript_by_id("ENST00000269305")
    print(transcript.name, transcript.protein_sequence[:30])
```

```text
TP53 17 7661779 7687550 -
TP53-201 MEEPQSDPSVEPPLSQETFSDLWKLLPENN
```

Coordinates are one-based and inclusive; TP53 is on the minus strand. Pin the
release you analyze: coordinates, IDs and sequences can change between
releases. A gene name can match several loci, so `data.genes_by_name("TP53")`
returns a list.

## Command line

The `pyensembl` command installs and manages local data:

```sh
pyensembl install --release 93 --species human  # download and index a release
pyensembl list                                  # installed genomes and their status
pyensembl inspect --release 93                  # files, sizes and provenance, offline
pyensembl available                             # supported species, assemblies and releases
pyensembl delete-all-files --release 93         # remove one release's local data
pyensembl prune --dry-run                       # shared reference DNA that no release uses
```

The [command-line reference](https://openvax.github.io/pyensembl/reference/cli/)
lists every option. For Tab completion, run `activate-global-python-argcomplete`
once, or add `eval "$(register-python-argcomplete pyensembl)"` to your shell's
startup file.

## Documentation

<a id="usage-tips"></a>The [documentation](https://openvax.github.io/pyensembl/)
starts with installation and a first lookup, then covers one task per page. In
reading order:

1. <a id="load-genome-in-python"></a><a id="data-structures"></a><a id="gene"></a><a id="transcript"></a>[Find genes and transcripts](https://openvax.github.io/pyensembl/guides/features/)
   by name, ID, position or biotype.
2. <a id="protein-information"></a>[Get protein and transcript sequences](https://openvax.github.io/pyensembl/guides/transcripts/),
   coding sequences, UTRs and exons.
3. <a id="reference-dna-optional"></a><a id="reading-sequences"></a><a id="choosing-dna"></a><a id="local-fasta-files"></a><a id="upgrading-from-2110"></a>[Read genomic DNA](https://openvax.github.io/pyensembl/guides/reference-dna/)
   for introns and flanking sequence.
4. <a id="annotation-coverage"></a><a id="list-supported-species"></a><a id="new-ensembl-platform"></a>[Choose a reference](https://openvax.github.io/pyensembl/guides/assembly-selection/):
   species, assembly, and numbered or dated Ensembl releases.
5. <a id="cache-location"></a><a id="list-installed-genomes"></a><a id="inspect-data-without-installing"></a><a id="managing-disk-space"></a><a id="how-the-shared-dna-cache-works"></a>[Install and manage data](https://openvax.github.io/pyensembl/guides/cache/)
   from Python, check what is installed and free disk space.
6. <a id="non-ensembl-data"></a>[Use your own annotations](https://openvax.github.io/pyensembl/guides/custom-genomes/)
   from GTF and FASTA files.

<a id="look-up-gene-name-aliases"></a><a id="api"></a><a id="genes"></a><a id="transcripts"></a><a id="exons"></a><a id="reference-dna"></a>Also
see [gene name aliases](https://openvax.github.io/pyensembl/guides/aliases/),
[genetic codes](https://openvax.github.io/pyensembl/guides/genetic-codes/) and
the [API reference](https://openvax.github.io/pyensembl/reference/).

## Development

<a id="development-setup"></a>From a clone of this repository:

```sh
python -m pip install -e '.[dev]'
./lint.sh   # ruff and mypy
./test.sh   # pytest with coverage
./docs.sh   # strict docs build and link checks; needs '.[docs]'
```

Most tests need the Ensembl releases that
[the Tests workflow](https://github.com/openvax/pyensembl/blob/main/.github/workflows/tests.yml)
installs. Every PR bumps `pyensembl/version.py`, including documentation
changes. See [contributing](https://openvax.github.io/pyensembl/dev/contributing/)
and [releasing](https://openvax.github.io/pyensembl/dev/releasing/).

## License

[Apache License 2.0](https://github.com/openvax/pyensembl/blob/main/LICENSE)
