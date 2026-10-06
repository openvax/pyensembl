[![Tests](https://github.com/openvax/pyensembl/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/pyensembl/actions/workflows/tests.yml)
[![Coverage](https://coveralls.io/repos/github/openvax/pyensembl/badge.svg?branch=main)](https://coveralls.io/github/openvax/pyensembl?branch=main)
[![PyPI](https://img.shields.io/pypi/v/pyensembl.svg)](https://pypi.org/project/pyensembl/)

# PyEnsembl

PyEnsembl queries [Ensembl](https://www.ensembl.org/) genes, transcripts, exons
and sequences from local indexes. It supports numbered releases, dated
new-platform annotations and custom matched GTF/FASTA files.

## Installation

Use Python 3.9 or later. Install the package and one explicit annotation:

```sh
python -m pip install pyensembl
pyensembl install --release 93 --species human
```

This selects human GRCh38 and downloads and indexes its annotation, transcript
and peptide files. Setup takes several minutes and disk space; whole-genome
DNA is optional. See [installation and first lookup](https://openvax.github.io/pyensembl/).

## Example Usage

After installing the data above:

```python
from pyensembl import EnsemblRelease

with EnsemblRelease(93, species="human") as data:
    gene = data.gene_by_id("ENSG00000141510")
    print(gene.name, gene.contig, gene.start, gene.end, gene.strand)
```

```text
TP53 17 7661779 7687550 -
```

Coordinates are one-based and inclusive; TP53 is on the minus strand. Pin the
release you analyze: coordinates, IDs and sequences can change between
releases, and a gene name can match several loci. The
[documentation](https://openvax.github.io/pyensembl/) continues with transcript
and protein sequences.

## Documentation <a id="usage-tips"></a>

The [documentation home page](https://openvax.github.io/pyensembl/) includes installation and
getting-started examples. Use these guides for other tasks:

| Task | Guide |
| --- | --- |
| <a id="load-genome-in-python"></a><a id="data-structures"></a><a id="gene"></a><a id="transcript"></a>Find genes by name, position or biotype | [Genes and transcripts](https://openvax.github.io/pyensembl/guides/features/) |
| <a id="protein-information"></a>Exons, coding sequences, UTRs and codon positions | [Exons and coding sequences](https://openvax.github.io/pyensembl/guides/transcripts/) |
| <a id="look-up-gene-name-aliases"></a>Search additional gene names | [Alias lookup](https://openvax.github.io/pyensembl/guides/aliases/) |
| <a id="reference-dna-optional"></a><a id="quick-start"></a><a id="reading-sequences"></a><a id="choosing-dna"></a><a id="local-fasta-files"></a><a id="upgrading-from-2110"></a>Read genomic DNA | [Reference DNA](https://openvax.github.io/pyensembl/guides/reference-dna/) |
| <a id="annotation-coverage"></a><a id="list-supported-species"></a>Choose an assembly, species and coverage | [Assembly selection](https://openvax.github.io/pyensembl/guides/assembly-selection/) |
| <a id="cache-location"></a><a id="list-installed-genomes"></a><a id="inspect-data-without-installing"></a><a id="managing-disk-space"></a><a id="how-the-shared-dna-cache-works"></a>Install from Python, check or delete data | [Install and manage data](https://openvax.github.io/pyensembl/guides/cache/) |
| <a id="non-ensembl-data"></a>Use your own annotations | [Custom genomes](https://openvax.github.io/pyensembl/guides/custom-genomes/) |
| <a id="new-ensembl-platform"></a>Use dated Ensembl downloads | [Dated annotations](https://openvax.github.io/pyensembl/ensembl-platform/) |
| <a id="api"></a><a id="genes"></a><a id="transcripts"></a><a id="exons"></a><a id="reference-dna"></a>Find a class or query method | [API reference](https://openvax.github.io/pyensembl/reference/) |

<a id="development-setup"></a>
For contributions, see [development setup](https://openvax.github.io/pyensembl/dev/contributing/),
[documentation style](https://openvax.github.io/pyensembl/dev/documentation-style/) and
[releasing](https://openvax.github.io/pyensembl/dev/releasing/). Every PR includes a version bump, lint and
tests; maintainers merge through a PR and publish with `./deploy.sh` from clean main.
