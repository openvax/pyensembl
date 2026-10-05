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
DNA is optional. See [installation and first lookup](https://github.com/openvax/pyensembl/blob/main/docs/getting-started.md).

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

Coordinates are one-based and inclusive; TP53 is on the minus strand. The
[tutorial](https://github.com/openvax/pyensembl/blob/main/docs/getting-started.md) continues with a selected transcript and
explains identifiers, versions and sequence lengths. Pin an annotation version:
coordinates and sequences can differ between releases, and names may match
several loci.

## Usage tips

The [documentation](https://github.com/openvax/pyensembl/blob/main/docs/index.md) separates task guides and complete reference.
Existing README anchors below lead to the corresponding guidance.

| Task | Guide |
| --- | --- |
| <a id="annotation-coverage"></a><a id="list-supported-species"></a>Choose an assembly, species and coverage | [Assembly selection](https://github.com/openvax/pyensembl/blob/main/docs/guides/assembly-selection.md) |
| <a id="new-ensembl-platform"></a>Use dated Ensembl downloads | [New platform](https://github.com/openvax/pyensembl/blob/main/docs/ensembl-platform.md) |
| <a id="non-ensembl-data"></a>Load custom matched files | [Custom genomes](https://github.com/openvax/pyensembl/blob/main/docs/guides/custom-genomes.md) |
| <a id="cache-location"></a><a id="list-installed-genomes"></a><a id="inspect-data-without-installing"></a>Inspect data and cache readiness | [Cached data](https://github.com/openvax/pyensembl/blob/main/docs/guides/cache.md) |
| <a id="look-up-gene-name-aliases"></a>Load human or custom aliases | [Alias lookup](https://github.com/openvax/pyensembl/blob/main/docs/guides/aliases.md) |
| <a id="load-genome-in-python"></a><a id="data-structures"></a><a id="gene"></a><a id="transcript"></a><a id="protein-information"></a>Work with feature objects | [Genes and transcripts](https://github.com/openvax/pyensembl/blob/main/docs/guides/features.md) |
| <a id="reference-dna-optional"></a><a id="quick-start"></a><a id="reading-sequences"></a><a id="choosing-dna"></a><a id="local-fasta-files"></a><a id="managing-disk-space"></a><a id="how-the-shared-dna-cache-works"></a><a id="upgrading-from-2110"></a>Read genomic intervals and manage DNA | [Reference DNA](https://github.com/openvax/pyensembl/blob/main/docs/guides/reference-dna.md) |
| <a id="api"></a><a id="genes"></a><a id="transcripts"></a><a id="exons"></a><a id="reference-dna"></a>Find complete interfaces and lookup methods | [API reference](https://github.com/openvax/pyensembl/blob/main/docs/reference/index.md) |

<a id="development-setup"></a>
For contributions, see [development setup](https://github.com/openvax/pyensembl/blob/main/docs/dev/contributing.md),
[documentation style](https://github.com/openvax/pyensembl/blob/main/docs/dev/documentation-style.md) and
[releasing](https://github.com/openvax/pyensembl/blob/main/docs/dev/releasing.md). Every PR includes a version bump, lint and
tests; maintainers merge through a PR and publish with `./deploy.sh` from clean main.
