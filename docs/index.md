# PyEnsembl <a id="install-data-and-look-up-a-gene"></a>

PyEnsembl lets you find where genes are located, which transcripts they produce,
and their sequences in Python. It uses [Ensembl](https://www.ensembl.org/)
annotations or your own files. Once the data is installed, queries run locally.

## Install <a id="install-pyensembl-and-reference-data"></a>

Use Python 3.9 or later. Run these commands in a terminal:

```sh
python -m pip install pyensembl
pyensembl install --release 93 --species human
```

The second command downloads and indexes human annotation, transcript and
protein data. It can take several minutes. These examples use GRCh38 and
Ensembl release 93 so you can reproduce the output; for your own analysis,
[choose the reference that matches your data](guides/assembly-selection.md).

## Look up a gene <a id="look-up-tp53"></a>

Run this in Python after installation. It retrieves TP53 by its gene ID:

```python
from pyensembl import EnsemblRelease

data = EnsemblRelease(93, species="human")
gene = data.gene_by_id("ENSG00000141510")
print(gene.id, gene.name, gene.contig, gene.start, gene.end, gene.strand)
```

```text
ENSG00000141510 TP53 17 7661779 7687550 -
```

The result gives the gene ID, name, chromosome, start, end and strand.
TP53 is on chromosome 17, on the minus strand. Coordinates are one-based and
include both ends; `start` is the lower coordinate on either strand.

You can also search with `data.genes_by_name("TP53")`. It returns a list because
a name can match several genes. [Gene aliases](guides/aliases.md) let you search
additional names such as p53.

## Read transcript and protein sequences <a id="select-one-transcript"></a>

A gene can have several transcripts. Continue in the same Python session with
one selected transcript:

```python
transcript = data.transcript_by_id("ENST00000269305")
print(transcript.sequence[:30])
print(transcript.protein_sequence[:30])
print(len(transcript.sequence), len(transcript.protein_sequence))
data.close()
```

```text
GTTTTCCCCTCCCATGTGCTCAAGACTGGC
MEEPQSDPSVEPPLSQETFSDLWKLLPENN
2579 393
```

The first two lines show the first 30 bases of cDNA and the first 30 amino acids
of its protein. The full sequences have 2,579 bases and 393 amino acids.
The cDNA is spliced and already oriented 5′ to 3′, including
for minus-strand transcripts. Noncoding or incomplete transcripts may have no
protein sequence. These are annotated reference sequences.

## Next steps

- [Find genes by name or genomic position and explore their transcripts](guides/features.md).
  Use the
  [API reference](reference/index.md) to find a specific method.
- [Read genomic DNA](guides/reference-dna.md) for introns or flanking regions.
  Whole-genome DNA is a separate, optional download.
- [Choose another species, assembly or annotation](guides/assembly-selection.md).
- [Manage downloaded data](guides/cache.md) to check installation, choose a cache
  location or free disk space.
