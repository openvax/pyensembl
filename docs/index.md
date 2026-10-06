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

## Read protein and transcript sequences <a id="read-transcript-and-protein-sequences"></a><a id="select-one-transcript"></a>

A gene can have several transcripts. Continue in the same Python session with
one selected transcript:

```python
transcript = data.transcript_by_id("ENST00000269305")
print(transcript.protein_sequence[:30])
print(transcript.sequence[:30])
print(len(transcript.protein_sequence), len(transcript.sequence))
data.close()
```

```text
MEEPQSDPSVEPPLSQETFSDLWKLLPENN
GTTTTCCCCTCCCATGTGCTCAAGACTGGC
393 2579
```

The first two lines show the first 30 amino acids of the protein and the first
30 bases of the transcript's cDNA, which begins with the 5′ UTR. The full
sequences have 393 amino acids and 2,579 bases.
The cDNA is spliced and already oriented 5′ to 3′, including
for minus-strand transcripts. Noncoding or incomplete transcripts may have no
protein sequence. Both come from the reference annotation, not from your
samples.

`data.close()` releases the open index files. You can instead write
`with EnsemblRelease(93, species="human") as data:` to close them
automatically.

## Next steps

- [Find genes and transcripts](guides/features.md) by name, ID, position or
  biotype, and choose among a gene's transcripts.
- [Get protein and transcript sequences](guides/transcripts.md): proteins by
  ID, coding sequences, UTRs, exons and genomic coordinates.
- [Read genomic DNA](guides/reference-dna.md) for introns or flanking regions.
  Whole-genome DNA is a separate, optional download.
- [Choose another species, assembly or annotation](guides/assembly-selection.md),
  or [install and manage data](guides/cache.md) from Python, check what is
  installed and free disk space.
- [Find a method](reference/genome.md#find-a-method) in the API reference.
