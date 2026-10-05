# Install data and look up a gene

This tutorial uses human GRCh38, [Ensembl](https://www.ensembl.org/) release 93.
Pinning the annotation makes the identifiers and coordinates below reproducible.
Other releases can produce different results.

## Install PyEnsembl and reference data

Use Python 3.9 or later and [pip](https://pip.pypa.io/en/stable/):

```sh
python -m pip install pyensembl
pyensembl install --release 93 --species human
```

The second command downloads matched GTF, transcript and peptide FASTA files
and builds local indexes. Installation can take several minutes and needs disk
space for compressed files and indexes. It does not download whole-genome DNA.
Check readiness with `pyensembl list`; see [cached data](guides/cache.md).

## Look up TP53

Run this in Python after installation:

```python
from pyensembl import EnsemblRelease

data = EnsemblRelease(93, species="human")
gene = data.gene_by_id("ENSG00000141510")
print(gene.id, gene.name, gene.contig, gene.start, gene.end, gene.strand)
```

```text
ENSG00000141510 TP53 17 7661779 7687550 -
```

| Field | Meaning |
| --- | --- |
| `id` | Stable gene identifier in this annotation |
| `name` | Annotated gene symbol |
| `contig` | Chromosome or other sequence name |
| `start`, `end` | One-based inclusive genomic coordinates |
| `strand` | Transcriptional direction: `"+"` or `"-"` |

Here the gene lies on chromosome 17 on the minus strand. The lower coordinate
is still `start`, regardless of strand. A genomic interval includes both ends;
its length is `end - start + 1` bases.

`data.genes_by_name("TP53")` returns all exact symbol matches. Names can be
ambiguous, including on alternative loci; choose by a stable ID or an explicit
selection rule. [Alias lookup](guides/aliases.md) accepts a separate local source.

## Select one transcript

A gene can have several transcripts. This example selects one by ID rather
than assuming the first returned transcript is a preferred isoform:

```python
transcript = data.transcript_by_id("ENST00000269305")
print(transcript.id, transcript.version, transcript.gene_id)
print(transcript.contig, transcript.start, transcript.end, transcript.strand)
```

```text
ENST00000269305 8 ENSG00000141510
17 7668402 7687538 -
```

The transcript belongs to the gene above. `version` is the annotation's
transcript version, separate from the stable ID. Its genomic span includes
introns. [Feature objects](guides/features.md) explains transcripts and exons;
the [complete reference](reference/features.md) lists their attributes.

## Read transcript and protein sequences

```python
print(len(transcript.sequence), len(transcript.protein_sequence))
data.close()
```

```text
2579 393
```

The transcript FASTA supplies 2,579 bases of spliced cDNA, already oriented
5′ to 3′; do not reverse-complement it again for a minus-strand transcript.
The peptide FASTA supplies 393 amino acids. These are annotated reference
sequences, not a sample-specific reconstruction. Noncoding or incomplete
transcripts may lack a peptide sequence.

For introns, intergenic intervals or flanking sequence, install optional
[reference DNA](guides/reference-dna.md). For a different coordinate system
or annotation, [choose an assembly and release](guides/assembly-selection.md)
before querying.
