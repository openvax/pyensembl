# Work with genes, transcripts and proteins

The [first lookup tutorial](../getting-started.md) uses pinned human GRCh38 data.
A gene may have several transcripts. Select a transcript by a stated rule or
stable ID; list order does not identify a preferred isoform. Protein sequence
is available only when a matching peptide record exists.

## Load genome in Python

This separate example selects fly release 100. Install it before running the
Python fragment:

```sh
pyensembl install --release 100 --species drosophila_melanogaster
```

```python
from pyensembl import EnsemblRelease

data = EnsemblRelease(100, species="drosophila_melanogaster")
gene = data.gene_by_id("FBgn0011747")
transcripts = gene.transcripts
```

## Data structures

| API object | Represents | Detail |
| --- | --- | --- |
| `Gene` | An annotated locus with one gene ID | [Gene reference](../reference/features.md#pyensembl.Gene) |
| `Transcript` | One transcript and its ordered exons | [Transcript reference](../reference/features.md#pyensembl.Transcript) |
| `Exon` | One annotated exon interval | [Exon reference](../reference/features.md#pyensembl.Exon) |
| `Protein` | An annotated translation | [Protein reference](../reference/features.md#pyensembl.Protein) |

### Gene

`gene.id`, `gene.name`, `gene.contig`, `gene.start`, `gene.end` and `gene.strand`
describe the selected annotation. Name lookup may return multiple loci; stable
IDs avoid choosing a gene by list position. See [aliases](aliases.md) for an
explicit source of additional names.

### Transcript

`gene.transcripts` returns annotated transcripts. `transcript.exons` follows
transcriptional order, including on the minus strand. The transcript's locus
includes introns; its cDNA sequence is spliced and already oriented 5′ to 3′.

### Protein information

`transcript.protein_id` and `transcript.protein_sequence` refer to its annotated
translation. Noncoding or incomplete transcripts may lack a translation or
sequence. Do not treat a missing peptide as evidence that a locus is absent.
Use the [Transcript reference](../reference/features.md#pyensembl.Transcript)
for coding intervals, completeness and sequence behavior.

## Query a genomic interval

This example requires the human release 77 data installed by the following
command. It preserves the original HLA-A lookup:

```sh
pyensembl install --release 77 --species human
```

```python
from pyensembl import EnsemblRelease

data = EnsemblRelease(77, species="human")
gene_names = data.gene_names_at_locus(contig=6, position=29945884)
exon_ids = data.exon_ids_of_gene_name("HLA-A")
print(gene_names)
```

The gene-name result is `["HLA-A"]`. Locus queries use one-based inclusive
coordinates on the selected assembly and can return overlapping genes.
[Lookup methods](../reference/lookups.md) describes filters and nearest-locus queries.
