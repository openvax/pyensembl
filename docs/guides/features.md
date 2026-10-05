# Find genes and transcripts <a id="work-with-genes-transcripts-and-proteins"></a>

These examples use the human GRCh38 / Ensembl release 93 data installed on the
[home page](../index.md#install). They cover name, ID and position lookups,
then show how to explore a gene's transcripts.

## Find a gene by name or ID

A name lookup returns every matching gene:

```python
from pyensembl import EnsemblRelease

with EnsemblRelease(93, species="human") as data:
    genes = data.genes_by_name("TP53")
    for gene in genes:
        print(gene.id, gene.name, gene.contig)
```

```text
ENSG00000141510 TP53 17
```

Names can match multiple loci. When you know the intended gene ID, use
`data.gene_by_id("ENSG00000141510")` to select it directly.
[Alias lookup](aliases.md) adds names from a separate source, such as HGNC.

## Query a genomic interval

Find gene names overlapping a position, or add `end` to search an interval:

```python
from pyensembl import EnsemblRelease

with EnsemblRelease(93, species="human") as data:
    print(data.gene_names_at_locus(contig="17", position=7668402))
```

```text
['TP53']
```

Coordinates are one-based and inclusive, and must use the selected assembly.
Use `genes_at_locus(...)` for gene objects instead of names. The optional
`strand` argument restricts the search to `"+"` or `"-"`.

## List a gene's transcripts

`gene.transcripts` returns the transcripts annotated for that gene:

```python
from pyensembl import EnsemblRelease

with EnsemblRelease(93, species="human") as data:
    gene = data.gene_by_id("ENSG00000141510")
    transcripts = gene.transcripts
```

List order does not identify a preferred isoform. Select a transcript by ID or
by a rule appropriate to your analysis.

## Inspect a transcript

After [installing human release 93](../index.md#install), you can inspect a
transcript's gene and genomic span:

```python
from pyensembl import EnsemblRelease

with EnsemblRelease(93, species="human") as data:
    transcript = data.transcript_by_id("ENST00000269305")
    print(transcript.id, transcript.version, transcript.gene_id)
    print(transcript.contig, transcript.start, transcript.end, transcript.strand)
```

```text
ENST00000269305 8 ENSG00000141510
17 7668402 7687538 -
```

The stable transcript ID and its version are separate fields. Its genomic
span includes introns; its spliced cDNA sequence does not.

## Get transcript and protein sequences

Use `transcript.sequence` for spliced cDNA and `transcript.protein_sequence`
for its annotated translation. The [home page](../index.md#read-transcript-and-protein-sequences)
shows an example and output. To include introns or flanking regions instead,
read [genomic DNA](reference-dna.md) from the same assembly.

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

## Examples with other reference data

### Fly annotation <a id="load-genome-in-python"></a>

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

### Human release 77: HLA-A

This example uses human release 77. Install its data before querying HLA-A:

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
