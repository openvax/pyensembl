# Find genes and transcripts <a id="work-with-genes-transcripts-and-proteins"></a>

These examples use the human GRCh38 / Ensembl release 93 data installed on the
[home page](../index.md#install). They run in one Python session.

## Find a gene by name or ID

A name lookup returns a list of every matching gene:

```python
from pyensembl import EnsemblRelease

data = EnsemblRelease(93, species="human")
for gene in data.genes_by_name("TP53"):
    print(gene.id, gene.name, gene.contig, gene.biotype)
```

```text
ENSG00000141510 TP53 17 protein_coding
```

A name can match several loci, such as copies on patch or haplotype contigs.
When you know the gene ID, `data.gene_by_id("ENSG00000141510")` selects it
directly. [Alias lookup](aliases.md) adds names from a separate source, such as
HGNC.

## Find genes at a position <a id="query-a-genomic-interval"></a>

Give a position, or add `end` to search an interval:

```python
print(data.gene_names_at_locus(contig="17", position=7668402))
for gene in data.genes_at_locus(contig="17", position=7660000, end=7690000):
    print(gene.name, gene.strand, gene.start, gene.end)
```

```text
['TP53']
WRAP53 + 7686071 7703502
TP53 - 7661779 7687550
AC087388.1 - 7685260 7686371
```

Coordinates are one-based and inclusive, and must use the selected assembly.
Overlapping genes on both strands are returned; pass `strand="+"` or `"-"` to
restrict the search. `transcripts_at_locus` and `exons_at_locus` work the same
way, and `nearest_gene` finds the closest gene when none overlaps.

## Filter by contig or biotype

Whole-genome lists accept optional `contig`, `strand` and `biotype` filters:

```python
coding_genes = data.genes(contig="17", biotype="protein_coding")
print(len(coding_genes))
```

```text
1183
```

Biotype names come from the annotation and can differ between releases; list
them with `{gene.biotype for gene in data.genes()}`. `gene_ids`, `gene_names`
and `transcript_ids` return identifiers without building objects, and
`data.contigs()` lists the chromosome and contig names.

## Choose among a gene's transcripts <a id="list-a-gene-s-transcripts"></a>

A gene usually has several transcripts. List order does not identify a
preferred isoform, so choose by ID or by attributes that suit your analysis:

```python
gene = data.gene_by_id("ENSG00000141510")
coding = [t for t in gene.transcripts
          if t.biotype == "protein_coding" and t.complete]
print(len(gene.transcripts), len(coding))
coding.sort(key=lambda t: len(t.protein_sequence), reverse=True)
for t in coding[:3]:
    print(t.id, t.name, t.support_level, len(t.protein_sequence))
```

```text
28 19
ENST00000269305 TP53-201 1 393
ENST00000445888 TP53-205 1 393
ENST00000615910 TP53-221 5 382
```

TP53 has 28 transcripts; 19 are protein coding with annotated start and stop
codons. `support_level` is Ensembl's transcript support level, from 1 (best
supported by mRNA evidence) to 5, or `None` when the annotation omits it.
Two transcripts can encode the same protein from different UTRs.

## Exons and sequences <a id="inspect-a-transcript"></a><a id="get-transcript-and-protein-sequences"></a><a id="data-structures"></a><a id="gene"></a><a id="transcript"></a><a id="protein-information"></a>

[Exons and coding sequences](transcripts.md) continues with a transcript's
exons, coding sequence, UTRs, codon positions and protein. The
[method overview](../reference/genome.md#find-a-method) lists every lookup.
Call `data.close()` when you are finished; a
`with EnsemblRelease(93, species="human") as data:` block closes it
automatically.

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
print(data.gene_names_at_locus(contig=6, position=29945884))
```

The result is `["HLA-A"]`. Release 77 also uses GRCh38, but its annotation
differs from release 93, so pin the release you analyze.
