pyensembl
=======

Python interface to Ensembl reference genome metadata (exons, transcripts, &c)

```python
from pyensembl import EnsemblRelease

# release 77 uses human reference genome GRCh38
data = EnsemblRelease(77)

# will return ['HLA-A']
gene_names = data.gene_names_at_locus(contig_name=6, position=29945884)

# get all exons associated with HLA-A
exon_ids  = data.exon_ids_for_gene_name('HLA-A')
```

# API

The `EnsemblRelease` object has methods to let you access all possible
combinations of the annotation features *gene\_name*, *gene\_id*,
*transcript\_name*, *transcript\_id*, *exon\_id* as well as the location of
these genomic elements (contig, start position, end position).

## Gene Names

`gene_names()`
: returns all gene names in the annotation database

`gene_names_at_locus(contig, position)`
: names of genes overlapping with the given locus
(returns a list to account for overlapping genes)

`gene_names_at_loci(contig, start, end)`
: names of genes overlapping with the given range of loci

`gene_name_of_gene_id(gene_id)`
: name of gene with given ID

`gene_name_of_transcript_id(transcript_id)`
: name of gene associated with given transcript ID

`gene_name_of_transcript_name(transcript_name)`
: name of gene associated with given transcript name

`gene_name_of_exon_id(exon_id)`
: name of gene associated with given exon ID

## Gene IDs

`gene_ids()`
: returns all gene IDs in the annotation database

## Transcript Names

`transcript_names()`
: returns all transcript names in the annotation database

## Transcript IDs

`transcript_ids()`
: returns all transcript IDs in the annotation database

## Exon IDs

`exon_ids()`
: returns all transcript IDs in the annotation database


## Locations

These functions currently assume that each gene maps to a single unique
location, which is invalid both with heavily copied genes
(e.g. [U1](http://en.wikipedia.org/wiki/U1_spliceosomal_RNA)) and with
polymorphic regions (e.g. HLA genes).

`location_of_gene_name(gene_name)`

`location_of_gene_id(gene_id)`

`location_of_transcript_id(transcript_id)`

`location_of_exon_id(exon_id)`
