pyensembl
=======

Python interface to Ensembl reference genome metadata (exons, transcripts, &c)

```python
from pyensembl import EnsemblRelease

# release 77 uses human reference genome GRCh38
data = EnsemblRelease(77)

# will return ['HLA-A']
gene_names = data.gene_names_at_locus(contig=6, position=29945884)

# get all exons associated with HLA-A
exon_ids  = data.exon_ids_of_gene_name('HLA-A')
```

# API

The `EnsemblRelease` object has methods to let you access all possible
combinations of the annotation features *gene\_name*, *gene\_id*,
*transcript\_name*, *transcript\_id*, *exon\_id* as well as the location of
these genomic elements (contig, start position, end position, strand).

## Gene Names

`gene_names()`
: returns all gene names in the annotation database

`gene_names_on_contig(contig)`
: all gene names on a particular chromosome/contig

`gene_names_at_locus(contig, position, end=None, strand=None)`
: names of genes overlapping with the given locus
(returns a list to account for overlapping genes)

`gene_name_of_gene_id(gene_id)`
: name of gene with given ID

`gene_name_of_transcript_id(transcript_id)`
: name of gene associated with given transcript ID

`gene_name_of_transcript_name(transcript_name)`
: name of gene associated with given transcript name

`gene_name_of_exon_id(exon_id)`
: name of gene associated with given exon ID


## Gene IDs

`gene_ids(contig=None, strand=None)`
: all gene IDs in the annotation database

`gene_id_of_gene_name(gene_name)`
: translate Ensembl gene ID to its corresponding name


## Transcript Names

`transcript_names(contig=None, strand=None)`
: all transcript names in the annotation database

## Transcript IDs

`transcript_ids(contig=None, strand=None)`
: returns all transcript IDs in the annotation database

`transcript_ids_of_gene_id(gene_id)`
: return IDs of all transcripts associated with given gene ID

`transcript_ids_of_gene_name(gene_name)`
: return IDs of all transcripts associated with given gene name

`transcript_id_of_transcript_name(transcript_name)`
: translate transcript name to its ID

`transcript_ids_of_exon_id(exon_id)`
: return IDs of all transcripts associatd with given exon ID


## Exon IDs

`exon_ids(contig=None, strand=None)`
: returns all transcript IDs in the annotation database

`exon_ids_of_gene_id(gene_id)`

`exon_ids_of_gene_name(gene_name)`

`exon_ids_of_transcript_name(transcript_name)`

`exon_ids_of_transcript_id(transcript_id)`


## Locations

These functions currently assume that each gene maps to a single unique
location, which is invalid both with heavily copied genes
(e.g. [U1](http://en.wikipedia.org/wiki/U1_spliceosomal_RNA)) and with
polymorphic regions (e.g. HLA genes).

`location_of_gene_name(gene_name)`

`location_of_gene_id(gene_id)`

`location_of_transcript_id(transcript_id)`

`location_of_exon_id(exon_id)`


## Start Codons

`start_codon_of_transcript_id`

`start_codon_of_transcript_name`


## Stop Codons

`stop_codon_of_transcript_id`

`stop_codon_of_transcript_name`
