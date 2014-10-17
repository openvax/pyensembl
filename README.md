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
: returns all possible gene names

`gene_names_at_locus(contig, position)`
: names of genes overlapping with the given locus

`gene_names_at_loci(contig, start, end)`
: names of genes overlapping with the given range of loci

