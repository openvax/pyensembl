ensembl
=======

Python interface to Ensembl reference genome metadata (exons, transcripts, &c)

```python
from ensembl import EnsemblRelease

# release 77 uses human reference genome GRCh38
data = EnsemblRelease(77)

# will return ['HLA-A']
gene_names = data.gene_names_at_locus(contig_name=6, position=29945884)

# get all exons associated with HLA-A
exon_ids  = data.exon_ids_for_gene_name('HLA-A')
```
