PyEnsembl
=======

Python interface to Ensembl reference genome metadata (exons, transcripts, &c)

# Example Usage

```python
from pyensembl import EnsemblRelease

# release 77 uses human reference genome GRCh38
data = EnsemblRelease(77)

# will return ['HLA-A']
gene_names = data.gene_names_at_locus(contig=6, position=29945884)

# get all exons associated with HLA-A
exon_ids  = data.exon_ids_of_gene_name('HLA-A')
```

# Installation

You can install PyEnsembl using [pip](https://pip.pypa.io/en/latest/quickstart.html):

```sh
pip install pyensembl
```

This should also install any required packages, such as [datacache](https://github.com/hammerlab/datacache) and
[BioPython](http://biopython.org/).

Before using PyEnsembl, run the following command to download and install
Ensembl data:

```
pyensembl install --release <list of Ensembl release numbers>
```

For example, `pyensembl install --release 75 76` will download and install all
data for Ensembl releases 75 and 76.

Alternatively, you can create the `EnsemblRelease` object with
`auto_download=True`. PyEnsembl will then download your data as you
need it, and there will be a delay of several minutes after your first
command.

# API

The `EnsemblRelease` object has methods to let you access all possible
combinations of the annotation features *gene\_name*, *gene\_id*,
*transcript\_name*, *transcript\_id*, *exon\_id* as well as the location of
these genomic elements (contig, start position, end position, strand).

## Genes

`genes(contig=None, strand=None)`
: returns list of Gene objects, optionally restricted to a particular contig
or strand.

`genes_at_locus(contig, position, end=None, strand=None)`
: returns list of Gene objects overlapping a particular position on a contig,
optionally extend into a range with the `end` parameter and restrict to
forward or backward strand by passing `strand='+'` or `strand='-'`.

`gene_by_id(gene_id)`
: return Gene object for given Ensembl gene ID (e.g. "ENSG00000068793")

`gene_names(contig=None, strand=None)`
: returns all gene names in the annotation database, optionally restricted
to a particular contig or strand.

`genes_by_name(gene_name)`
 : get all the unqiue genes with the given name (there might be multiple
due to copies in the genome), return a list containing a Gene object for each
distinct ID.

`gene_by_protein_id(protein_id)`
: find Gene associated with the given Ensembl protein ID (e.g. "ENSP00000350283")

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

`gene_ids(contig=None, strand=None)`
: all gene IDs in the annotation database

`gene_ids_of_gene_name(gene_name)`
: all Ensembl gene IDs with the given name


## Transcripts

`transcripts(contig=None, strand=None)`
: returns list of Transcript objects for all transcript entries in the
Ensembl database, optionally restricted to a particular contig or strand.

`transcript_by_id(transcript_id)`
: construct Transcript object for given Ensembl transcript ID (e.g. "ENST00000369985")

`transcripts_by_name(transcript_name)`
: returns list of Transcript objects for every transcript matching the given name.

`transcript_names(contig=None, strand=None)`
: all transcript names in the annotation database

`transcript_ids(contig=None, strand=None)`
: returns all transcript IDs in the annotation database

`transcript_ids_of_gene_id(gene_id)`
: return IDs of all transcripts associated with given gene ID

`transcript_ids_of_gene_name(gene_name)`
: return IDs of all transcripts associated with given gene name

`transcript_ids_of_transcript_name(transcript_name)`
: find all Ensembl transcript IDs with the given name

`transcript_ids_of_exon_id(exon_id)`
: return IDs of all transcripts associatd with given exon ID


## Exons

`exon_ids(contig=None, strand=None)`
: returns all transcript IDs in the annotation database

`exon_ids_of_gene_id(gene_id)`

`exon_ids_of_gene_name(gene_name)`

`exon_ids_of_transcript_name(transcript_name)`

`exon_ids_of_transcript_id(transcript_id)`

