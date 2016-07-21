[![Build Status](https://travis-ci.org/hammerlab/pyensembl.svg?branch=master)](https://travis-ci.org/hammerlab/pyensembl) [![Coverage Status](https://coveralls.io/repos/hammerlab/pyensembl/badge.svg?branch=master&service=github)](https://coveralls.io/github/hammerlab/pyensembl?branch=master) [![DOI](https://zenodo.org/badge/18834/hammerlab/pyensembl.svg)](https://zenodo.org/badge/latestdoi/18834/hammerlab/pyensembl)


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
pyensembl install --release <list of Ensembl release numbers> --species <species-name>
```

For example, `pyensembl install --release 75 76 --species human` will download and install all
human reference data from Ensembl releases 75 and 76.

Alternatively, you can create the `EnsemblRelease` object from inside a Python
process and call `ensembl_object.download()` followed by `ensembl_object.index()`.

## Cache Location
By default, PyEnsembl uses the platform-specific `Cache` folder
and caches the files into the `pyensembl` sub-directory.
You can override this default by setting the environment key `PYENSEMBL_CACHE_DIR`
as your preferred location for caching:

```sh
export PYENSEMBL_CACHE_DIR=/custom/cache/dir
```

or

```python
import os

os.environ['PYENSEMBL_CACHE_DIR'] = '/custom/cache/dir'
# ... PyEnsembl API usage
```

# Non-Ensembl Data

PyEnsembl also allows arbitrary genomes via the specification
of local file paths or remote URLs to both Ensembl and non-Ensembl GTF
and FASTA files. (Warning: GTF formats can vary, and handling of
non-Ensembl data is still very much in development.)

For example:

```
data = Genome
    reference_name='GRCh38',
    annotation_name='my_genome_features',
    gtf_path_or_url='/My/local/gtf/path_to_my_genome_features.gtf'))
# parse GTF and construct database of genomic features
data.index()
gene_names = data.gene_names_at_locus(contig=6, position=29945884)
```

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
: returns list of exons IDs in the annotation database, optionally restricted
by the given chromosome and strand

`exon_ids_of_gene_id(gene_id)`
: returns list of exon IDs associated with a given gene ID

`exon_ids_of_gene_name(gene_name)`
: returns list of exon IDs associated with a given gene name

`exon_ids_of_transcript_id(transcript_id)`
: returns list of exon IDs associated with a given transcript ID

`exon_ids_of_transcript_name(transcript_name)`
: returns list of exon IDs associated with a given transcript name
