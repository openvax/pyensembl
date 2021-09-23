<a href="https://app.travis-ci.com/github/openvax/pyensembl">
    <img src="https://app.travis-ci.com/openvax/pyensembl.svg?branch=master" alt="Build Status" />
</a>
<a href="https://coveralls.io/github/openvax/pyensembl?branch=master">
    <img src="https://coveralls.io/repos/openvax/pyensembl/badge.svg?branch=master&service=github" alt="Coverage Status" />
</a>
<a href="https://pypi.python.org/pypi/pyensembl/">
    <img src="https://img.shields.io/pypi/v/pyensembl.svg?maxAge=1000" alt="PyPI" />
</a>


PyEnsembl
=======
PyEnsembl is a Python interface to [Ensembl](http://www.ensembl.org) reference genome metadata such as exons and transcripts. PyEnsembl downloads [GTF](https://en.wikipedia.org/wiki/Gene_transfer_format) and [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files from the [Ensembl FTP server](ftp://ftp.ensembl.org) and loads them into a local database. PyEnsembl can also work with custom reference data specified using user-supplied GTF and FASTA files. 

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

This should also install any required packages, such as [datacache](https://github.com/openvax/datacache) and
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

```python
data = Genome(
    reference_name='GRCh38',
    annotation_name='my_genome_features',
    gtf_path_or_url='/My/local/gtf/path_to_my_genome_features.gtf')
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

<dl>
<dt>genes(contig=None, strand=None)</dt>
<dd>Returns a list of Gene objects, optionally restricted to a particular contig
or strand.</dd>

<dt>genes_at_locus(contig, position, end=None, strand=None)</dt>
<dd>Returns a list of Gene objects overlapping a particular position on a contig,
optionally extend into a range with the end parameter and restrict to
forward or backward strand by passing strand='+' or strand='-'.</dd>

<dt>gene_by_id(gene_id)</dt>
<dd>Return a Gene object for given Ensembl gene ID (e.g. "ENSG00000068793").</dd>

<dt>gene_names(contig=None, strand=None)</dt>
<dd>Returns all gene names in the annotation database, optionally restricted
to a particular contig or strand.</dd>

<dt>genes_by_name(gene_name)</dt>
<dd>Get all the unqiue genes with the given name (there might be multiple
due to copies in the genome), return a list containing a Gene object for each
distinct ID.</dd>

<dt>gene_by_protein_id(protein_id)</dt>
<dd>Find Gene associated with the given Ensembl protein ID (e.g. "ENSP00000350283")</dd>

<dt>gene_names_at_locus(contig, position, end=None, strand=None)
</dt>
<dd>Names of genes overlapping with the given locus, optionally restricted by strand.
(returns a list to account for overlapping genes)</dd>

<dt>gene_name_of_gene_id(gene_id)
</dt>
<dd>Returns name of gene with given genen ID.</dd>

<dt>gene_name_of_transcript_id(transcript_id)
</dt><dd>Returns name of gene associated with given transcript ID.</dd>

<dt>gene_name_of_transcript_name(transcript_name)
</dt>
<dd>Returns name of gene associated with given transcript name.</dd>

<dt>gene_name_of_exon_id(exon_id)
</dt><dd>Returns name of gene associated with given exon ID.</dd>

<dt>gene_ids(contig=None, strand=None)
</dt>
<dd>Return all gene IDs in the annotation database, optionally restricted by
chromosome name or strand.</dd>

<dt>gene_ids_of_gene_name(gene_name)
</dt>
<dd>Returns all Ensembl gene IDs with the given name.</dd>

</dl>

## Transcripts

<dl>
<dt>transcripts(contig=None, strand=None)</dt>
<dd>Returns a list of Transcript objects for all transcript entries in the
Ensembl database, optionally restricted to a particular contig or strand.</dd>

<dt>transcript_by_id(transcript_id)</dt>
<dd>Construct a Transcript object for given Ensembl transcript ID (e.g. "ENST00000369985")</dd>

<dt>transcripts_by_name(transcript_name)</dt>
<dd>Returns a list of Transcript objects for every transcript matching the given name.</dd>

<dt>transcript_names(contig=None, strand=None)</dt>
<dd>Returns all transcript names in the annotation database.</dd>

<dt>transcript_ids(contig=None, strand=None)</dt>
<dd>Returns all transcript IDs in the annotation database.</dd>

<dt>transcript_ids_of_gene_id(gene_id)</dt>
<dd>Return IDs of all transcripts associated with given gene ID.</dd>

<dt>transcript_ids_of_gene_name(gene_name)</dt>
<dd>Return IDs of all transcripts associated with given gene name.</dd>

<dt>transcript_ids_of_transcript_name(transcript_name)</dt>
<dd>Find all Ensembl transcript IDs with the given name.</dd>

<dt>transcript_ids_of_exon_id(exon_id)</dt>
<dd>Return IDs of all transcripts associatd with given exon ID.</dd>
</dl>

## Exons

<dl>
<dt>exon_ids(contig=None, strand=None)</dt>
<dd>Returns a list of exons IDs in the annotation database, optionally restricted
by the given chromosome and strand.</dd>

<dt>exon_by_id(exon_id)</dt>
<dd>Construct an Exon object for given Ensembl exon ID (e.g. "ENSE00001209410")</dd>

<dt>exon_ids_of_gene_id(gene_id)</dt>
<dd>Returns a list of exon IDs associated with a given gene ID.</dd>

<dt>exon_ids_of_gene_name(gene_name)</dt>
<dd>Returns a list of exon IDs associated with a given gene name.</dd>

<dt>exon_ids_of_transcript_id(transcript_id)</dt>
<dd>Returns a list of exon IDs associated with a given transcript ID.</dd>

<dt>exon_ids_of_transcript_name(transcript_name)</dt>
<dd>Returns a list of exon IDs associated with a given transcript name.</dd>
</dl>
