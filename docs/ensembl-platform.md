# The new Ensembl platform

PyEnsembl reads downloadable annotations and sequences into local indexes.
The transition from [Ensembl Beta](https://www.ensembl.org/) to the new website
does not change an installed dataset. Numbered releases remain available through
`EnsemblRelease`; new assembly/date datasets use `EnsemblAnnotation`.

## Choose the dataset explicitly

| Data source | Selection | Python entry point |
| --- | --- | --- |
| [Numbered Ensembl releases](https://ftp.ensembl.org/pub/) | Species and integer release | `EnsemblRelease(93, species="human")` |
| [New Ensembl downloads](https://ftp.ebi.ac.uk/pub/ensemblorganisms/) | Versioned assembly accession, provider and annotation date | `EnsemblAnnotation("GCA_000001405.29", "2023_03")` |
| [Custom GTF and FASTA](https://github.com/openvax/pyensembl#non-ensembl-data) | Matched files selected by the caller | `Genome(...)` |

As checked on 2026-10-05, Ensembl's [transition announcement](https://www.ensembl.info/2025/12/02/updates-to-programmatic-access-to-ensembl-and-transitioning-to-the-new-ensembl-platform/)
states that legacy FTP and API services remain available but stop receiving
updates after Ensembl 116; new datasets are delivered through the new platform.
PyEnsembl's existing numbered-release URLs and selection rules are preserved.

The new platform's [FTP layout](https://www.ensembl.info/2026/06/26/updates-to-ftp-site-of-the-new-ensembl-website/)
groups data by GCA/GCF assembly accession, provider and annotation date.
Old species-name directories were scheduled for retirement in August 2026.
Select an existing directory from the downloads page. A website release label
such as 2026-07 is different from the annotation directory `2023_03`.

## Install a dated annotation

This example selects the documented human GRCh38 assembly accession and dated
Ensembl geneset. It downloads GTF, cDNA and peptide files from
[one dataset directory](https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2023_03/).

```python
from pyensembl import EnsemblAnnotation

with EnsemblAnnotation(
    assembly_accession="GCA_000001405.29",
    annotation_date="2023_03",
    provider="ensembl",
    species="human",
    reference_name="GRCh38",
    include_alt=True,
) as data:
    data.download()
    data.index()
    gene = data.gene_by_id("ENSG00000141510")
    print(gene.name, gene.contig, gene.start, gene.end, gene.strand)
```

`include_alt=True` selects `genes-including_alt.gtf.gz`; the default selects
`genes.gtf.gz`. Alternative loci can add name matches. Confirm that the chosen
coverage file exists for your dataset, and use stable gene IDs when a name is
ambiguous. Accession, provider, date and coverage have separate default caches.
If overriding `cache_directory_path`, use a distinct directory for each dataset.
The optional `species` label guards species-specific alias sources; it does not
select or verify the assembly's species. Keep it consistent with the accession.

For reference DNA, add `genome_fasta=True`. `genome_fasta_mask="none"`, `"soft"`
or `"hard"` selects the corresponding combined genome FASTA from the same
directory. Genome DNA can require several gigabytes. BGZF FASTA files are read
as gzip; reference DNA is materialized and indexed locally for random access.
Coordinates use one-based inclusive intervals; `strand="-"` reverse-complements
DNA, and `mask="raw"` preserves soft masking.

Community providers may have different available files. Check the directory
before installation. Use `Genome` with explicitly matched files when a dataset
does not provide the standard GTF, cDNA and peptide filenames. GFF3 and other
formats require conversion to GTF first. PyEnsembl does not select a rolling
latest annotation or integrate the new GraphQL/refget services.
