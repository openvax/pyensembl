# Use dated Ensembl annotations

Ensembl's [new platform](https://www.ensembl.org/) publishes annotations by
genome assembly and date instead of numbered releases. PyEnsembl treats each
annotation date as a release: `EnsemblRelease("2026_04", species="human")`
works like `EnsemblRelease(116, species="human")` and selects GRCh38.
[Choose a reference](guides/assembly-selection.md) compares both with custom
files. Installed datasets are not affected by the transition.

## Install a dated release

This example downloads the GTF, cDNA and peptide files of the human GRCh38
annotation dated 2026_04 from
[its dataset directory](https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/).

```python
from pyensembl import EnsemblRelease

data = EnsemblRelease("2026_04", species="human")
data.download()
data.index()
gene = data.gene_by_id("ENSG00000141510")
print(gene.name, gene.contig, gene.start, gene.end, gene.strand)
```

From the command line:

```sh
pyensembl install --release 2026_04 --species human
```

A date selects the annotation of the species' current assembly, the assembly
of its last numbered release, under PyEnsembl's name for that assembly. Dated
and numbered releases of an assembly share its cache directory, for example
`GRCh38/ensembl2026_04/` beside `GRCh38/ensembl116/`.

For human GRCh38, `2026_04` has the same gene and transcript IDs and versions
as release 116; `2025_12` matches release 115 and `2023_03` matches release 110.
`2024_11` is a GENCODE subset without non-canonical lncRNA transcripts; release
114 has the complete annotation.

Human GRCh38 uses `genes-including_alt.gtf.gz`, the counterpart of the
complete patch and haplotype GTF of numbered GRCh38 releases. Other species
use `genes.gtf.gz`; dated zebrafish GRCz11 therefore omits the alternative
loci of its numbered releases. The cDNA FASTA covers every transcript biotype.
For reference DNA, add `genome_fasta=True`; `genome_fasta_mask="none"`,
`"soft"` or `"hard"` selects the corresponding combined genome FASTA. Only
toplevel DNA is published. [Read genomic DNA](guides/reference-dna.md)
describes coordinates, strand and masking.

## Find a date <a id="choose-the-dataset-explicitly"></a>

`pyensembl available` lists each species' annotation dates beside its numbered
releases. In Python, `available_dated_releases` returns them oldest first:

```python
from pyensembl import available_dated_releases

print(available_dated_releases("human"))
```

```text
['2023_03', '2024_11', '2025_12', '2026_04']
```

The dates come from the directory of the species' current assembly and
provider on the new platform's [FTP site](https://www.ensembl.info/2026/06/26/updates-to-ftp-site-of-the-new-ensembl-website/),
for example [human GRCh38](https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/),
which is where dated releases download from. They are cached on first use and
read offline afterwards; pass `refresh=True` to check for new dates.
`pyensembl available` checks again when its cached dates are a day old, and
uses the cached dates offline. Downloading a date that Ensembl doesn't publish
raises `UnpublishedDateError`, which lists the dates that exist, once the
download finds the file missing.

An annotation date records when the annotation was built, not when it was
published: Arabidopsis's current annotation is dated `2010_09`, and fly's
`2022_07` annotation is older than release 116's. A website release label such
as 2026-07 is different from an annotation date and is rejected.

## Other assemblies and providers

`EnsemblAnnotation` selects any dataset by assembly accession, provider and
date: an older assembly such as GRCh37 (`2013_09`), a species PyEnsembl does
not list, or another provider's annotation.

```python
from pyensembl import EnsemblAnnotation

with EnsemblAnnotation(
    assembly_accession="GCA_000001405.14",
    annotation_date="2013_09",
    provider="ensembl",
    species="human",
    reference_name="GRCh37",
) as data:
    data.download()
    data.index()
```

`available_annotation_dates("GCA_000001405.14", "ensembl")` lists the dates of
an accession and provider, as `available_dated_releases` does for a species.

`include_alt=True` selects `genes-including_alt.gtf.gz`; the default selects
`genes.gtf.gz`. Alternative loci can add name matches. Confirm that the chosen
coverage file exists for your dataset, and use stable gene IDs when a name is
ambiguous. Accession, provider, date and coverage have separate default caches.
If overriding `cache_directory_path`, use a distinct directory for each dataset.
The optional `species` label guards species-specific alias sources; it does not
select or verify the assembly's species. Keep it consistent with the accession.
`genome_fasta=True` and `genome_fasta_mask` select reference DNA as above.

Community providers may have different available files. Check the directory
before installation. Use `Genome` with explicitly matched files when a dataset
does not provide the standard GTF, cDNA and peptide filenames. GFF3 and other
formats require conversion to GTF first. PyEnsembl does not select a rolling
latest annotation or integrate the new GraphQL/refget services.

## Transition from numbered releases

As checked on 2026-10-05, Ensembl's [transition announcement](https://www.ensembl.info/2025/12/02/updates-to-programmatic-access-to-ensembl-and-transitioning-to-the-new-ensembl-platform/)
states that legacy FTP and API services remain available but stop receiving
updates after Ensembl 116; new datasets are delivered through the new platform.
PyEnsembl's numbered-release URLs and selection rules are preserved.
Old species-name directories on the new platform were scheduled for retirement
in August 2026.
