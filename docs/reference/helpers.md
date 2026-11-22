# Reference and species helpers

These helpers normalize known species/reference names or select a matching
dataset. For reproducible analysis, prefer an explicit assembly and annotation
version. See [assembly selection](../guides/assembly-selection.md).

::: pyensembl.cached_release

::: pyensembl.available_dated_releases

::: pyensembl.available_annotation_dates

::: pyensembl.UnpublishedDateError

::: pyensembl.find_species_by_name

::: pyensembl.find_species_by_reference

::: pyensembl.genome_for_reference_name

::: pyensembl.which_reference

::: pyensembl.check_species_object

::: pyensembl.normalize_reference_name

::: pyensembl.normalize_species_name

::: pyensembl.species.Species

## Pinned rapid-release helpers

::: pyensembl.pinned_assemblies

::: pyensembl.pinned_reference_aliases

::: pyensembl.resolve_pinned_assembly

::: pyensembl.find_pinned_reference_name

::: pyensembl.is_rcrs_reference

::: pyensembl.contig_aliases

::: pyensembl.normalize_pinned_contig

::: pyensembl.pinned_contig_spellings

::: pyensembl.pinned_contig_lengths

::: pyensembl.mitochondrion_note

## Package values and convenience genomes <a id="exported-constants"></a>

`MAX_ENSEMBL_RELEASE` is the newest numbered release supported by this package.
`RCRS_ACCESSION` (`"NC_012920.1"`) is the RefSeq accession of the revised
Cambridge Reference Sequence used for mitochondrial annotation in pinned
human T2T-CHM13v2.0 releases.
`__version__` is its package version. `ensembl_grch36`, `ensembl_grch37` and
`ensembl_grch38` are convenience `EnsemblRelease` objects selected the first
time each is used, using `genome_for_reference_name`; `import pyensembl` does
not read the cache. Their selections depend on installed/downloaded data at
that moment and the supported release range. Use an explicit
`EnsemblRelease` for reproducibility. They do not represent the new platform's
dated annotations.
